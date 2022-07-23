# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 16:28:16 2021

@author: BenJanevic

"""

import math
import pandas as pd
import numpy as np
import time
import random
import matplotlib.pyplot as plt
import csv
import threading
from multiprocessing import Pool
from itertools import repeat, cycle, chain
from source import Source
# from lf_geometry import WasteBoundary as WB
import check_inbound as ci # new check inbound module
import GA_matlab as GA # new matlab GA module
from configparser import ConfigParser, ExtendedInterpolation


# read with ini file
configur = ConfigParser(interpolation=ExtendedInterpolation())
configur.read('metaData.ini')

#############
# CONSTANTS
#############

R = float(configur.get('genetic algorithm', 'R'))
T = float(configur.get('genetic algorithm', 'T'))
P = float(configur.get('genetic algorithm', 'P'))
M = float(configur.get('genetic algorithm', 'M'))
V = float(configur.get('genetic algorithm', 'V'))

# R = 8.31446261815324  # universal gas law constant (J K^-1 mol^-1)
# T = 273.15 # temp (K)
# P = 101325 # standard atmosphere pressure (Pa) 
# M = 16.04 # molar mass of CH4 (g/mol)
# V = 22.414 # molar volume of CH4 @ 0C and 1 atm (L)

################################
# Matlab version of GA VARIABLES
################################

inputMapJson = configur.get('genetic algorithm', 'inputMapJson')
parent_percent = float(configur.get('genetic algorithm', 'parent_percent'))
elite_percent = float(configur.get('genetic algorithm', 'elite_percent'))
crossover_fraction = float(configur.get('genetic algorithm', 'crossover_fraction'))
mutation_scale = float(configur.get('genetic algorithm', 'mutation_scale'))
mutation_shrink = float(configur.get('genetic algorithm', 'mutation_shrink'))

# inputMapJson = "result_Raster.json"
# parent_percent = 0.75
# elite_percent = 0.1
# crossover_fraction = 0.8
# mutation_scale = 6
# mutation_shrink = 0.01


##################
# Old GA VARIABLES
##################

cpu_count = int(configur.get('genetic algorithm', 'cpu_count'))
incsv = configur.get('genetic algorithm', 'incsv')
ch4_field = configur.get('genetic algorithm', 'ch4_field')
wspeed_field = configur.get('genetic algorithm', 'wspeed_field')
wspeed_units = configur.get('genetic algorithm', 'wspeed_units')
wdir_field = configur.get('genetic algorithm', 'wdir_field')
x_field = configur.get('genetic algorithm', 'x_field')
y_field = configur.get('genetic algorithm', 'y_field')
insolation = int(configur.get('genetic algorithm', 'insolation'))
temp_c = float(configur.get('genetic algorithm', 'temp_c'))
p_Pa = float(configur.get('genetic algorithm', 'p_Pa'))
num_solution_sets = int(float(configur.get('genetic algorithm', 'num_solution_sets')))
num_peaks = int(float(configur.get('genetic algorithm', 'num_peaks')))
iterations = int(float(configur.get('genetic algorithm', 'iterations')))

speed = configur.get('genetic algorithm', 'speed')
remove_zero_ppms = configur.getboolean('genetic algorithm', 'remove_zero_ppms')
fixed_sources = configur.getboolean('genetic algorithm', 'fixed_sources')

# cpu_count = 3 # 3 is optimal until multiprocessing is reworked
# incsv = "test_few.csv"
# ch4_field = "ppm"
# wspeed_field = "windspeed"
# wspeed_units = "Meters/Second" #"Miles/Hour" #"Meters/Second"
# wdir_field = "winddir" #"Wind Direction (direction of origin, geographic degrees)"
# x_field = "x"
# y_field = "y"
# insolation = 1
# temp_c = 28.0
# p_Pa = 99322
# num_solution_sets = 50
# num_peaks = 3
# iterations = 500

# speed = "High"
# remove_zero_ppms = True
# fixed_sources = False
# mutation_prob = 0.25
# mut = 2
# eliteSize = 10


# change cwind_pos to dwind_pos
def create_solution_csvs(df, best_fit_set):
    df.drop(columns=["u", "a", "stb1", "stb2", "stb3", "s", "c",
                          "x2", "y2", "crosswind", "downwind", "dwind_pos",
                          "o_y", "o_z"], inplace=True, errors="ignore")
    df.rename(columns={"pred_ch4": "Predicted CH4 (ppm)", "stability": "Stability Class"}, inplace=True)
    df.to_csv("solution.csv")
    
    with open("solution_sources.csv", 'w', newline='') as csvfile: 
        writer = csv.writer(csvfile)
        
        writer.writerow(["ID", x_field, y_field, "Q (g/s)", "Error"])
        for s in best_fit_set[1]["sources"]:
            writer.writerow([s.index, s.x, s.y, round(s.Q / 1000000, 4),
                             best_fit_set[1]["fitness"]])
            

# change cwind_pos to dwind_pos
def calc_solution(df, best_fit_set, temp_k):
    df["pred_ch4"] = 0.0
    for s in best_fit_set[1]["sources"]:
        solve_dists_3(df, s.x, s.y)
        df["dwind_pos"] = df["downwind"].clip(lower=0, upper=1)
        df["dwind_pos"] = np.ceil(df["dwind_pos"])
        briggs_df_2(df)
        calc_concentration_df_2(s.Q, df)
    
    df["Residual"] = df[ch4_field] - df["pred_ch4"]
    mean = df["Residual"].mean()
    std = df["Residual"].std()
    df["Std. Residual"] = (df["Residual"] - mean)/std
    
    # convert back to PPM
    df["pred_ch4"] = df.apply(lambda r: ugm3_to_ppm(r["pred_ch4"], M, temp_k, p_Pa), axis=1)
    df["pred_ch4"] = round(df["pred_ch4"], 2)
    df[ch4_field] = df.apply(lambda r: ugm3_to_ppm(r[ch4_field], M, temp_k, p_Pa), axis=1)
    df[ch4_field] = round(df[ch4_field], 2)


# adjust unit conversion function
def ppm_to_ugm3(ppm, m, t, p):
    """
    NEED TO TEST/EVALUATE FOR ACCURACY
    
    Convert ppm to micrograms/m^3

    Parameters
    ----------
    ppm : Float or Int
        Parts per million.
    M : Float
        Molecular weight of the gas.
    T : Float
        Temperature in Kelvin
    P : Float
        Pressure in Pa.

    Returns
    -------
    Float
        Micrograms/m^3 of the gas.

    """
    return ((ppm * m * p) / (R * t))


# adjust unit conversion function
def ugm3_to_ppm(gm3, m, t, p):
    return ((gm3 * R * t) / (p * m))


def calc_concentration(Q, u, σ_y, σ_z, y):
    """
    Calculates gas concentration in micrograms/meter^3

    Parameters
    ----------
    Q : Float
        Gas constant emission rate in micrograms/second.
    u : Float
        Wind speed in meters/second.
    σ_y : Float
        Horizontal spread parameter (m).
    σ_z : Float
        Vertical spread parameter (m).
    y : Float
        Crosswind distance from ground level to source receptor (m)

    Returns
    -------
    C : Float
        Gas concentration in micrograms/meter^3.

    """
    C = (Q / (math.pi * u * σ_y * σ_z)) * math.exp( -1 * (y**2) / (2 * σ_y**2) )
    return C


# downwind changed to crosswind, binary cwind_pos changed to dwind_pos
def calc_concentration_df(Q, df, pred_ch4, o_y, o_z, crosswind, dwind_pos):
    '''
    Calculate gas concentration in micrograms/m^3. Store result in Pandas Series.

    Parameters
    ----------
    Q : Float
        Gas constant emission rate in micrograms/second.
    df : DataFrame
        DataFrame of all receptors, used for getting windspeed.
    pred_ch4 : Series
        Series of predicted ch4 concentrations associated with each receptor.
    o_y : Series
        Series of horizontal spread parameters associated with each receptor.
    o_z : Series
        Serires of vertical spread parameters asssociated with each receptor.
    crosswind : Series
        Series of crosswind distance from current source to each receptor.
    dwind_pos : Series
        Series of either 1 or 0 for downwind value positive or negative.

    '''
    pred_ch4 += (Q / (math.pi * df["u"] * o_y * o_z) *
        np.exp( -1 * (crosswind**2) / (2 * o_y **2))) * dwind_pos
    

# change downwind to crosswind, change cwind_pos to dwind_pos
def calc_concentration_df_2(Q, df):
    '''
    Calculate gas concentration in micrograms/m^3 and store in dataframe column.

    Parameters
    ----------
    Q : Float
        Gas constant emission rate in micrograms/second.
    df : DataFrame
        DataFrame of all receptors.

    '''
    df["pred_ch4"] += (Q / (math.pi * df["u"] * df["o_y"] * df['o_z']) *
        np.exp(-1 * (df["crosswind"]**2) / (2 * df["o_y"]**2))) * df["dwind_pos"]


def pasquill_class(windspeed, insolation=0):
    windspeed = int(windspeed-2)
    if windspeed < 0:
        windspeed = 0
    if windspeed > 4:
        windspeed = 4
        
    table = [["A", "A", "B"],
             ["A", "B", "C"],
             ["B", "B", "C"],
             ["C", "C", "D"],
             ["C", "D", "D"]]

    return table[windspeed][insolation]


def pasquill_constants(df):
    df.loc[df["stability"] == "A", ["stb1", "stb2", "stb3"]] = [0.22, 0.2, 0.0]
    
    df.loc[df["stability"] == "B", ["stb1", "stb2", "stb3"]] = [0.16, 0.12, 0.0]

    df.loc[df["stability"] == "C", ["stb1", "stb2", "stb3"]] = [0.11, 0.08, 0.0002]

    df.loc[df["stability"] == "D", ["stb1", "stb2", "stb3"]] = [0.08, 0.06, 0.0015]
    
    df.loc[df["stability"] == "E", ["stb1", "stb2", "stb3"]] = [0.06, 0.03, 0.0003]

    df.loc[df["stability"] == "F", ["stb1", "stb2", "stb3"]] = [0.04, 0.016, 0.0003]        
    

# change crosswind to downwind
def briggs_df(df, downwind):
    o_y = df["stb1"] * downwind * (1 + 0.0001 * downwind) ** -0.5
    o_z = df["stb2"] * downwind * (1 + df["stb3"] * downwind) ** -0.5
    return o_y, o_z


# fix downwind and crosswind
def solve_dists_2(df, peak_X, peak_Y):  
    '''
    Solve downwind and crosswind distance from source to receptors.

    Parameters
    ----------
    df : DataFrame
        DESCRIPTION.
    peak_X : Float
        Lon/Easting of source position.
    peak_Y : Float
        Lat/Northing of source position.

    Returns
    -------
    crosswind : Series
        Series of crosswind distances from source to each receptor.
    downwind : Series
        Series of downwind distances from source to each receptor.

    '''
    s = np.sin(df["a"])
    c = np.cos(df["a"])
    x = peak_X - df[x_field]
    y = peak_Y - df[y_field]
    downwind = (x * s + y * c)
    crosswind = (x * c - y * s) 
    
    return crosswind, downwind


# change crosswind to downwind
def briggs_df_2(df):
    df["o_y"] = df["stb1"] * df["downwind"] * (1 + 0.0001 * df["downwind"]) ** -0.5
    df["o_z"] = df["stb2"] * df["downwind"] * (1 + df["stb3"] * df["downwind"]) ** -0.5

# modify downwind and crosswind as solve_dists_2
def solve_dists_3(df, peak_X, peak_Y):
    df["s"] = np.sin(df["a"])
    df["c"] = np.cos(df["a"])
    df["x2"] = peak_X - df[x_field]
    df["y2"] = peak_Y - df[y_field]
    
    df["downwind"] = df["x2"] * df["s"] + df["y2"] * df["c"]
    df["crosswind"] = df["x2"] * df["c"] - df["y2"] * df["s"]

    
def check_coord(x, y, coords):
    return (x, y) in coords


def fitness_2(df, pred_ch4):
    error = (df[ch4_field] - pred_ch4) ** 2
    ch4sq = (df[ch4_field]) ** 2
    error = error.sum() / ch4sq.sum()
    return error


def selection(source_array, eliteSize):
    '''
    Parameters
    ----------
    source_array : TYPE
        List of dicts {source, fitness}.
    eliteSize : TYPE
        DESCRIPTION.
    total_fitness : TYPE
        DESCRIPTION.

    Returns
    -------
    selection_results : TYPE
        DESCRIPTION.

    '''
    selection_results = []
                
    for i in range(eliteSize):
        selection_results.append(source_array[i])
    
    for i in range(len(source_array) - eliteSize):
        sample = random.sample(source_array, 3)
        selection_results.append(min(sample, key=lambda x: x["fitness"]))

    return selection_results


def tournament_selection(source_array, eliteSize, total_fitness):
    pass


def breed(p1, p2, method="swap"):
    child = {"sources": []}
    if method == "splice":
        # List comprehension might create copies, research this
        childp1 = [None for _ in range(len(p1["sources"]))]
        childp2 = [None for _ in range(len(p1["sources"]))]
        
        gene1 = int(random.random() * len(p1["sources"]))
        gene2 = int(random.random() * len(p1["sources"]))
        
        start = min(gene1, gene2)
        end = max(gene1, gene2)
        
        for i in range(start, end):
            childp1[i] = p1["sources"][i]
        
        for i in range(len(p1["sources"])):
            if i < start or i >= end:
                childp2[i] = p2["sources"][i]
                
        for pair in zip(childp1, childp2):
            if pair[0] is not None:
                child["sources"].append(pair[0])
            else:
                child["sources"].append(pair[1])

    if method == "swap":
        for i in range(len(p1["sources"])):
            if random.random() < 0.5:
                child["sources"].append(p1["sources"][i])
            else:
                child["sources"].append(p2["sources"][i])
            
    return child


def crossover(mating_pool, eliteSize, breedmethod):
    children = []
    length = len(mating_pool) - eliteSize
    pool = random.sample(mating_pool, len(mating_pool))
    
    for i in range(eliteSize):
        children.append(mating_pool[i])
        
    for i in range(length):#eliteSize, len(mating_pool)):
        child = breed(pool[i], pool[len(mating_pool)-i-1], breedmethod)
        children.append(child)
        
    return children


def mutate_xy(x, y, wb):
    new_x = x + random.uniform(-100, 100)
    new_y = y + random.uniform(-100, 100)
    
    # NEED TO REMOVE: restricts sources to within landfill bounds
    if new_x < 0:
        new_x = 0
    if new_x > 800:
        new_x = 800
        
    if new_y < 0:
        new_y = 0
    if new_y > 800:
        new_y = 800
        
    return x, y


def mutate_q(q, max_q, wb):
    new_q = q + random.uniform(-5000000, 5000000)
    
    if new_q < 0:
        new_q = 0
        
    if new_q > max_q:
        new_q = max_q

    return new_q


def mutate(population, eliteSize, wb, max_q, mutationRate=0.01):
    mutated_pop = []

    for i, row in enumerate(population):
        new_row = {"sources": []}
        chance1 = random.random()
        chance2 = random.random()
        
        if chance1 < mutationRate and i < eliteSize:
            for s in row["sources"]:
                new_x = s.x
                new_y = s.y
                new_q = s.Q
                
                if chance2 < 0.45 and not s.is_fixed():
                    new_x, new_y = mutate_xy(s.x, s.y, wb)
                    
                if chance2 < 0.9:
                    new_q = mutate_q(s.Q, max_q, wb)
                    
                else:
                    if not s.is_fixed():
                        new_x, new_y = mutate_xy(s.x, s.y, wb)
                    new_q = mutate_q(s.Q, max_q, wb)
                
                new_row["sources"].append(Source(s.index,
                                                  new_x,
                                                  new_y,
                                                  new_q,
                                                  s.fixed))
        else:
            new_row = row
                
        mutated_pop.append(new_row)
        
    return mutated_pop


# change cwind_pos to dwind_pos, flip crosswind and downwind.
def parallel_calc(df, source_array, i, lookup):
    pred_ch4 = pd.Series(np.zeros(len(df.index)))
    for s in source_array[i]["sources"]:
        crosswind, downwind = solve_dists_2(df, s.x, s.y)
        dwind_pos = downwind.clip(lower=0, upper=1)
        dwind_pos = np.ceil(dwind_pos)
        o_y, o_z = briggs_df(df, downwind)
        calc_concentration_df(s.Q, df, pred_ch4, o_y, o_z, crosswind, dwind_pos)
    
    error = fitness_2(df, pred_ch4)
    source_array[i]["fitness"] = error    
    

# change cwind_pos to dwind_pos, flip crosswind and downwind
def parallel_calc_mp(gr, df, lookup):
    '''
    Takes in a list of source configurations calculates the fitness for each
    set of sources.

    Parameters
    ----------
    gr : List
        List of source Dicts.
    df : Dataframe
        Dataframe of receptors.

    Returns
    -------
    result : List
        List of source sets with fitness calculation.

    '''
    result = []
    pred_ch4 = pd.Series(np.zeros(len(df.index)))
    
    for row in gr:
        for s in row["sources"]:
            if not s.is_fixed():
                crosswind, downwind = solve_dists_2(df, s.x, s.y)
                dwind_pos = downwind.clip(lower=0, upper=1)
                dwind_pos = np.ceil(dwind_pos)
                o_y, o_z = briggs_df(df, downwind)
                calc_concentration_df(s.Q, df, pred_ch4, o_y, o_z, crosswind, dwind_pos)
            
            else:
                calc_concentration_df(s.Q, df, pred_ch4, lookup[s.index]["o_y"],
                                      lookup[s.index]["o_z"], lookup[s.index]["crosswind"],
                                      lookup[s.index]["dwind_pos"])
        
        error = fitness_2(df, pred_ch4)
        row["fitness"] = error
        result.append(row)
        
        pred_ch4.values[:] = 0.0
        
    return result


def cycle_baskets(items, maxbaskets):
    baskets = [[] for _ in range(min(maxbaskets, len(items)))]
    for item, basket in zip(items, cycle(baskets)):
        basket.append(item)
    return baskets


def calc_time_remaining(itertime_start, iterations, generation, f): 
    itertime_end = time.time()
    sec_per_iter = (itertime_end - itertime_start) / 10.0
    time_remaining = round((sec_per_iter * (iterations - generation)) / 60, 2)
    return f"Gen: {generation:<4} Fitness: {f:<20} Est. Time Remaining: {time_remaining}m"


def cull_and_randomize_sources(source_array, peaks_dict, wb):
    # UPDATE TO HANDLE FIXED SOURCES
    i = len(source_array) // 2
    source_array = source_array[:i]
    
    for i in range(num_solution_sets - i):
        row = {"sources": [], "fitness": 0.0}
        for k, v in peaks_dict.items():
            x, y = wb.rand_point() 
            row["sources"].append(Source(k, x, y, random.random() * 5000000))
        source_array.append(row)


# change cwind_pos to dwind_pos, flip crosswind and downwind
def precompute_source_parameters(source_array, df):
    '''
    Maps precomputed parameters for eaech receptor to each source.

    Parameters
    ----------
    source_array : TYPE
        DESCRIPTION.
    df : TYPE
        DESCRIPTION.

    Returns
    -------
    lookup : TYPE
        DESCRIPTION.
    '''
    lookup = {}
    
    for s in source_array[0]["sources"]:
        crosswind, downwind = solve_dists_2(df, s.x, s.y)
        dwind_pos = downwind.clip(lower=0, upper=1)
        dwind_pos = np.ceil(dwind_pos)
        o_y, o_z = briggs_df(df, downwind)
        
        lookup[s.index] = {"o_y": o_y, "o_z": o_z, "crosswind": crosswind,
                           "dwind_pos": dwind_pos}
        
    return lookup
                

def main():
    # 0, read in map json to make the first population and provide other check inbound functions.
    map_info = ci.read_in_map_info(inputMapJson)
    mut_scale = mutation_scale
    # get all points within polygon
    all_points_inbound = ci.generate_all_points_within_polygon(map_info)

    # 1, read in receptor csv, convert temperature unit
    temp_k = temp_c + 273.15
    df = pd.read_csv(incsv)
    
    # 2, convert concentration unit.
    df[ch4_field] = df.apply(lambda r: ppm_to_ugm3(r[ch4_field], M, temp_k, p_Pa), axis=1)
    # max_ugm3 = df[ch4_field].max()
    # calculate max emission rate based on max concentration measured
    # NOT YET IMPLEMENTED
    
    # 3, create waste boundary based on receptor points
    # Need to fix
    # wb = WB(list(zip(df[x_field], df[y_field])))

    # 4, calculate peak numbers.
    df['peak'] = df[ch4_field][(df[ch4_field].shift(1) < df[ch4_field]) & (df[ch4_field].shift(-1) < df[ch4_field])]
    peaks = df[df["peak"] > 0]
    peaks = peaks.nlargest(num_peaks, ch4_field)
    
    # 5, generate sources with random locations within hull and random emissions
    #   from 0 to the max emission rate
    # MAX EMISSION RATE NOT YET IMPLEMENTED
    source_array = []
    peaks_dict = peaks.to_dict('index')
    
    for i in range(num_solution_sets):
        row = {"sources": [], "fitness": 0.0}
        for k, v in peaks_dict.items():
            if fixed_sources:
                x, y = v[x_field], v[y_field]
            else:
                # x, y = wb.rand_point()
                random_source = random.choice(all_points_inbound)
                x = random_source[0]
                y = random_source[1]

            row["sources"].append(Source(k, x,
                                         y,
                                         random.random() * 25000000,
                                         fixed_sources))
        source_array.append(row)
        
    # 6, for test use
    # coords = set()
    # for v in peaks_dict.values():
    #     coords.add((v[x_field], v[y_field]))

    # 7, filter out 0 ppm if needed.
    if remove_zero_ppms:
        df = df[df[ch4_field] != 0.0]
    
    # 8, Drop where windspeed is 0, model unsuitable for these cases
    df = df[df[wspeed_field] != 0.0]
    if wspeed_units != "Meters/Second":
        df['u'] = df[wspeed_field] / 2.237
    else:
        df['u'] = df[wspeed_field]
    
    # 9, set the initial predicted concentration for all receptors to 0
    df["pred_ch4"] = 0.0

    # 10, convert the unit of wind direction
    # df['a'] = 90 + df[wdir_field]
    # df['a'] = np.deg2rad(df['a'])
    df['a'] = np.deg2rad(df[wdir_field])

    # 11, calculate atmospheric stability
    df["stability"] = df.apply(lambda r: pasquill_class(r['u'], insolation), axis=1)
    pasquill_constants(df)

    lookup = None
    if fixed_sources:
        lookup = precompute_source_parameters(source_array, df)

    # 12, criteria and variables for the genetic iteration.
    best_fit = math.inf
    best_fit_set = [0, None]
    fitness_tracker, worst_tracker = [], []
    start = time.time()

    # 13, genetic algorithm
    for generation in range(iterations):
        try:
            if generation % 10 == 0:
                itertime_start = time.time()

            if speed == "Medium":
                threads = []
                for i in range(len(source_array)):
                    process = threading.Thread(target=parallel_calc,
                                               args=[df, source_array, i, lookup])
                    process.start()
                    threads.append(process)

                for process in threads:
                    process.join()

                source_array.sort(key=lambda x: x["fitness"])

            elif speed == "High" :
                grouped_sources = cycle_baskets(source_array, cpu_count)
                with Pool(processes=cpu_count) as p:
                    result = p.starmap(parallel_calc_mp, zip(grouped_sources, repeat(df), repeat(lookup)))

                source_array = list(chain(*result))
                source_array.sort(key=lambda x: x["fitness"])
                

            f = source_array[0]["fitness"]
            fitness_tracker.append(f)
            worst_f = source_array[-1]["fitness"]
            worst_tracker.append(worst_f)

            print(np.mean([i["fitness"] for i in source_array]))

            if generation < 10:
                print(f"Gen: {generation:<4} Fitness: {f:<20}")

            if generation % 10 == 0 and generation != 0:
                print(calc_time_remaining(itertime_start, iterations, generation, f))
                #cull_and_randomize_sources(source_array, peaks_dict, wb)

            if f < best_fit:
                best_fit = f
                best_fit_set[1] = source_array[0]
                best_fit_set[0] = generation

            # mating_pool = selection(source_array, eliteSize)
            # children = crossover(mating_pool, eliteSize, breedmethod="swap")
            # source_array = mutate(children, eliteSize, wb, max_q=23000000, mutationRate=0.01)

            # New GA:
            # 1, create parent set.
            parent_size = int(parent_percent * num_solution_sets)
            parent_set = GA.select_parents_by_fitness(source_array, parent_size)

            # 2, first offspring component: select elite children from population directly
            elite_size = int(elite_percent * num_solution_sets)
            elite_child = source_array[0:elite_size]

            # 3, second offspring component: select crossover children from parent set
            crossover_size = int(crossover_fraction * (num_solution_sets - elite_size))
            crossover_child = GA.swap_crossover(parent_set, crossover_size)

            # 4, third offspring component: select mutation children from parent set
            mutation_size = num_solution_sets - elite_size - crossover_size
            mutation_child = GA.mutation(parent_set, mut_scale, mutation_size, map_info)
            # matlab version of mutation scale decrease:
            mut_scale = mut_scale * (1 - mutation_shrink*((i+1)/iterations))
            
            # 5, replace the old population with new generation.
            source_array = elite_child + crossover_child + mutation_child

        except KeyboardInterrupt:
            break

    stop = time.time()

    # 14, related output
    print(f"\nProcessing time: {round(stop-start)}s ")
    if generation < iterations:
        print(f"Generation reached: {generation}\n")
        
    print("Calculating solution...")
    calc_solution(df, best_fit_set, temp_k)
        
    print("Creating CSV products...")
    create_solution_csvs(df, best_fit_set)
    
    print("Plotting results...")
    plt.plot(fitness_tracker)
    plt.show()
    plt.plot(worst_tracker)
    plt.show()
    df.plot.scatter(x=x_field, y=y_field, c="Predicted CH4 (ppm)", cmap="viridis")

            
if __name__ == "__main__":
    main()