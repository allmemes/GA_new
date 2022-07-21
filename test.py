import check_inbound as ci
import matplotlib.pyplot as plt
import pandas as pd
import source
import numpy as np
from GA_matlab import *
import multiprocessing
import time

R = 8.31446261815324  # universal gas law constant (J K^-1 mol^-1)
T = 273.15 # temp (K)
P = 101325 # standard atmosphere pressure (Pa)
M = 16.04 # molar mass of methane (g/mol)
V = 22.414 # gas molacular volume

def solve_dists_2(df, peak_X, peak_Y):
    '''
    Originally, 
    df: df to record all field data, where "a" refer to the angle that is transferred from field wind angle
        df['a'] = 90 + df[wdir_field]
        df['a'] = np.deg2rad(df['a'])
    However,
    looks like there is no need to transfer field wind angle to the coordinate angle. Now the "a" field is 
    treated as field angle for the demo
    
    peak_X: x coordinate for the source
    peak_Y: y coordinate for the source
    
    Returns    
    -------
    crosswind : Series
        Series of crosswind distances from source to each receptor.
    downwind : Series
        Series of downwind distances from source to each receptor.
    '''
    s = np.sin(df["a"])
    c = np.cos(df["a"])
    x = peak_X - df["x"]
    y = peak_Y - df["y"]
    downWind = (x * s + y * c)
    crossWind = (x * c - y * s) 
    return downWind, crossWind

def solve_distance_Ben(df, peak_X, peak_Y):
    s = np.sin(df["a"])
    c = np.cos(df["a"])
    x = df["x"] - peak_X
    y = df["y"] - peak_Y
    
    crosswind = x * c - y * s
    downwind = x * s + y * c
    
    return crosswind, downwind

def test(source, receptor, input_direction):
    print("----------------------------------------")
    print("The input degree is ", input_direction)
    squared_distance = (receptor[0] - source[0])**2 + (receptor[1] - source[1])**2
    print("The squared distance is ", squared_distance)

    # my method:
    print("My degree is: ", input_direction)
    my_direciton = np.deg2rad(input_direction)
    mySin = np.sin(my_direciton)
    myCos = np.cos(my_direciton)
    my_x = source[0] - receptor[0]
    my_y = source[1] - receptor[1]
    my_downWind = (my_x * mySin + my_y * myCos)
    my_crossWind = (my_x * myCos - my_y * mySin)
    print("My downwind and crosswind are ", my_downWind, my_crossWind)
    print("My calculated distance is ", my_downWind**2 + my_crossWind**2)

    # Ben's method:
    print("Ben's degree is 90 more: ", input_direction + 90)
    Ben_direciton = np.deg2rad(input_direction + 90)
    benSin = np.sin(Ben_direciton)
    benCos = np.cos(Ben_direciton)    
    ben_x = receptor[0] - source[0]
    ben_y = receptor[1] - source[1]
    ben_downWind = ben_x * benSin + ben_y * benCos
    ben_crossWind = ben_x * benCos - ben_y * benSin
    print("Ben's downwind and crosswind are ", ben_downWind, ben_crossWind)
    print("Ben's calculated distance is ", ben_downWind**2 + ben_crossWind**2)
    print("----------------------------------------")



def ppm_to_ugm3(ppm, m, t, p):
    """
    NEED TO TEST/EVALUATE FOR ACCURACY
    
    Convert ppm to micrograms/m^3

    Parameters
    ----------
    ppm : Float or Int
        Parts per million.
    M (or input m) : Float
        Molecular weight of the gas.
    T (or input t) : Float
        Temperature in Kelvin
    P (or input p): Float
        Pressure in Pa.

    R : Constant

    Returns
    -------
    Float
        Micrograms/m^3 of the gas.

    """
    return ((ppm * m * p) / (R * t))

def ppmConversion(ppm, m, t, p):
    return ((ppm * m * 1000) * (T * p)) / (V * P * t)

def ugm3_to_ppm(gm3, m, t, p):
    return ((gm3 * R * t) / (p * m))

def copy(lista):
    newList = []
    for i in lista:
        newList.append(i+1)
    return newList

def calc_time_remaining(itertime_start, iterations, generation, f):
    itertime_end = time.time()
    sec_per_iter = (itertime_end - itertime_start) / 10.0
    time_remaining = round((sec_per_iter * (iterations - generation)) / 60, 2)
    return f"Gen: {generation:<4} Fitness: {f:<20} Est. Time Remaining: {time_remaining}m"

# test
if __name__ == "__main__":
    # source1 = (2,2)
    # receptor1 = (2,-2)
    # test(source1, receptor1, 315)

    # source2 = (1,0)
    # receptor2 = (4,2)
    # test(source2, receptor2, 210)

    # source3 = (2,1)
    # receptor3 = (-2,-2)
    # test(source3, receptor3, 30)

    # source4 = (3,-1)
    # receptor4 = (-2,1)
    # test(source4, receptor4, 150)

    # source5 = (2,2)
    # receptor5 = (0,-4)
    # test(source5, receptor5, 225)
     
    print(ppm_to_ugm3(1, M, T, P))
    print(ppmConversion(1, M ,T, P))
    print(ppm_to_ugm3(0.3, 46, T, P))
    print(ppmConversion(0.3, 46, T, P))
    print(ppm_to_ugm3(1, 64.07, T, P))
    print(ppmConversion(1, 64.07, T, P))

    # listA = [[1,1,1,1,1,1,1,1,2,2,3,3,3,3,3,3,3,3,3,0,0,0],
    #            [1,1,2,2,2,2,2,2,3,3,3,3,3,3,3,3,0,0,0,0,0,0],
    #            [2,2,3,3,3,3,3,3,3,3,3,3,3,3,0,0,0,0,0,0,0,0],
    #            [3,3,3,3,3,3,3,3,3,3,3,3,3,0,0,0,0,0,0,0,0,0],
    #            [3,3,3,3,3,3,3,3,3,3,3,0,0,0,0,0,0,0,0,0,0,0],
    #            [3,3,3,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    #            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
    # dictA = {}
    # mapList = np.array(listA)
    # for i, j in enumerate(mapList):
    #     newList = j.tolist()
    #     for k, l in enumerate(newList):
    #         if l == 3:
    #             newList[k] = window_search(mapList, (i, k))
    #     dictA[i] = newList
    
    # print(dictA)

    # inputPath = r"D:\internship\files\code\GasFlux\result_Raster.json"
    # map_info = read_in_map_info(inputPath)
    # lowerLeftX = map_info[0]
    # lowerLeftY = map_info[1]
    # upperRightX = map_info[2]
    # upperRightY = map_info[3]
    # mapList = map_info[4]
    # cellSize = (upperRightX - lowerLeftX) / len(mapList[0])

    # 1, Test for points inside
    # points_inside = pd.read_csv(r"D:\internship\files\input\station.csv")
    # for i,j in zip(points_inside["X"], points_inside["Y"]):
    #     source_all_in = source.Source(0, i, j, 10)
    #     print(check_in_polygon(mapList, source_all_in, lowerLeftX, upperRightY, cellSize))

 
    scale = 5
    c = 5
    b = 0.01
    iteration = 500
    
    # for i in range(1500):
    #     print(scale)
    #     scale = c - (c/1500)*(i+1)
    
    for i in range(iteration):
       print(scale)
       scale = scale*(1-b*((1+i)/iteration))
    
    # inputPath = r"D:\internship\files\code\GasFlux\result_Raster.json"
    # map_info = ci.read_in_map_info(inputPath)
    # allPoints = ci.generate_all_points_within_polygon(map_info)
    # first_population = create_first_population(30, 10, allPoints)
    # result = Genetic_algorithm(map_info, 1, first_population, 20)
    # print(result)

    # print(multiprocessing.cpu_count())
    
    # a = [1,2,3,4,5,6]
    # for i in range(100):
    #     start = time.time()
    #     time.sleep(2)
    #     print(calc_time_remaining(start, 100, i, i))

