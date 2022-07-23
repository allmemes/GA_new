# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 13:47:26 2021

@author: BenJanevic
"""
from logging import config
import random
import math
import csv
import math
import pandas as pd
from matplotlib import pyplot
import check_inbound as ci
from configparser import ConfigParser, ExtendedInterpolation


def calc_concentration(Q, u, σ_y, σ_z, y):
    C = (Q / (math.pi * u * σ_y * σ_z)) * math.exp( -1 * (y**2) / (2 * σ_y**2) )
    return C

def ppm_to_ugm3(ppm, m, t, p):
    """
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

def ugm3_to_ppm(gm3, m, t, p):
    return ((R * t * gm3) / (m * p))

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

def briggs(stability, x):
    if stability == "A":
        return (0.22 * x * (1 + 0.0001 * x) ** -0.5,
                0.2 * x)
    
    if stability == "B":
        return (0.16 * x * (1 + 0.0001 * x) ** -0.5,
                0.12 * x)
    
    if stability == "C":
        return (0.11 * x * (1 + 0.0001 * x) ** -0.5,
                0.08 * x * (1 + 0.0002 * x) ** -0.5)
    
    if stability == "D":
        return (0.08 * x * (1 + 0.0001 * x) ** -0.5,
                0.06 * x * (1 + 0.0015 * x) ** -0.5)
    
    if stability == "E":
        return (0.06 * x * (1 + 0.0001 * x) ** -0.5,
                0.03 * x * (1 + 0.0003 * x) ** -0.5)
    
    if stability == "F":
        return (0.04 * x * (1 + 0.0001 * x) ** -0.5,
                0.016 * x * (1 + 0.0003 * x) ** -0.5)

def rotate(pivot_x, pivot_y, point_x, point_y, angle):
    # s = math.sin(angle)
    # c = math.cos(angle)
    # point_x -= pivot_x
    # point_y -= pivot_y
    # x = point_x * c - point_y * s
    # y = point_x * s + point_y * c
    # return x, y  # cross, along
    s = math.sin(angle)
    c = math.cos(angle)
    x = pivot_x - point_x
    y = pivot_y - point_y
    downwind = (x * s + y * c)
    crosswind = (x * c - y * s)

    return crosswind, downwind


class Point:
    def __init__(self, x, y, ppm=0.0, is_source=False, q=0.0, wind_dir=0.0,
                 wind_speed=0.0):
        self.x = x
        self.y = y
        self.ppm = ppm
        self.is_source = is_source
        self.q = q
        self.downwind = 0
        self.crosswind = 0
        self.wind_dir = wind_dir
        self.wind_speed = wind_speed


    def __repr__(self):
        return f"({self.x}, {self.y}): {self.ppm}"
    
    
    def __str__(self):
        return self.__repr__

# R = 8.31446261815324  # universal gas law constant (J K^-1 mol^-1)
# M = 16.04
# size = 40
# gridsize = 20
# max_q = 25
# base_wind_dir = 45
# base_wind_speed = 2
# jitter = 1


configur = ConfigParser(interpolation=ExtendedInterpolation())
configur.read('metaData.ini')
R = float(configur.get('genetic algorithm', 'R'))
M = float(configur.get('genetic algorithm', 'M'))
precision = int(configur.get('test', 'precision'))
num_sources = int(configur.get('genetic algorithm', 'num_peaks'))
insolation = int(configur.get('genetic algorithm', 'insolation'))
temp_k = 273.15 + float(configur.get('genetic algorithm', 'temp_c'))
p_Pa = float(configur.get('genetic algorithm', 'p_Pa'))
outReceptors = configur.get('test', 'outReceptors')
outSource = configur.get('test', 'outSource')
inputCsv = configur.get('test', 'inputCsv')
mapJson = configur.get('genetic algorithm', 'inputMapJson')


mapInfo = ci.read_in_map_info(mapJson)
allPointsWithin = ci.generate_all_points_within_polygon(mapInfo)
receptors = pd.read_csv(inputCsv)

# random.seed(1)
m = []
for i, j, k, l in zip(receptors["X"], receptors["Y"], receptors["Direction"], receptors["Speed"]):
    m.append(Point(round(i,precision), round(j,precision), wind_dir=k, wind_speed=l))

sources = []
for i in range(num_sources):
    source = random.choice(allPointsWithin)
    x = round(source[0], precision)
    y = round(source[1], precision)
    # x = round(random.uniform(0, gridsize * size), 2)
    # y = round(random.uniform(0, gridsize * size), 2)
    q = round(random.random() * 25000000, precision)
    sources.append((x, y, q))

# for row in m:
#     for r in row:
#         for s in sources:
#             Q = s[2]
#             #x, y = solve_dists(s[0], s[1], r.x, r.y, wind_dir) 52508.27770620771
#             cross, along = rotate(s[0], s[1], r.x, r.y, math.radians(r.wind_dir))
            
#             cross = round(cross, precision)
#             along = round(along, precision)
#             if along > 0:
#                 σ_y, σ_z = briggs(pasquill_class(r.wind_speed, 2), along)
#                 try:
#                     C = round(calc_concentration(Q, r.wind_speed, σ_y, σ_z, cross), precision)
#                     C = ugm3_to_ppm(C, M, temp_k, p_Pa)
#                 except ZeroDivisionError:
#                     C = 0
#                 r.ppm += C
#             r.downwind = along
#             r.crosswind = cross

for i in m:
    for s in sources:
        Q = s[2]
        #x, y = solve_dists(s[0], s[1], r.x, r.y, wind_dir) 52508.27770620771
        cross, along = rotate(s[0], s[1], i.x, i.y, math.radians(i.wind_dir))

        cross = round(cross, precision)
        along = round(along, precision)
        if along > 0:
            σ_y, σ_z = briggs(pasquill_class(i.wind_speed, 2), along)
            try:
                C = calc_concentration(Q, i.wind_speed, σ_y, σ_z, cross)
                C = ugm3_to_ppm(C, M, temp_k, p_Pa)
            except ZeroDivisionError:
                C = 0
            i.ppm += C
        i.downwind = along
        i.crosswind = cross

xs = []
ys = []
ppms = []

with open(outReceptors, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    
    csvwriter.writerow(["x", "y", "ppm", "windspeed", "winddir", "is_source", "Q", "Alongwind", "Crosswind"])
    
    # for row in m:
    #     for p in row:
    #         csvwriter.writerow([p.x, p.y, p.ppm, p.wind_speed, p.wind_dir, p.is_source, p.q, p.downwind, p.crosswind])
    #         xs.append(p.x)
    #         ys.append(p.y)
    #         ppm = p.ppm
    #         if ppm > 100:
    #             ppm = 100
    #         ppms.append(ppm)

    for p in m:
        csvwriter.writerow([p.x, p.y, round(p.ppm, precision), p.wind_speed, p.wind_dir, p.is_source, p.q, p.downwind, p.crosswind])
        xs.append(p.x)
        ys.append(p.y)
        ppm = p.ppm
        if ppm > 100:
            ppm = 100
        ppms.append(ppm)

with open(outSource, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["x", "y", "Q"])
    for row in sources:
        csvwriter.writerow(row)

pyplot.scatter(xs, ys, c=ppms, cmap='viridis')
pyplot.colorbar()
pyplot.show()