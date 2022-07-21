import numpy as np
import random
import source
import json
import pandas as pd
import time

# raster method
def read_in_map_info(inputFilePath):
    mapList = []
    with open(inputFilePath, "r") as file:
        data = json.load(file)
        lowerLeftX = data["extent"]["llx"]
        lowerLeftY = data["extent"]["lly"]
        upperRightX = data["extent"]["urx"]
        upperRightY = data["extent"]["ury"]
        for i in data.keys():
            if i.isnumeric():
                mapList.append(data[i])
    return lowerLeftX, lowerLeftY, upperRightX, upperRightY, mapList

def check_in_polygon(mapList, point:source.Source, lowerLeftX, upperRightY, cell_size):
    # convert to coordinates by treating the lower left corner as the origin.
    x = int((point.x - lowerLeftX) / (cell_size))
    y = int((upperRightY - point.y) / (cell_size))
    if x < 0:
        return False
    if y < 0:
        return False
    try:
        if mapList[y][x] >= 1:
            return True
        else:
            return False
    except:
        return False

def move_into_boundary(mapList, point:source.Source, lowerLeftX, upperRightY, cell_size):
    x = int((point.x - lowerLeftX) / (cell_size))
    y = int((upperRightY - point.y) / (cell_size))
    if x < 0:
        print("Coordiniate outside boudary")
        return None
    if y < 0:
        print("Coordiniate outside boudary")
        return None
    try:
        if isinstance(mapList[y][x], list):
            return mapList[y][x]
        else:
            print("Points already inside boundary")
            return None
    except:
        print("Coordiniate outside boudary")
        return None

def generate_all_points_within_polygon(map_info):
    lowerLeftX = map_info[0]
    upperRightX = map_info[2]
    upperRightY = map_info[3]
    mapList = map_info[4]
    cell_size = (upperRightX - lowerLeftX) / len(mapList[0])

    all_points_within = []
    for i, j in enumerate(mapList):
        for k, l in enumerate(j):
            try:
                if l >= 1:
                    projected_x = lowerLeftX + cell_size * k + 1/2 * cell_size
                    projected_y = upperRightY - cell_size * i - 1/2 * cell_size
                    all_points_within.append((projected_x, projected_y))
            except:
                continue
    return all_points_within

if __name__ == "__main__":
    '''
    Some variables required to run the functions
    '''
    inputPath = r"D:\internship\files\code\GasFlux\result_Raster.json"
    map_info = read_in_map_info(inputPath)
    lowerLeftX = map_info[0]
    lowerLeftY = map_info[1]
    upperRightX = map_info[2]
    upperRightY = map_info[3]
    mapList = map_info[4]
    cellSize = (upperRightX - lowerLeftX) / len(mapList[0])


    '''
    Uncomment them to run test.
    '''

    # 1, Test for points inside
    # points_inside = pd.read_csv(r"D:\internship\files\input\station.csv")
    # for i,j in zip(points_inside["X"], points_inside["Y"]):
    #     source_all_in = source.Source(1, i, j, 10)
    #     print(check_in_polygon(mapList, source_all_in, lowerLeftX, upperRightY, cellSize))

    # # 2, Test for points outside
    # points_outside = pd.read_csv(r"D:\internship\files\input\points_outside.csv")
    # for i,j in zip(points_outside["x"], points_outside["y"]):
    #     source_all_out = source.Source(1, i, j, 10)
    #     print(check_in_polygon(mapList, source_all_out, lowerLeftX, upperRightY, cellSize))

    # # 3, Test for points near boundary.
    # bounudary_source1 = source.Source(1, 805518.95, 4705340.93, 10)
    # print(check_in_polygon(mapList, bounudary_source1, lowerLeftX, upperRightY, cellSize))
    # bounudary_source2 = source.Source(1, 807660.58, 4694950.56, 10)
    # print(check_in_polygon(mapList, bounudary_source2, lowerLeftX, upperRightY, cellSize))
    # slightly_out = source.Source(1, 807639.41, 4694926.85, 10)
    # print(check_in_polygon(mapList, slightly_out, lowerLeftX, upperRightY, cellSize))

    # 4, Test for 1000000 random points within polygon
    # test_points = generate_all_points_within_polygon(map_info)
    # print(len(test_points))
    # start = time.time()
    # for i in range(1000000):
    #     data = random.choice(test_points)
    #     test_source = source.Source(1, data[0], data[1], 10)
    #     check_in_polygon(mapList, test_source, lowerLeftX, upperRightY, cellSize)
    # end = time.time()
    # print(end - start)

    # # Another test to output random points to csv for visualization.
    # dicta = {"x":[], "y":[]}
    # test_points = generate_all_points_within_polygon(map_info)
    # for i in range(100):
    #     data = random.choice(test_points)
    #     dicta["x"].append(data[0])
    #     dicta["y"].append(data[1])
    # outputDF = pd.DataFrame(data = dicta)
    # outputDF.to_csv(r"D:\internship\files\input\random_inside_points.csv")

    # # 5, Check whether there is 0 in new mapList
    # print(cellSize)
    # for i in mapList:
    #     for j in i:
    #         try:
    #             if j == 0:
    #                 print("found 0")
    #         except:
    #             continue

    # # 6, Test for moving points inbound
    # dicta = {"x":[], "y":[]}
    # nearby_points = pd.read_csv(r"D:\internship\files\input\nearby_points.csv")
    # for i, j in zip(nearby_points["x"], nearby_points["y"]):
    #     points = source.Source(1, i, j, 10)
    #     coordinate = move_into_boundary(mapList, points, lowerLeftX, upperRightY, cellSize)
    #     try:
    #         dicta["x"].append(coordinate[0])
    #         dicta["y"].append(coordinate[1])
    #     except:
    #         continue
    # outputDF = pd.DataFrame(data = dicta)
    # outputDF.to_csv(r"D:\internship\files\input\moved_points.csv")