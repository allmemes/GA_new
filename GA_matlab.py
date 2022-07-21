import numpy as np
import random
import source
import check_inbound as ci

def select_parents_by_fitness(source_array, parentSize):
    fitness_sum = 0
    for i in source_array:
        fitness_sum += i["fitness"]

    probability = []
    for i in source_array:
        probability.append(1 / (i["fitness"] / fitness_sum))

    selected_individuals = random.choices(source_array, weights=probability, k=parentSize)
    return selected_individuals

def select_parents_by_tournament():
    pass 

def comparator(source):
    return source["fitness"]

def select_elite(source_array, eliteSize):
    sortedList = sorted(source_array, key = comparator)
    return sortedList[0:eliteSize]

def single_crossover(parent, outputSize):
    pass 

def two_crossover(parent, outputSize):
    pass 

def swap_crossover(parent, outputSize):
    children = []
    for i in range(outputSize):
        father = random.choice(parent)
        mother = random.choice(parent)
        child = {"sources":[], "fitness": 0}
        for j in range(len(father["sources"])):
            flag = random.randint(0, 1)
            if flag:
                child["sources"].append(father["sources"][j])
            else:
                child["sources"].append(mother["sources"][j])
        children.append(child)
    return children

def mutation(parent, scale, outputSize, mapInfo):
    # unpack the map information for move in bound function
    lowerLeftX = mapInfo[0]
    upperRightX = mapInfo[2]
    upperRightY = mapInfo[3]
    mapList = mapInfo[4]
    cellSize = (upperRightX - lowerLeftX) / len(mapList[0])

    children = []
    for i in range(outputSize):
        single_parent = random.choice(parent)
        child = {"sources":[], "fitness":0}
        for j in single_parent["sources"]:
            newX = j.x
            newY = j.y

            # if the sources are not fixed, then mutate the x and y. if they are fixed, only mutate the concentration.
            if not j.is_fixed():
                delta_x = np.random.normal(0, scale)
                if abs(delta_x) > 4*scale:
                    if delta_x > 0:
                        delta_x = 4*scale
                    else:
                        delta_x = -4*scale
                newX = newX + delta_x

                delta_y = np.random.normal(0, scale)
                if abs(delta_y) > 4*scale:
                    if delta_y > 0:
                        delta_y = 4*scale
                    else:
                        delta_y = -4*scale
                newY = newY + delta_y
                
            newQ = j.Q + np.random.normal(0, scale*100000)
            

            # check whether inside the polygon, move if necessary
            newSource = source.Source(j.index, newX, newY, newQ, j.fixed)
            if not ci.check_in_polygon(mapList, newSource, lowerLeftX, upperRightY, cellSize):
                movedCoord = ci.move_into_boundary(mapList, newSource, lowerLeftX, upperRightY, cellSize)
                newSource.x = movedCoord[0]
                newSource.y = movedCoord[1]

            child["sources"].append(newSource)
        children.append(child)
    return children

def Genetic_algorithm(mapInfo, iteration, population, parent_size,
                        elite_percent = 0.05, parent_method = "fit_scaled",
                        crossover_fraction = 0.8, crossover_method = "swap",
                        mutation_scale = 1, mutation_shrink = 1):
    population_size = len(population)
    for i in range(1, iteration + 1):
        # calcuate predicted concentration.

        # evaluate fitness

        # start breeding
        # first component: elite children
        elite_size = int(elite_percent * population_size)
        elite_child = select_elite(population, elite_size)

        # select parent set for cross over and mutation
        if parent_method == "fit_scaled":
            parents = select_parents_by_fitness(population, parent_size)
        elif parent_method == "tournament":
            parents = select_parents_by_tournament()
        else:
            print("parent selection method not supported")
            return None

        # determine the size of crossover children
        crossover_size = int(crossover_fraction * (population_size - elite_size))
        # second component: crossover child
        if crossover_method == "swap":
            crossover_child = swap_crossover(parents, crossover_size)
        elif crossover_method == "single":
            crossover_child = single_crossover(parents, crossover_size)
        elif crossover_method == "two":
            crossover_child = two_crossover(parents, crossover_size)
        else:
            print("corossover method not supported")
            return None

        # determine the size of mutation children
        mutation_size = population_size - elite_size - crossover_size
        # third component: mutation child
        mutation_child = mutation(parents, mutation_scale, mutation_size, mapInfo)
        mutation_scale = mutation_scale * (1 - mutation_shrink*(i/iteration))

        # add three components up
        population = elite_child + crossover_child + mutation_child
    
    return population

# create test population
def create_first_population(population_size, source_num, all_points_inbound):
    source_list = []
    for i in range(population_size):
        source_config = {"sources":[], "fitness":0}
        source_config["fitness"] = round(random.uniform(0, 3),3)
        for j in range(source_num):
            random_source = random.choice(all_points_inbound)
            x = random_source[0]
            y = random_source[1]
            concentration = random.randrange(0, 100)

            source_point = source.Source(j, x, y, concentration)
            source_config["sources"].append(source_point)
        source_list.append(source_config)
    return source_list

if __name__ == "__main__":
    # random.seed(10)
    # inputPath = r"D:\internship\files\code\GasFlux\result_Raster.json"
    # inputPath = "/Users/apple/Desktop/python/internship/GasFlux/GasFlux/result_raster.json"
    # map_info = ci.read_in_map_info(inputPath)
    # allPoints = ci.generate_all_points_within_polygon(map_info)
    # first_population = create_first_population(30, 10, allPoints)
    # result = Genetic_algorithm(map_info, 1, first_population, 20)
    # print(result)
    pass