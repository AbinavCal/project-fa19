# from __future__ import print_function
import os
import sys
sys.path.append('..')
sys.path.append('../..')
import argparse
import utils
from student_utils import *
# import numpy


import random
from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp

# from pyscipopt import Model, quicksum, multidict
# from cvxopt import matrix, solvers

"""
======================================================================
  Complete the following function.
======================================================================
"""

def solve(list_of_locations, list_of_homes, starting_car_location, adjacency_matrix, params=[]):
    """
    Write your algorithm here.
    Input:
        list_of_locations: A list of locations such that node i of the graph corresponds to name at index i of the list
        list_of_homes: A list of homes
        starting_car_location: The name of the starting location for the car
        adjacency_matrix: The adjacency matrix from the input file
    Output:
        A list of locations representing the car path
        A dictionary mapping drop-off location to a list of homes of TAs that got off at that particular location
        NOTE: both outputs should be in terms of indices not the names of the locations themselves
    """
    # SHORTHANDS:
    # nx = networkx
    locations = list_of_locations
    homes = list_of_homes

    # print(locations)
    # print(homes)

    # Find index i where the car starts
    start_index = locations.index(starting_car_location)

    # Convert adjacency matrix representation to graph in networkx [Return Type is (graph, msg=_)]
    graph, _ = adjacency_matrix_to_graph(adjacency_matrix)
    
    # Store number of "vertices" V
    V = len(locations)

    # Run networkx's Floyd Warshall Alg (V^3 runtime, V^2 space)
        #   distance = {(source, target), dist} dictionary of shortest path distance
        #   pred     = {(source, target), ?} dictionary of predecessors
    distance = nx.floyd_warshall(graph)

    # print(distance[0][1])
    # Calculate an approximate optimal D = len(dropoffs)
    D = random.randint(1, V)

    # Find dropoff points D
    # dropoff_points = range(D)       # TODO

    x = set()
    
    while len(x) < (int) (len(locations)/5):
        y = random.randint(1, len(locations))
        if y not in x:
            x.add(y)
    dropoff_points = list(x)
    
    # dropoff_points = []    
    # while len(dropoff_points) < len(homes):
    #     dropoff_points.append(random.randint(1, len(locations)))

    # print("dropoff: ", dropoff_points)
    # print("\n")
    
    # print("adjacency matrix:", adjacency_matrix)

    # def flp(I,J,d,M,f,c):
    #     model = Model("flp")
    #     x,y = {},{}
    #     for j in J:
    #         y[j] = model.addVar(vtype="B", name="y(%s)"%j)
    #         for i in I:
    #             x[i,j] = model.addVar(vtype="C", name="x(%s,%s)"%(i,j))
    #     for i in I:
    #         model.addCons(quicksum(x[i,j] for j in J) == d[i], "Demand(%s)"%i)
    #     for j in M:
    #         model.addCons(quicksum(x[i,j] for i in I) <= M[j]*y[j], "Capacity(%s)"%i)
    #     for (i,j) in x:
    #         model.addCons(x[i,j] <= d[i]*y[j], "Strong(%s,%s)"%(i,j))
    #     model.setObjective(
    #         quicksum(f[j]*y[j] for j in J) +
    #         quicksum(c[i,j]*x[i,j] for i in I for j in J),
    #         "minimize")
    #     model.data = x,y
    #     return model

    # # c =     

    # model = flp(I, J, d, M, f, c)
    # model.optimize()
    # EPS = 1.e-6
    # x,y = model.__data
    # edges = [(i,j) for (i,j) in x if model.GetVal(x[i,j]) > EPS]
    # facilities = [j for j in y if model.GetVal(y[j]) > EPS]
    # print ("Optimal value=", model.GetObjVal())
    # print ("Facilities at nodes:", facilities)
    # print ("Edges:", edges)



    # A = matrix([ [-1.0, -1.0, 0.0, 1.0], [1.0, -1.0, -1.0, -2.0] ])
    # b = matrix([ 1.0, -2.0, 0.0, 4.0 ])
    # c = matrix([ 2.0, 1.0 ])
    # sol=solvers.lp(c,A,b)
    # print(sol['x'])
    # dropoff_points = sol['x']


    # Compute the TSP on dropoffs
    induced = graph.subgraph(dropoff_points)
    induced_adj_matrix = nx.to_numpy_matrix(induced)
    # print("graph", graph)
    # print("induced", induced)
    # print("adj matrix", induced_adj_matrix)
    # print("\n")
    data = {}
    data['distance_matrix'] = induced_adj_matrix
    data['num_vehicles'] = 1
    data['depot'] = 0

    # Helper fn 1: returns total path (with last point = start point)
    def get_path(manager, routing, assignment):
        index = routing.Start(0)
        path = []
        while not routing.IsEnd(index):
            path.append(manager.IndexToNode(index))
            index = assignment.Value(routing.NextVar(index))
        path.append(path[0])
        return path

    # Helper fn 2: used by OR-Tools to find distances between nodes
    def distance_callback(from_index, to_index):
        from_node = manager.IndexToNode(from_index)
        to_node = manager.IndexToNode(to_index)
        return data['distance_matrix'][from_node, to_node]

    # Create OR-Tools' routing index manager and Routing Model
    manager = pywrapcp.RoutingIndexManager(len(data['distance_matrix']), data['num_vehicles'], data['depot'])
    routing = pywrapcp.RoutingModel(manager)

    transit_callback_index = routing.RegisterTransitCallback(distance_callback)

    # Define cost of each arc.
    routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

    # Setting first solution heuristic.
    search_parameters = pywrapcp.DefaultRoutingSearchParameters()
    search_parameters.first_solution_strategy = (routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)

    # Solve TSP and determine path/total distance
    assignment = routing.SolveWithParameters(search_parameters)
    if assignment:
        car_path = get_path(manager, routing, assignment)

    # Create dropoff location dictionary
    home_indices = convert_locations_to_indices(homes, locations)
    dropoff_dictionary = {}
    for index in car_path:
        if index in dropoff_points:                                                               # TODO
            x = random.choice(home_indices)
            home_indices.remove(x)
            dropoff_dictionary[index] = [x]                                       # TODO

    # DEBUG COST
    c, m = cost_of_solution(graph, car_path, dropoff_dictionary)
    # print("cost:", c)
    # print(m)
    # print("car_path:", car_path)
    # print("\n")

    # Return two dictionaries
    # print(dropoff_dictionary)
    return car_path, dropoff_dictionary    

"""
======================================================================
   No need to change any code below this line
======================================================================
"""

"""
Convert solution with path and dropoff_mapping in terms of indices
and write solution output in terms of names to path_to_file + file_number + '.out'
"""
def convertToFile(path, dropoff_mapping, path_to_file, list_locs):
    string = ''
    for node in path:
        string += list_locs[node] + ' '
    string = string.strip()
    string += '\n'

    dropoffNumber = len(dropoff_mapping.keys())
    string += str(dropoffNumber) + '\n'
    for dropoff in dropoff_mapping.keys():
        strDrop = list_locs[dropoff] + ' '
        for node in dropoff_mapping[dropoff]:
            strDrop += list_locs[node] + ' '
        strDrop = strDrop.strip()
        strDrop += '\n'
        string += strDrop
    utils.write_to_file(path_to_file, string)

def solve_from_file(input_file, output_directory, params=[]):
    print('Processing', input_file)

    input_data = utils.read_file(input_file)
    num_of_locations, num_houses, list_locations, list_houses, starting_car_location, adjacency_matrix = data_parser(input_data)
    car_path, drop_offs = solve(list_locations, list_houses, starting_car_location, adjacency_matrix, params=params)

    basename, filename = os.path.split(input_file)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    output_file = utils.input_to_output(input_file, output_directory)

    convertToFile(car_path, drop_offs, output_file, list_locations)


def solve_all(input_directory, output_directory, params=[]):
    input_files = utils.get_files_with_extension(input_directory, 'in')

    for input_file in input_files:
        solve_from_file(input_file, output_directory, params=params)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Parsing arguments')
    parser.add_argument('--all', action='store_true', help='If specified, the solver is run on all files in the input directory. Else, it is run on just the given input file')
    parser.add_argument('input', type=str, help='The path to the input file or directory')
    parser.add_argument('output_directory', type=str, nargs='?', default='.', help='The path to the directory where the output should be written')
    parser.add_argument('params', nargs=argparse.REMAINDER, help='Extra arguments passed in')
    args = parser.parse_args()
    output_directory = args.output_directory
    if not args.all:
        input_directory = args.input
        solve_all(input_directory, output_directory, params=args.params)
    else:
        input_file = args.input
        solve_from_file(input_file, output_directory, params=args.params)
