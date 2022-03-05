# CPSC 450 - Bioinformatics
# Isayha Raposo - 230133508
# Assignment 2 - String Reconstruction from Read-Pairs Problem

# I should note that I have created and used setters/getters in this solution despite not being Pythonic
# I have done so (as opposed to directly accessing member variables) primarily for ease of reading/marking

from sys import argv, exit
from os import path

bases = ['A', 'T', 'C', 'G']

# A simple node class that holds:
    # A "paired (pre/suf)fix" value
    # References to nodes previous and nodes next
        # The values of the "edges" connecting said nodes
class Node:
    def __init__(self, paired_fix):
        self.paired_fix = paired_fix
        self.next = dict()
        self.prev = dict()
    
    def edit_next(self, edge, next_node):
        self.next[edge] = next_node

    def edit_prev(self, edge, prev_node):
        self.prev[edge] = prev_node

    def get_next(self):
        return(self.next.items())

    def get_prev(self):
        return(self.prev.items())

    def get_paired_fix(self):
        return(self.paired_fix)

# A class that converts/"glues" specified k-mer pairs into a de Bruijn Graph constructed using Node objects
class de_Bruijn_graph:
    def __init__(self, k_mer_pairs):
        # A master list of the nodes found within the graph
        self.nodes = []

        for k_mer_pair in k_mer_pairs:
            # Creation of paired prefix/suffix nodes
            paired_prefix = [k_mer[:-1] for k_mer in k_mer_pair]
            paired_suffix = [k_mer[1:] for k_mer in k_mer_pair]
            prefix_node = Node(paired_prefix)
            suffix_node = Node(paired_suffix)
            prefix_node.edit_next(k_mer_pair, suffix_node)
            suffix_node.edit_prev(k_mer_pair, prefix_node)
            self.nodes.append(prefix_node)
            self.nodes.append(suffix_node)
        
        # "Glueing" of nodes with identical "paired (pre/suf)fixes" 
        node_hashmap = dict()

        for current_node in self.nodes:
            # Paired fixes are stored in nodes as lists for flexibility
            # Lists are not hashable, thus the conversion below from a list to a simple concatenated string
            current_paired_fix = ''.join(paired_fix for paired_fix in current_node.get_paired_fix())

            # If no other node with the current "paired (pre/suf)fix" has not yet been considered/does not hash to some other (existing) node...
            if current_paired_fix not in node_hashmap:
                # Let the current "paired (pre/suf)fix" map to the current node
                node_hashmap[current_paired_fix] = current_node
            # If the current "paired (pre/suf)fix" has already been considered/hashes to some other (existing) node...
            else:
                existing_node = node_hashmap[current_paired_fix]

                # Label the current and existing nodes as "dominant" or "subservient", where the former has more "next" nodes than the latter
                # This "labelling" is not necessary (either of the nodes could arbitrary be chosen for deletion)
                    # However, theoretically, this should add a minor increase in inefficiency
                if len(existing_node.get_next()) >= len(current_node.get_next()):
                    dominant_node = existing_node
                    subservient_node = current_node
                else:
                    dominant_node = current_node
                    subservient_node = existing_node
                    # If the existing node is labelled as subservient/to be deleted, let the current "paired (pre/suf)fix" now hash to the current node
                    node_hashmap[current_paired_fix] = current_node

                # Rearrange edges in order to "delete" the subservient node (leave it orphaned)
                for edge, next_node in subservient_node.get_next():
                    dominant_node.edit_next(edge, next_node)
                    subservient_node.edit_next(edge, None)
                    next_node.edit_prev(edge, dominant_node)
                for edge, prev_node in subservient_node.get_prev():
                    dominant_node.edit_prev(edge, prev_node)
                    subservient_node.edit_prev(edge, None)
                    prev_node.edit_next(edge, dominant_node)

        # Update the master list of the nodes found within the graph
        self.nodes = list(node_hashmap.values())

    def get_nodes(self):
        return self.nodes

# Returns the (ints) k and d and (list of tuples) paired k-mers found within the user-specified data file
def process_data_file(data_file_name):
    data_file = open(data_file_name, 'r')
    k, d = [int(x) for x in data_file.readline().strip().split(' ')]

    if k < 1 or k > 10:
        print("ERROR: Invalid k value specified. See assignment details.")
        exit(0)
    
    if d < 1 or d > 10:
        print("ERROR: Invalid d value specified. See assignment details.")
        exit(0)

    k_mer_pairs = []

    data_file_line = data_file.readline()
    while(data_file_line):
        k_mers = data_file_line.replace('\n', '').strip().split('|')
        if len(k_mers) != 2:
            print("ERROR: k-mers must be specified in pairs. See assignment details.")
            exit(0)

        for k_mer in k_mers:
            if len(k_mer) != k:
                print("ERROR: Invalid k-mer specified (k-mer length not equal to specified k value).")
                exit(0)

            for base in k_mer.upper():
                if base not in bases:
                    print("ERROR: Invalid k-mer specified (Invalid base/nucleotide specified).")
                    exit(0)

        k_mer_pairs.append(tuple(k_mers))
        data_file_line = data_file.readline()

    return k, d, k_mer_pairs

# Prints string to the console in addition to writing string to output_file
def print_and_write(output_file, string):
    print(string)
    output_file.write(string + '\n')

# A quick means of manually checking the de Bruijn graph constructor (used exclusively to debug alongside proper debug tools)
def check_de_Bruijn_graph(graph):
    print('\n')
    for node in graph.nodes:
        print('Node:', node.get_paired_fix(), 'Next Node(s):', [next_node[1].get_paired_fix() for next_node in node.get_next()], 'Previous Node(s):', [prev_node[1].get_paired_fix() for prev_node in node.get_prev()])
    print('\n')

def euler(graph, current_node, visited, strings=[]):
    visited_copy = visited.copy()
    strings_copy = strings.copy()

    if all(visited_copy.values()):
        return strings_copy

    for next in current_node.get_next():
        next_node = next[1]
        edge = next[0]
        edge_to_take = (current_node, edge, next_node)

        if visited_copy[edge_to_take] == 0:
            visited_copy[edge_to_take] = 1
            strings_copy.append(edge_to_take[1])
            return euler(graph, next_node, visited_copy, strings_copy)

# Processes the solution provided (converts it from an ordered list of k-mer pairs into a single string) and returns None if the solution is found to be incorrect upon verification
def process_solution(raw_solution, k, d):
    first_string = ''
    second_string = ''

    for pair in raw_solution:
        first_k_mer = pair[0]
        second_k_mer = pair[1]

        if first_string == '':
            first_string = first_k_mer
        else:
            if first_string[len(first_string) - k + 1:] == first_k_mer[:k-1]:
                first_string = first_string + first_k_mer[-1]
            else:
                return None
        
        if second_string == '':
            second_string = second_k_mer
        else:
            if second_string[len(second_string) - k + 1:] == second_k_mer[:k-1]:
                second_string = second_string + second_k_mer[-1]
            else:
                return None

    if first_string[k + d:] == second_string[:len(second_string) - k - d]:
        return first_string[:k + d] + second_string
    else:
        return None

def main():
    arg_count = len(argv)
    if arg_count < 2:
        print("ERROR: No command line argument provided. See README.md.")
        exit(1)

    data_file_name = argv[1]
    if not path.isfile(data_file_name):
        print("ERROR: Specified data file " + data_file_name + " not found. See README.md.")
        exit(1)
    
    k, d, k_mer_pairs = process_data_file(data_file_name)

    print("Specified data file found and processed:", "\nk =", k, "\nd =", d, "\nPaired Reads:", k_mer_pairs, '\n')

    #
    # Solution Logic:
    #

    solution = None

    graph = de_Bruijn_graph(k_mer_pairs)
    print("A Paired de Bruijn Graph has been successfully generated from the contents of the specified data file.\n")

    start_node = None

    # Attempt to find a "start" node; one which has no nodes prior...
    for node in graph.get_nodes():
        if len(node.get_prev()) == 0:
            start_node = node

    # If such a "start" node exists, use it as a starting point for the recursive path algorithm called and get the solution
    if start_node is not None:
        # Create the base case for a dict "visited" to ultimately keep track of visited edges recursively (through copying)
        visited = dict()
        # The following 5 lines simply "initialize" the dict above
        for node in graph.get_nodes():
            for next in node.get_next():
                next_node = next[1]
                edge = next[0]
                visited[(node, edge, next_node)] = 0
        
        raw_solution = euler(graph, start_node, visited)
        solution = process_solution(raw_solution, k, d)

    # If such a "start" node DOESN'T exist, try using each node as a starting point for the recursive path algorithm until a verified solution is returned
    else:
        for node in graph.get_nodes():
            start_node = node
            # Create the base case for a dict "visited" to ultimately keep track of visited edges recursively (through copying)
            visited = dict()
            # The following 5 lines simply "initialize" the dict above
            for node in graph.get_nodes():
                for next in node.get_next():
                    next_node = next[1]
                    edge = next[0]
                    visited[(node, edge, next_node)] = 0
            
            potential_raw_solution = euler(graph, start_node, visited)
            potential_solution = process_solution(potential_raw_solution, k, d)
            # potential_solution will be None if potential_raw_solution failed verification. See process_solution().
            if potential_solution is not None:
                solution = potential_solution
                break

    #
    # Output Logic:
    #

    output_file_name = data_file_name.split(".txt")[0] + "_output.txt"
    output_file = open(output_file_name, "w")

    print_and_write(output_file, "Solution: " + solution)
    print("\nThis solution has been written to " + output_file_name)

if __name__ == "__main__":
    main()