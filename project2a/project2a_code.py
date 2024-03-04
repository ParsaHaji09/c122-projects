##########################################################################################
# GOAL: recreate original reference genome from spectrum of 480 reads using eulerian cycle
##########################################################################################
""" functions """
def overlap_graph(patterns):
    result = {}
    for pattern1 in patterns:
        if pattern1 in result:
            continue
        result[pattern1] = []
        for pattern2 in patterns:
            if pattern1[1:] == pattern2[:-1] and pattern2 not in result[pattern1]:
                result[pattern1].append(pattern2)
        if not result[pattern1]:
            del result[pattern1]
    return result

def find_start_end(adj_list):
    count_in = {}

    start = ''
    end = ''

    for value in adj_list.values():
        for i in value:
            if i in count_in:
                count_in[i] += 1
            count_in[i] = 1
    
    potential_list = []
    for key in count_in.keys():
        if key not in adj_list:
            end = key
        elif count_in[key] != len(adj_list[key]):
            potential_list.append(key)
    
    for key in adj_list.keys():
        if key not in count_in:
            start = key
        elif count_in[key] != len(adj_list[key]):
            potential_list.append(key)
    
    return start, end

def eulerian_cycle(adj, start):
    edge_count = dict()
    for key in adj.keys():
        edge_count[key] = len(adj[key])
    curr_path = []
    circuit = []
    curr_path.append(start)
    curr_v = start
    while len(curr_path):
        if edge_count[curr_v]:
            curr_path.append(curr_v)
            next_v = adj[curr_v][-1]
            edge_count[curr_v] -= 1
            adj[curr_v].pop()
            curr_v = next_v
        else:
            circuit.append(curr_v)
            curr_v = curr_path[-1]
            curr_path.pop()
    result = []
    for i in range(len(circuit) - 1, -1, -1):
        result.append(circuit[i])
    return result
 
##########################################################################################
""" main """

spectrum_reads = []
file_path = './project2a_spectrum.fasta'
with open(file_path, 'r') as f:
    for line in f:
        line = line.strip()
        if line and line[0] != '>':
            spectrum_reads.append(line)

adjacency_matrix = overlap_graph(spectrum_reads)
start,end = find_start_end(adjacency_matrix)

adjacency_matrix[end] = [start]

result = eulerian_cycle(adjacency_matrix, start)

indices = []
for i in result:
    indices.append(spectrum_reads.index(i))

indices.pop()

with open('result.txt', 'a') as file:
    for i in indices:
        file.write(">read_" + str(i) + f"\n")