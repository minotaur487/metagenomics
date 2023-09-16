import collections

def find_hamming_distance(s1, s2):
    mismatches = 0
    s1_length = len(s1)
    s2_length = len(s2)
    length = min(s1_length, s2_length)
    for i in range(length):
        if s1[i] != s2[i]:
            mismatches += 1
    return mismatches

def build_genome_index_map(reference_genome, length):
    genome_index_map = collections.defaultdict(list)
    for i in range(len(reference_genome) - length):
        genome_index_map[reference_genome[i:i+length]].append(i)
    return genome_index_map

def generate_output(reads_to_genome_map, output_path):
    f = open(output_path, 'w')
    read_id_list = [int(x) for x in list(reads_to_genome_map.keys())]
    sorted_read_ids = sorted(read_id_list)
    for read_id in sorted_read_ids:
        read_id = str(read_id)
        genome_id = reads_to_genome_map[read_id]
        line = '>read_' + str(read_id) + '\tGenome_Number_' + str(genome_id) + '\n'
        f.write(line)
    return
