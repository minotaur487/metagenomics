import collections
import random

import alignment
import read_util

GENOME_PATH_SUFFIX = '.fasta'
GENOME_PATH_PREFIX = './data/genome_'

def find_match_frequency(read_map, read_length):
    match_frequency_map = collections.defaultdict(int)
    i = 0
    while True:
        try:
            match_frequency_map = count_matches(i, read_map, match_frequency_map, read_length)
        except FileNotFoundError:
            print(f'Last valid genome id: {i-1}')
            break
        i += 1
    return match_frequency_map

def count_matches(i, read_map, match_frequency_map, read_length):
    genome_id, genome_length, sequence = get_genome_info(i)
    _, num_of_aligned_reads = alignment.align_reads(read_map, sequence, read_length)
    match_frequency_map[genome_id] = num_of_aligned_reads
    return match_frequency_map

def find_most_probable_genomes(match_frequency_map, threshold):
    keys = list(match_frequency_map.keys())
    for genome_ids in keys:
        if match_frequency_map[genome_ids] < threshold:
            del match_frequency_map[genome_ids]
    return list(match_frequency_map.keys())

def map_reads_to_genomes(most_probable_genome_ids, id_to_read_map, read_to_id_map, read_length):
    read_to_genome_id_map = {}
    for genome_id in most_probable_genome_ids:
        genome_id, genome_length, sequence = get_genome_info(genome_id)
        aligned_reads, _ = alignment.align_reads(id_to_read_map, sequence, read_length)
        for read in aligned_reads:
            read_id = read_to_id_map[read]
            if read_id in read_to_genome_id_map:
                print('duplicate mapping')
            read_to_genome_id_map[read_id] = genome_id
    
    for i in range(len(id_to_read_map)):
        if str(i) not in read_to_genome_id_map:
            read_to_genome_id_map[str(i)] = random.choice(most_probable_genome_ids)
        
    return read_to_genome_id_map

def get_genome_info(i):
    genome_path = GENOME_PATH_PREFIX + str(i) + GENOME_PATH_SUFFIX
    f = open(genome_path, 'r')
    genome_id, genome_length, sequence = read_util.get_genome(f)
    return genome_id, genome_length, sequence
