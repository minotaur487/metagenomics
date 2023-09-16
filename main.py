import match
import read_util
import subprocess
import utils

THRESHOLD = 3000
OUTPUT_PATH = './predictions.csv'
READS_PATH = './data/reads.fasta'

def main():
    f_reads = open(READS_PATH, 'r')
    id_to_read_map, read_to_id_map, read_length = read_util.format_single_reads(f_reads)
    match_frequency_map = match.find_match_frequency(id_to_read_map, read_length)
    most_probable_genome_ids = match.find_most_probable_genomes(match_frequency_map, THRESHOLD)
    reads_to_genome_map = match.map_reads_to_genomes(most_probable_genome_ids, id_to_read_map, read_to_id_map, read_length)
    utils.generate_output(reads_to_genome_map, OUTPUT_PATH)
    subprocess.run(['zip', '-r', 'predictions.zip', 'predictions.csv'])

main()
