import collections
from read_util import AlignedRead
from utils import build_genome_index_map, find_hamming_distance

SIZE_OF_PARTITION = 10

def align_reads(read_map, reference_genome, length_of_read):
    def get_possible_read_positions(genome_index_map):
        possible_alignment_map = collections.defaultdict(list)
        for read_id, read in read_map.items():
            possible_positions = []
            for i in range(0, len(read), SIZE_OF_PARTITION):
                subread = read[i:i+SIZE_OF_PARTITION]
                if subread in genome_index_map:
                    possible_positions.append(genome_index_map[subread])
                else:
                    possible_positions.append([None])
            possible_alignment_map[read] = possible_positions
        return possible_alignment_map

    def get_number_of_partitions(length_of_read):
        number_of_partitions = length_of_read / SIZE_OF_PARTITION
        return int(number_of_partitions)
    
    def align_reads_with_buckets(possible_read_positions, num_of_partitions):
        def cast_subread_position_vote(beginning, end, read, majority_map):
            key = (beginning, end)
            key_if_deletion = (beginning - 1, end - 1)  # negating shift from deletion
            key_if_insertion = (beginning + 1, end + 1) # negating shift from insertion

            if read not in majority_map:
                majority_map[read] = { key : 1 }
            else:
                if key_if_deletion in majority_map[read]:
                    majority_map[read][key_if_deletion] += 1
                elif key in majority_map[read]:
                    majority_map[read][key] += 1
                elif key_if_insertion in majority_map[read]:
                    majority_map[read][key_if_insertion] += 1
                elif key not in majority_map[read]:
                    majority_map[read][key] = 1
            return majority_map
        
        def add_subread_to_bucket(beginning, end, subread_beginning, subread_end, read, buckets):
            def handle_missing_first_pos(buckets, key, read):
                subread_start, _ = buckets[read][key][0]
                first_pos = key[0]
                if subread_start != first_pos:
                    buckets[read][key].insert(0, (first_pos, first_pos+SIZE_OF_PARTITION))
                return buckets
            
            key = (beginning, end)
            key_if_deletion = (beginning - 1, end - 1)  # negating shift from deletion
            key_if_insertion = (beginning + 1, end + 1) # negating shift from insertion
            
            actual_key = key
            subread_start_end = (subread_beginning, subread_end)
            if read not in buckets:
                buckets[read] = { key : [subread_start_end] }
            else:
                if key_if_deletion in buckets[read]:
                    actual_key = key_if_deletion
                elif key_if_insertion in buckets[read]:
                    actual_key = key_if_insertion

                if actual_key in buckets[read]:
                    buckets[read][actual_key].append(subread_start_end)
                elif actual_key not in buckets[read]:
                    buckets[read][actual_key] = [subread_start_end]
            buckets = handle_missing_first_pos(
                buckets,
                actual_key,
                read
            )
            return buckets

        def calculate_votes_and_generate_subread_position_list(i, read, subread_positions_lists, majority_map, buckets):
            for candidate_position in subread_positions_lists:
                if candidate_position is None:
                    continue
                
                read_beginning = candidate_position - (i * SIZE_OF_PARTITION)
                read_end = candidate_position + ((num_of_partitions - i) * SIZE_OF_PARTITION)
                subread_beginning = candidate_position
                subread_end = candidate_position + SIZE_OF_PARTITION

                majority_map = cast_subread_position_vote(read_beginning, read_end, read, majority_map)
                buckets = add_subread_to_bucket(
                    read_beginning,
                    read_end,
                    subread_beginning,
                    subread_end,
                    read,
                    buckets
                )

            return buckets, majority_map

        def generate_alignment_buckets(possible_read_positions):
            potential_buckets = {}    # read : list of lists of start and end of subreads
            majority_map = {} # read : list of votes per position
            for read, possible_positions_list in possible_read_positions.items():
                for i, subread_positions_lists in enumerate(possible_positions_list):
                    potential_buckets, majority_map = calculate_votes_and_generate_subread_position_list(
                        i,
                        read,
                        subread_positions_lists,
                        majority_map,
                        potential_buckets
                    )
            most_probable_aligned_reads = choose_most_likely_positions(potential_buckets, majority_map)
            return most_probable_aligned_reads
        
        def choose_most_likely_positions(all_buckets, majority_map):
            most_probable_aligned_reads = {}
            for read, votes in majority_map.items():
                most_probable_range = max(votes.keys(), key=votes.get)
                most_probable_aligned_reads[read] = all_buckets[read][most_probable_range]
            return most_probable_aligned_reads
        most_probable_aligned_reads = generate_alignment_buckets(possible_read_positions)
        return most_probable_aligned_reads

    def format_aligned_reads(aligned_reads_range_map, reference_genome):
        alignment_map = collections.defaultdict(list)
        for read, (beginning, end) in aligned_reads_range_map.items():
            window = reference_genome[beginning:end]
            hamming_distance = find_hamming_distance(read, window)
            alignment_map[read] = [AlignedRead(beginning, hamming_distance)]
        return alignment_map

    num_of_partitions = get_number_of_partitions(length_of_read)
    if num_of_partitions is None:
        raise ValueError('Partition size doesn\'t evenly divide into length of read.')

    genome_index_map = build_genome_index_map(reference_genome, SIZE_OF_PARTITION)
    possible_read_positions = get_possible_read_positions(genome_index_map)

    aligned_reads = align_reads_with_buckets(possible_read_positions, num_of_partitions)
    num_of_aligned_reads = len(aligned_reads)
    return aligned_reads, num_of_aligned_reads
