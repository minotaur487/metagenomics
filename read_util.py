def get_genome(f):
    identifier_line = f.readline().strip()
    genome_id, genome_length = format_genome_identifier_line(identifier_line)

    sequence = ''
    while True:
        line = f.readline().strip()
        if line == '':
            break
        sequence += line
    return genome_id, genome_length, sequence

def format_genome_identifier_line(line: str):
    line = line.split()
    genome_id = line[0][15:]
    _, length = line[2].split(':')
    length = int(length)
    return genome_id, length

def format_paired_reads(f):
    def get_paired_read(f):
        identifier_line = f.readline().strip()
        if identifier_line == '':
            return '', '', ''
        read = f.readline().strip()

        _, ids = identifier_line.split('_')
        read_id, pos_in_pair = ids.split('/')

        return read_id, pos_in_pair, read
    
    paired_reads = {}
    length_of_read = None
    i = 0
    while True:
        read_id, pos_in_pair, read = get_paired_read(f)
        if read_id == '':
            break
        length_of_read = len(read)
        paired_reads[(read_id, pos_in_pair)] = read
        i += 1
    return paired_reads, length_of_read

def format_single_reads(f):
    def get_single_read(f):
        identifier_line = f.readline().strip()
        if identifier_line == '':
            return '', ''
        read = f.readline().strip()
        read_id = identifier_line[6:]
        return read_id, read
    
    id_to_read_map = {}
    read_to_id_map = {}
    read_length = 0
    while True:
        read_id, read = get_single_read(f)
        if read_id == '':
            break
        read_length = max(len(read), read_length)
        id_to_read_map[read_id] = read
        read_to_id_map[read] = read_id
    
    return id_to_read_map, read_to_id_map, read_length

class AlignedRead:
    def __init__(self, pos=-1, hamming_distance=-1):
        self.pos = pos
        self.hamming_distance = hamming_distance
    
    def __iter__(self):
        return iter((self.pos, self.hamming_distance))

    def __str__(self):
        return f'pos: {self.pos}, h distance: {self.hamming_distance}'
    
    def __repr__(self):
        return f'pos: {self.pos}, h distance: {self.hamming_distance}'
