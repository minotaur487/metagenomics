def generate_sorted_aligned_reads(alignment_list, path):
    f = open(path, 'w')
    write_alignment_list(f, alignment_list)
    f.close()
    return

def format_aligned_reads(alignment_map, text_to_id_map):
    alignment_list = []
    read_id_ordering = []
    for read, subreads_list in alignment_map.items():
        alignment_list.append((read, subreads_list[0]))
    alignment_list.sort(key=lambda x: x[1][0])
    for entry in alignment_list:
        read = entry[0]
        read_id = text_to_id_map[read]
        read_id_ordering.append(read_id)
    return read_id_ordering

def write_alignment_list(f, alignment_list):
    for read_id in alignment_list:
        s = '>read_' + read_id + '\n'
        f.write(s)
    return
