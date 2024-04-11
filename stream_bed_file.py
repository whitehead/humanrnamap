import csv


def filter_bedgraph_file(filename, start, end):
    filtered_data = []
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            pos = int(row[1])
            if start <= pos <= end:
                filtered_data.append({'pos': pos, 'D1': float(row[2])})
    return filtered_data

# Usage
chrom = 'chrX'  # Replace with the desired chromosome
f_s = 12345     # Replace with the start position
f_e = 23456     # Replace with the end position

filename = f'D1_{chrom}_mmRate_pos.bedGraph'
gene_D1 = filter_bedgraph_file(filename, f_s, f_e)
