# Copyright (C) 2023 Samson Weiner (samson.weiner@uconn.edu) and
# Mukul S. Bansal (mukul.bansal@uconn.edu).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/.

import os
import subprocess
import argparse
import time
import sys

import pandas as pd
import numpy as np

np.set_printoptions(threshold=sys.maxsize)

# Read data into pandas dataframe
def process_data(data_path, use_total):
    df = pd.read_csv(data_path, delimiter='\t')
    df = df.set_index(['CELL', 'chrom', 'start'])
    df = df.drop(columns=['end'])
    df = df.sort_index()
    if use_total:
        df.columns = [0]
    else:
        df = df['CN states'].str.split(',', expand=True)
    df = df.astype('int32')
    return df

# Manhattan distance between rows a and b.
def distance_manhattan(a, b):
    return np.sum(np.abs(np.subtract(a, b)))

# Euclidean distance between rows a and b.
def distance_euclidean(a, b):
    return np.dot(a,a) - 2*np.dot(a,b) + np.dot(b,b)

# Root distance between rows a and b.
def distance_root(a, b):
    return np.sum(np.sqrt(np.abs(np.subtract(a,b))))

# Log distance between rows a and b
def distance_log(a, b):
    return np.sum(np.log(np.abs(np.subtract(a, b)) + 1))

# Computes cell by cell distance matrix.
def compute_distmatrix(prof_df, use_breakpoint, metric, use_total):
    cell_names = prof_df.index.get_level_values(0).drop_duplicates().to_numpy()
    chrom_names = prof_df.index.get_level_values(1).drop_duplicates().to_numpy()
    n = cell_names.shape[0]

    dist_matrix = np.zeros((n, n))

    for chrom in chrom_names:
        if use_total:
            alleles = [0]
        else:
            alleles = [0, 1]
        for allele in alleles:
            cur_prof = prof_df.xs(chrom, level='chrom')[allele]
            cur_prof = cur_prof.unstack(level='start').to_numpy()
            if use_breakpoint:
                cur_prof = np.apply_along_axis(np.ediff1d, 1, cur_prof)

            for i in range(cell_names.size - 1):
                for j in range(i+1, cell_names.size):
                    dist_matrix[i,j] += metric(cur_prof[i], cur_prof[j])
    if metric == distance_euclidean:
        dist_matrix = np.sqrt(dist_matrix)

    dist_matrix = dist_matrix + dist_matrix.T - np.diag(np.diag(dist_matrix))
    return dist_matrix, cell_names

# Saves the distance matrix in PHYLIP format
def save_distmatrix(dist_matrix, cell_names, file_name):
    n = len(cell_names)
    max_char = max([len(cell) for cell in cell_names]) + 1
    with open(file_name, 'w+') as f:
        f.write(str(n) + '\n')
        for i in range(n):
            cell = cell_names[i]
            x = np.array2string(dist_matrix[i], formatter={'float_kind':lambda x: "%.5f" % x})[1:-1].replace('\n', '')
            sep_char = ' ' * (max_char-len(cell))
            f.write(cell + sep_char + x + '\n')

# Build the tree using the balanced ME of fastme
def construct_tree_fast(file_name, save_dm, tree_prefix, method, use_NNI, fastme_path, seed):
    if method == 'fastNJ':
        flags = ['N']
    elif method == 'fastuNJ':
        flags = ['U']
    elif method == 'balME':
        if use_NNI:
            flags = ['B', '-n']
        else:
            flags = ['B', '-s']
    elif method == 'olsME':
        if use_NNI:
            flags = ['O', '-n']
        else:
            flags = ['O', '-s']
    if seed:
        flags += ['-z', str(seed)]
    call = subprocess.run([fastme_path, '-i', file_name, '-o', tree_prefix, '-m'] + flags, capture_output=True, text=True)
    if not save_dm:
        os.remove(file_name)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to dataset.')
    parser.add_argument('-o', '--output', type=str, default=None, help='Output directory. Will create one if it does not exist. Default is current directory.')
    parser.add_argument('-p', '--prefix', type=str, default=None, help='Prefix to add to output files.')
    parser.add_argument('-s', '--save-dm', action='store_true', help='Toggle to save the distance matrix to a file.')
    parser.add_argument('-b', '--breakpoint', action='store_true', help='Toggle to use breakpoint profiles.')
    parser.add_argument('-t', '--total-cn', action='store_true', help='Use total copy numbers instead of allele-specific copy numbers.')
    parser.add_argument('-d', '--dist-type', type=str, default='root', help='Distance measure type.')
    parser.add_argument('-m', '--rec-method', type=str, default=None, help='Phylogenetic reconstruction algorithm. If not specified, will compute the distance matrix and save to a file. Options are \'fastNJ\', \'fastuNJ\', \'balME\', \'olsME\'.')
    parser.add_argument('-n', '--use-NNI', action='store_true', help='For ME methods, toggle to use NNI tree search. By default, SPR tree search is used.')
    parser.add_argument('-f', '--fastme-path', type=str, default='fastme', help='Path to \'fastme\' executable. By default, assumes the fastme executable is added to the user $PATH and is called directly.')
    parser.add_argument('-z', '--seed', type=int, default=None, help='Randomization seed used in fastme.')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    in_path, out_path, prefix, save_dm, rec_method, use_breakpoint, dist_type, use_total, use_NNI, fastme_path, seed = args.input, args.output, args.prefix, args.save_dm, args.rec_method, args.breakpoint, args.dist_type, args.total_cn, args.use_NNI, args.fastme_path, args.seed

    metrics = {
        'manhattan': distance_manhattan,
        'euclidean': distance_euclidean,
        'root': distance_root,
        'log': distance_log
    }

    start_time = time.time()
    part1_time = part2_time = part3_time = None

    if not out_path:
        out_path = os.getcwd()
    if not os.path.isdir(out_path):
        os.makedirs(out_path)
    if dist_type not in list(metrics.keys()):
        print('Invalid distance type.')
        return
    
    if not prefix:
        if use_breakpoint:
            prefix = f'breakpoint_{dist_type}'
        else:
            prefix = f'standard_{dist_type}'
        if rec_method:
            prefix += f'_{rec_method}'
            if rec_method == 'balME' or rec_method == 'olsME':
                if use_NNI:
                    prefix += '_nni'

    print('Reading data...', end="", flush=True)
    if in_path[-4:] == '.tsv':
        profiles_df = process_data(in_path, use_total)
        print('Done')
        part1_time = time.time()
        print('Computing distance matrix...', end="", flush=True)
        dist_matrix, cell_names = compute_distmatrix(profiles_df, use_breakpoint, metrics[dist_type], use_total)
        dm_prefix = os.path.join(out_path, prefix + '_dm.phy')
        save_distmatrix(dist_matrix, cell_names, dm_prefix)
        print('Done')
        part2_time = time.time()
    elif in_path[-4:] == '.phy':
        dm_prefix = in_path
        save_dm = True
        print('Done')
        part1_time = start_time
    if rec_method:
        print('Constructing tree...', end="", flush=True)
        if rec_method in ['NJ', 'uNJ', 'balME', 'olsME']:
            tree_prefix = os.path.join(out_path, prefix + '_tree.nwk')
            construct_tree_fast(dm_prefix, save_dm, tree_prefix, rec_method, use_NNI, fastme_path, seed)
            part3_time = time.time()
        else:
            print('Invalid reconstruction method.')
            return
    print('Done')

    info_prefix = os.path.join(out_path, prefix + '_info.txt')
        
    with open(info_prefix, 'w+') as f:
        f.write('Distance: ' + '\t' + dist_type + '\n')
        if rec_method:
            f.write('Method: ' + '\t' + rec_method + '\n')
        else:
            f.write('Method: ' + '\t' + 'None' + '\n')
        f.write('Breakpoint: ' + '\t' + str(use_breakpoint) + '\n')
        f.write('Total copy numbers: ' + '\t' + str(use_total) + '\n')
        if in_path[-4:] != '.phy':
            f.write('Reading time elapsed: ' + '\t' + str(part1_time - start_time) + '\n')
        if part2_time:
            f.write('Dm computation time elapsed: ' + '\t' + str(part2_time - part1_time) + '\n')
            if part3_time:
                f.write('Tree construction time elapsed: ' + '\t' + str(part3_time - part2_time) + '\n')
        else:
            f.write('Tree construction time elapsed: ' + '\t' + str(part3_time - part1_time) + '\n')
        f.write('Total time elapsed: ' + '\t' + str(time.time() - start_time) + '\n')
        
    if rec_method != None:
        f1, f2 = open(info_prefix, 'a+'), open(dm_prefix + '_fastme_stat.txt', 'r')
        f1.write('\n\nfastme stats\n')
        f1.write(f2.read())
        f1.seek(0)
        f1.close()
        f2.close()
        os.remove(dm_prefix + '_fastme_stat.txt')

if __name__ == '__main__':
    main()