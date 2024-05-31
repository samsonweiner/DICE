import os
import subprocess
import time

import shutil
from ete3 import Tree
#import pandas as pd
from skbio import DistanceMatrix
from skbio.tree import nj

#from tree import Node, Tree
#from willowtree import *
from collections import deque

def prepare_medicc2(prefix, totalCN=False):
    data_path = prefix + 'profiles.tsv'
    out_dir = prefix + 'medicc2/'
    if not os.path.isdir(out_dir):
        call = subprocess.call(['mkdir', out_dir])
    out_path = prefix + 'medicc2/medicc_input.tsv'

    data = {}
    if totalCN:
        headers = ['sample_id', 'chrom', 'start', 'end', 'cn_a']
    else:
        headers = ['sample_id', 'chrom', 'start', 'end', 'cn_a', 'cn_b']
    f = open(data_path)
    lines = f.readlines()
    f.close()
    for line in lines[1:]:
        cell, chrom, start, end, CN = line[:-1].split('\t')
        chrom = chrom[chrom.index('r')+1:]
        if totalCN:
            row = [cell, 'chrom' + chrom, start, end, CN]
        else:
            cn_a, cn_b = CN[:CN.index(',')], CN[CN.index(',')+1:]
            row = [cell, 'chrom' + chrom, start, end, cn_a, cn_b]
        if cell not in data:
            data[cell] = [row]
        else:
            data[cell].append(row)

    f = open(out_path, 'w+')
    f.write('\t'.join(headers) + '\n')
    for cell, lines  in data.items():
        for line in lines:
            f.write('\t'.join(line) + '\n')
    f.close()

def prepare_cnp2cnp(prefix, totalCN = False, asIs = True):
    data_path = prefix + 'profiles.tsv'
    out_dir = prefix + 'cnp2cnp/'
    if not os.path.isdir(out_dir):
        call = subprocess.call(['mkdir', out_dir])
    out_path = prefix + 'cnp2cnp/cnp2cnp_input.fa'

    f = open(data_path)
    lines = f.readlines()
    f.close()

    data = {}
    loci = []
    for line in lines:
        curline = line[:-1].split('\t')
        if curline[0] != 'CELL':
            cell = curline[0]
            chrom = curline[1]
            if totalCN or asIs:
                start = int(curline[2])
                end = int(curline[3])
            else:
                start = int(int(curline[2]) / 10)
                end = int(int(curline[3]) / 10)
            CNs = curline[-1]
            if totalCN:
                cn_a = int(CNs)
                cn_b = 0
            else:
                cn_a = int(CNs[:CNs.index(',')])
                cn_b = int(CNs[CNs.index(',')+1:])
            if cell not in data:
                data[cell] = {}
            locus = chrom[3:] + '_' + str(start+1) + '_' + str(end)
            if locus not in loci:
                loci.append(locus)
            data[cell][locus] = str(cn_a + cn_b)

    f = open(out_path, 'w+')
    for cell in data.keys():
        profile = [data[cell][locus] for locus in loci]
        data[cell] = profile
        f.write('>' + cell + '\n')
        f.write(','.join(profile) + '\n')
    f.close()

def prepare_MEDALT(prefix, totalCN = False):
    data_path = prefix + 'profiles.tsv'
    out_dir = prefix + 'MEDALT/'
    if not os.path.isdir(out_dir):
        call = subprocess.call(['mkdir', out_dir])
    out_path = prefix + 'MEDALT/MEDALT_input.tsv'

    f = open(data_path)
    lines = f.readlines()
    f.close()

    data, cellnames = {}, set()
    for line in lines:
        curline = line[:-1].split('\t')
        if curline[0] != 'CELL':
            cell = curline[0]
            chrom = int(curline[1][3:])
            start = int(curline[2])
            end = int(curline[3])
            CNs = curline[-1]
            if totalCN:
                cn_a = int(CNs)
                cn_b = 0
            else:
                cn_a = int(CNs[:CNs.index(',')])
                cn_b = int(CNs[CNs.index(',')+1:])
            cellnames.add(cell)
            if chrom not in data:
                data[chrom] = {}
            if start not in data[chrom]:
                data[chrom][start] = {}
            data[chrom][start][cell] = str(cn_a + cn_b)

    f = open(out_path, 'w+')
    headers = ['chrom', 'chrompos']
    for cell in cellnames:
        headers.append(cell)
    f.write('\t'.join(headers) + '\n')
    for chrom in data:
        for pos in data[chrom]:
            line = [str(chrom), str(pos)]
            for cell in cellnames:
                line.append(data[chrom][pos][cell])
            f.write('\t'.join(line) + '\n')
    f.close()

def prepare_sitka(prefix, totalCN=False):
    data_path = prefix + 'profiles.tsv'
    out_dir = prefix + 'sitka/'
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    out_path = prefix + 'sitka/sitka_input.csv'

    f = open(data_path)
    lines = f.readlines()
    f.close()

    data, cellnames = {}, set()
    size = -1
    for line in lines:
        curline = line[:-1].split('\t')
        if curline[0] != 'CELL':
            cell = curline[0]
            chrom = curline[1][3:]
            if chrom[0] == '_':
                chrom = int(chrom[1:])
            else:
                chrom = int(chrom)
            start = int(curline[2]) + 1
            end = int(curline[3])
            if size == -1:
                size = end - start
            CNs = curline[-1]
            if totalCN:
                cn_a = int(CNs)
                cn_b = 0
            else:
                cn_a = int(CNs[:CNs.index(',')])
                cn_b = int(CNs[CNs.index(',')+1:])
            cellnames.add(cell)
            if chrom not in data:
                data[chrom] = {}
            if start not in data[chrom]:
                data[chrom][start] = {}
            #data[chrom][start][0][cell] = str(cn_a)
            #data[chrom][start][1][cell] = str(cn_b)
            data[chrom][start][cell] = str(cn_a + cn_b)

    cells = list(cellnames)
    with open(out_path, 'w+') as f:
        headers = ','.join([c[4:] for c in cells])
        f.write(headers + '\n')
        #for allele in [0,1]:
        for chrom in data:
            keys = list(data[chrom])
            keys.sort()
            for pos in keys:
                line = [f'{chrom}_{pos}_{pos+size}'] + [data[chrom][pos][cell] for cell in cells]
                f.write(','.join(line) + '\n')
    
    out_path2 = prefix + 'sitka/binary_input.csv'
    binsize = str(size+1)
    params = f'Rscript ../../tools/sitkatree/example/cntob.R -i {out_path} -o {out_path2} -b {binsize}'
    subprocess.run(params, shell=True)



def prepare_lazac(prefix, totalCN = False):
    data_path = prefix + 'profiles.tsv'
    out_dir = prefix + 'lazac/'
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    out_path = prefix + 'lazac/lazac_input.csv'

    data = {}
    if totalCN:
        headers = ['node', 'chrom', 'start', 'end', 'cn_a']
    else:
        headers = ['node', 'chrom', 'start', 'end', 'cn_a', 'cn_b']
    f = open(data_path)
    lines = f.readlines()
    f.close()
    for line in lines[1:]:
        cell, chrom, start, end, CN = line[:-1].split('\t')
        chrom = chrom[chrom.index('r')+1:]
        if totalCN:
            row = [cell, chrom, start, end, CN]
        else:
            cn_a, cn_b = CN[:CN.index(',')], CN[CN.index(',')+1:]
            row = [cell, chrom, start, end, cn_a, cn_b]
        if cell not in data:
            data[cell] = [row]
        else:
            data[cell].append(row)

    f = open(out_path, 'w+')
    f.write(','.join(headers) + '\n')
    for cell, lines  in data.items():
        for line in lines:
            f.write(','.join(line) + '\n')
    f.close()

def format_MEDALT_tree(tree_path):
    f = open(tree_path)
    lines = f.readlines()
    f.close()

    nodes, cell_names = {}, set()
    leftnodes, rightnodes = set(), set()
    for line in lines[1:]:
        node1, node2, dist = line.split('\t')
        if node1 != 'root':
            cell_names.add(node1)
            leftnodes.add(node1)
        if node2 != 'root':
            cell_names.add(node2)
            rightnodes.add(node2)
        if node1 not in nodes:
            nodes[node1] = Node(name=node1)
        if node2 not in nodes:
            nodes[node2] = Node(name=node2)
        nodes[node1].set_child(nodes[node2])
        nodes[node2].set_len(dist)

    if 'root' in nodes:
        t = Tree(root=nodes['root'])
    else:
        temp_root = Node(name='root')
        start_cell = nodes[list(leftnodes - rightnodes)[0]]
        temp_root.set_child(start_cell)
        start_cell.set_len(0)
        t = Tree(root=temp_root)

    #t.mut_to_binary()
    leaf_names = [node.name for node in t.iter_descendants() if not node.is_root()]
    dq = deque()
    dq.append(t.root)

    while dq:
        node = dq.popleft()
        if not node.is_leaf():
            if node.name in leaf_names:
                parent = node.parent
                children = [child for child in node.children]
                node.detach()
                new_node1 = parent.add_child()
                new_node1.set_child(node)
                for child in children:
                    child.detach()
                    new_node1.set_child(child)
                dq.append(new_node1)
            else:
                for child in node.children:
                    dq.append(child)
    return t

def run_medicc2(method_path, prefix, num_proc, totalCN = False):
    in_path = prefix + 'medicc2/medicc_input.tsv'
    out_path = prefix + 'medicc2/'

    start = time.time()
    if totalCN:
        subprocess.call([method_path, in_path, out_path, '--topology-only', '--no-plot', '-j', str(num_proc), '--total-copy-numbers', '--input-allele-columns', 'cn_a'])
    else:    
        subprocess.call([method_path, in_path, out_path, '--topology-only', '--no-plot', '-j', str(num_proc)])
    medicc2_time = time.time() - start
    
    f = open(out_path + 'times.txt', 'w+')
    f.write(str(medicc2_time))
    f.close()

def run_cnp2cnp(method_path, prefix):
    in_path = prefix + 'cnp2cnp/cnp2cnp_input.fa'
    out_path = prefix + 'cnp2cnp/dist.phy'

    start = time.time()
    subprocess.call(['python3', method_path, '-m', 'matrix', '-i', in_path, '-o', out_path])
    cnp2cnp_time = time.time() - start

    f = open(out_path)
    lines = f.readlines()[1:]
    f.close()

    data = []
    ids = []
    for line in lines:
        row = line.strip().split(' ')
        ids.append(row[0])
        row = row[1:]
        while not row[0].isdigit():
            row.pop(0)
        #final_index = max(index for index, item in enumerate(row) if item == '')
        #row = row[final_index+1:]
        data.append(row)
    
    dm = DistanceMatrix(data, ids)
    tree_str = nj(dm, result_constructor=str)
    f = open(prefix + 'cnp2cnp/tree.nwk', 'w+')
    f.write(tree_str)
    f.close()

    f = open(prefix + 'cnp2cnp/times.txt', 'w+')
    f.write(str(cnp2cnp_time))
    f.close()

def run_MEDALT(method_path, prefix):
    in_path = prefix + 'MEDALT/MEDALT_input.tsv'
    out_dir = prefix + 'MEDALT/'


    start = time.time()
    subprocess.call(['python2', method_path + 'scTree.py', '-P', method_path, '-I', in_path, '-D', 'D', '-G', 'hg19', '-O', out_dir])
    medalt_time = time.time() - start

    t = format_MEDALT_tree(prefix + 'MEDALT/CNV.tree.txt')
    t.save(prefix + 'MEDALT/lineage_tree.nwk')

    f = open(out_dir + 'times.txt', 'w+')
    f.write(str(medalt_time))
    f.close()

def run_sitka(method_path, prefix):
    in_path = prefix + 'sitka/binary_input.csv'
    out_dir = prefix + 'sitka/'

    params = f'{method_path}corrupt-straighten --input {in_path}'
    subprocess.run(params, shell=True)
    os.rename('results/latest/output.csv', f'{out_dir}output.csv')

    params = f'{method_path}corrupt-filter --input {out_dir}output.csv'
    subprocess.run(params, shell=True)
    os.rename('results/latest/filtered.csv', f'{out_dir}filtered.csv')

    params = [
       f'{method_path}corrupt-infer-with-noisy-params',
        '--model.binaryMatrix', f'{out_dir}filtered.csv',
        '--model.globalParameterization', 'true',
        '--model.fprBound', '0.1',
        '--model.fnrBound', '0.5',
        '--engine', 'PT',
        '--engine.initialization', 'FORWARD',
        '--engine.nScans', '1000',
        '--engine.nPassesPerScan', '1',
        '--engine.nChains', '1',
        '--engine.nThreads.fraction', '0.9'
    ]
    params = ' '.join(params)
    subprocess.run(params, shell=True)
    os.rename('results/latest/samples/phylo.csv', f'{out_dir}phylo.csv')

    params = f'{method_path}corrupt-average --csvFile {out_dir}phylo.csv --logisticTransform false'
    subprocess.run(params, shell=True)
    os.rename('results/latest/average.csv', f'{out_dir}average.csv')

    params = f'{method_path}corrupt-greedy --tipInclusionProbabilities ReadOnlyCLMatrix {out_dir}average.csv'
    subprocess.run(params, shell=True)
    t = Tree('results/latest/tree.newick', format=1)
    keep = [l.name for l in t if 'cell' in l.name]
    t.prune(keep)
    for leaf in t:
        leaf.name = leaf.name.replace('_', '')
        leaf.name = leaf.name.replace('cell', 'leaf')
    t.write(format=9, outfile=f'{out_dir}sitka_tree.nwk')

    os.remove(f'{out_dir}output.csv')
    os.remove(f'{out_dir}filtered.csv')
    os.remove(f'{out_dir}phylo.csv')
    os.remove(f'{out_dir}average.csv')
    shutil.rmtree('results')
#def compute_consensus(prefix):
#    t1 = Tree(newick=prefix + 'standard/man_tree.nwk')
#    t2 = Tree(newick=prefix + 'breakpoint/man_tree.nwk')
#    ct = strict_consensus(t1, t2, rooted=False)
#    out_str = str(ct)

#    f = open(prefix + 'consensus_tree.nwk', 'w+')
#    f.write(out_str)
#    f.close()
def run_lazac(method_path, prefix):
    in_path = prefix + 'lazac/lazac_input.csv'
    out_path = prefix + 'lazac/out'

    start = time.time()
    params = f'lazac distance {in_path} -o {out_path}'
    subprocess.run(params, shell=True)

    params = f'python {method_path}scripts/nj.py {out_path}_dist_matrix.csv --output {out_path}_nj_tree.nwk'
    subprocess.run(params, shell=True)
    params = f'python {method_path}scripts/resolve_polytomies.py {out_path}_nj_tree.nwk --output {out_path}_nj_binary_tree.nwk'
    subprocess.run(params, shell=True)
    params = f'lazac nni {in_path} {out_path}_nj_binary_tree.nwk -a 2 -o {out_path}_rcnt'
    subprocess.run(params, shell=True)
    lazac_time = time.time() - start
    
    f = open(prefix + 'lazac/times.txt', 'w+')
    f.write(str(lazac_time))
    f.close()

def run_variants(method_path, prefix, totalCN = False):
    out_dir = os.path.join(prefix, 'variants')

    #print('deleting')
    allfiles = os.listdir(out_dir)
    for filename in allfiles:
        if 'balME' in filename or 'olsME' in filename:
            os.remove(os.path.join(out_dir, filename))

    dists = ['euclidean', 'manhattan', 'root', 'log']
    methods = ['balME', 'olsME']

    #print('making dm')
    for dist in dists:
        cmd = ['python3', method_path, '-i', prefix + 'profiles.tsv', '-o', out_dir, '-d', dist]
        if totalCN:
            cmd.append('-t')
        subprocess.call(cmd)
        #if dist != 'euclidean':
        cmd.append('-b')
        subprocess.call(cmd)

    #print('making trees')
    for dist in dists:
        for method in methods:
            inpath = os.path.join(out_dir, f'standard_{dist}_dm.phy')
            cmd = ['python3', method_path, '-i', inpath, '-o', out_dir, '-d', dist, '-m', method]
            if totalCN:
                cmd.append('-t')
            subprocess.call(cmd)
            #cmd.append('-n')
            #subprocess.call(cmd)
        
    #for dist in dists[1:]:
    #    for method in methods:
            inpath = os.path.join(out_dir, f'breakpoint_{dist}_dm.phy')
            cmd = ['python3', method_path, '-i', inpath, '-o', out_dir, '-d', dist, '-m', method, '-b']
            if totalCN:
                cmd.append('-t')
            subprocess.call(cmd)
            #cmd.append('-n')
            #subprocess.call(cmd)

    #subprocess.call(['python3', method_path, '-i', prefix + 'profiles.tsv', '-o', out_dir, '-d', dist, '-m', method, '-t'])
                #subprocess.call(['python3', method_path, '-i', prefix + 'profiles.tsv', '-o', out_dir, '-d', dist, '-m', method, '-n', '-t'])
            #else:
                #subprocess.call(['python3', method_path, '-i', prefix + 'profiles.tsv', '-o', out_dir, '-d', dist, '-m', method])
                #subprocess.call(['python3', method_path, '-i', prefix + 'profiles.tsv', '-o', out_dir, '-d', dist, '-m', method, '-n'])
    methods = ['fastNJ', 'fastuNJ']
    for method in methods:
        inpath = os.path.join(out_dir, f'breakpoint_euclidean_dm.phy')
        cmd = ['python3', method_path, '-i', inpath, '-o', out_dir, '-d', 'euclidean', '-m', method, '-b']
        if totalCN:
            cmd.append('-t')
        subprocess.call(cmd)

