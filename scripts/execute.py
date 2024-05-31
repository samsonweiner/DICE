import argparse
from utilities import *


method_paths = {
    'medicc2': 'medicc2',
    'medalt': '../../tools/MEDALT/',
    'cnp2cnp': '../../tools/cnp2cnp/cnp2cnp.py',
    'local': '../methods/ICE.py',
    'sitka': '../../tools/sitkatree/sitka/build/install/nowellpack/bin/',
    'lazac': '../../tools/lazac-copy-number/'
}

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, required=True, help='Path to dataset.')
parser.add_argument('-p', '--prepare_inputs', action='store_true', help='Prepare method inputs.')
parser.add_argument('-n', '--num-processor', type=int, default=1)
parser.add_argument('-t', '--totalCN', action='store_true')
parser.add_argument('-d', '--DICE', action='store_true')
parser.add_argument('-0', '--medicc2', action='store_true')
parser.add_argument('-1', '--medalt', action='store_true')
parser.add_argument('-2', '--cnp2cnp', action='store_true')
parser.add_argument('-3', '--sitka', action='store_true')
parser.add_argument('-4', '--lazac', action='store_true')
args = parser.parse_args()

def prepare_data(args):
    if args.medicc2:
        prepare_medicc2(args.input, totalCN = args.totalCN)
    if args.medalt:
        prepare_MEDALT(args.input, totalCN = args.totalCN)
    if args.cnp2cnp:
        prepare_cnp2cnp(args.input, totalCN = args.totalCN)
    if args.sitka:
        prepare_sitka(args.input, totalCN = args.totalCN)
    if args.lazac:
        prepare_lazac(args.input, totalCN = args.totalCN)
    
def run_methods(args):
    if args.dice:
        run_variants(method_paths['local'], args.input, totalCN = args.totalCN)
    if args.medicc2:
        run_medicc2(method_paths['medicc2'], args.input, args.num_processor, totalCN = args.totalCN)
    if args.medalt:
        run_MEDALT(method_paths['medalt'], args.input)
    if args.cnp2cnp:
        run_cnp2cnp(method_paths['cnp2cnp'], args.input)
    if args.sitka:
        run_sitka(method_paths['sitka'], args.input)
    if args.lazac:
        run_lazac(method_paths['lazac'], args.input)

if args.input[-1] != '/':
    args.input += '/'

if args.prepare_inputs:
    prepare_data(args)

run_methods(args)
