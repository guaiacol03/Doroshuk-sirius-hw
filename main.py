import argparse
import os

import numpy as np
import pandas as pd

from ExportHelper import ExportHelper
from StatWorker import StatWorker
from StripWorker import StripWorker

VALID_NUCLEOTIDES = np.array(["A", "T", "G", "C"], dtype="S1")


parser = argparse.ArgumentParser()
parser.add_argument('input', type=str,
                    help='Input fastq file')
parser.add_argument('--noask', action='store_true', help='Don\'t ask for confirmation')
parser.add_argument('--export_csv', action='store_true', help='Export csv file with results')
parser.add_argument('--export_merge', action='store_true', help='Export graphs as a single file')
parser.add_argument('--export', type=str,
                    help='Folder to export results to (all export args are ignored otherwise)', required=False)
parser.add_argument('--export_dpi', type=int, default=100,
                    help='Resolution to plot graphs at (default 100)')
parser.add_argument('--adapters', type=str,
                    help='File with adapters to remove (if none skip)', required=False)
parser.add_argument('--export_fastq', action='store_true', help='Export trimmed fastq file')
parser.add_argument('--adapters_len', type=int,
                    help='Minimal length of adapters (>=3)', default=3, required=False)
parser.add_argument('--mode', choices=['base', 'strip', 'plot'], nargs='+',
                    help='Modes to run the program (space-separated)', required=True)
args = parser.parse_args()

def run():
    global parsed
    # open and read file
    if not os.path.isfile(args.input):
        raise Exception('Input file does not exist')
    with open(args.input) as file:
        input = np.array(file.readlines())
    print("Reading fasta file at " + args.input)

    try:
        input = input.reshape((-1, 4)) # convert to records of 4 strings
    except ValueError as e:
        raise Exception("Fasta file has uneven number of strings, likely corrupted or incomplete")
    checkRow = input[:, 2]
    if not np.all(checkRow == '+\n'): # check if 3rd string is "+"
        raise Exception("Fasta file has records without \"+\", likely corrupted or incomplete")
    input = np.delete(input, 2, 1)
    input = np.strings.strip(input)
    print("Format correct, started parsing")

    # transform to dataframe
    parsed = pd.DataFrame(input, index=None, columns=["Data", "Sequence", "Quality"])
    # process sequences
    parsed.Sequence = parsed.Sequence.apply(lambda seq: np.fromiter(seq, dtype="S1")) # nucleotides as byte array
    checkNucs = np.append(VALID_NUCLEOTIDES, b'N')
    nucCheckRow = parsed.Sequence.apply(lambda seq: np.isin(seq, checkNucs).all()) # check if all nucleotides are valid
    failedRows = nucCheckRow[nucCheckRow == False].index.values
    if len(failedRows) > 0:
        raise Exception("Rows " + str(failedRows) + " contain invalid nucleotides (not A,T,G,C or N)")
    # process quality
    parsed.Quality = parsed.Quality.apply(lambda seq: np.fromiter((ord(sym) for sym in seq), dtype=int))
    # def conv_meta_to_dict(meta):
    #     spl = meta.split(" ")
    #     if spl[0][0] != "@":
    #         return None
    #
    #     exp = dict()
    #     exp["Index"] = spl[0]
    #     for e in spl[1:]:
    #         espl = e.split("=")
    #         if len(espl) != 2:
    #             return None
    #         exp[espl[0]] = espl[1]
    #     return exp
    #
    # parsed.Data = parsed.Data.apply(conv_meta_to_dict)
    # failedRows = parsed[parsed["Data"].isna()].index.values
    # if len(failedRows) > 0:
    #     raise Exception("Rows " + str(failedRows) + " contain invalid metadata")
    print("Fasta file parsed without errors\n")

arg_exp_frame = False
export_path = None
if args.export is not None:
    export_path = args.export
    if os.path.exists(export_path):
        if not os.path.isdir(export_path):
            raise Exception("Export path is not a directory")
    else:
        os.mkdir(export_path)

    arg_exp_frame = args.export_csv
    if arg_exp_frame: exp = pd.DataFrame(columns=["Property", "Value"])

run()
if 'base' in args.mode:
    v = StatWorker(parsed)
    print("-- General statistics --")
    print("Input file: " + args.input.split('/')[-1])
    x = ExportHelper.print_stat(parsed, v, export_frame = arg_exp_frame)
    if export_path is not None:
        if arg_exp_frame:
            exp = pd.concat((exp, x), ignore_index=True)
if 'strip' in args.mode:
    if args.adapters is None:
        raise Exception("Mode \"strip\" requires adapters")

    if not os.path.isfile(args.adapters): # read adapters file
        raise Exception('Adapters file does not exist')
    with open(args.adapters) as strip_file:
        strip_strings = np.array(strip_file.readlines())
    stripper = StripWorker(strip_strings)

    _adapters_len = args.adapters_len
    if _adapters_len < 3:
        print('Adapters length is too small, defauling to 3')

    s_ret = stripper.check_adapters(parsed.Sequence.values, _adapters_len)
    while True:
        ret = stripper.filter_results(s_ret, _adapters_len)
        ExportHelper.print_slicing_data(stripper, ret)
        if args.noask:
            break

        print("\nOK or change minimal length? [y/n]")
        ans = input()
        if ans in ['y', 'Y']:
            break
        elif ans in ['n', 'N']:
            print("Enter new minimal length:")
            _adapters_len = int(input())
            if _adapters_len < 3:
                print('Adapters length is too small, defauling to 3')
                _adapters_len = 3
        else:
            print('Invalid input')
    if arg_exp_frame: exp.loc[len(exp)] = {"Property": "Min. inclusion", "Value": _adapters_len}

    f = stripper.slice_seq([parsed.Sequence.values, parsed.Quality.values], ret, return_finals=arg_exp_frame)
    if arg_exp_frame:
        parsed["Sequence"] = f[0][0]
        parsed["Quality"] = f[0][1]
        exp = pd.concat((exp, ExportHelper.export_slicing_res(stripper, ret, (f[1], f[2]))), ignore_index=True)
    else:
        parsed["Sequence"] = f[0][0]
        parsed["Quality"] = f[0][1]
    print("Slicing finished")

    if args.export_fastq and (export_path is not None):
        print("Composing fastq file...")
        ExportHelper.compose_fastq(parsed, export_path + "/fastq_sliced.fastq")
        print("Saved successfully")
if ('plot' in args.mode) and (export_path is not None):
    ex_dpi = args.export_dpi
    if not args.export_merge:
        ExportHelper.plot_save(v.plot_len_distribution, export_path + "/pl_len_distribution.png", dpi=ex_dpi)
        ExportHelper.plot_save(v.plot_gc_distribution, export_path + "/pl_gc_distribution.png", dpi=ex_dpi)
        ExportHelper.plot_save(v.plot_nucl_positions, export_path + "/pl_nucl_positions.png", dpi=ex_dpi)
        ExportHelper.plot_save(v.plot_qual_distribution, export_path + "/pl_qual_distribution.png", dpi=ex_dpi)
        ExportHelper.plot_save(v.plot_qual_positions, export_path + "/pl_qual_positions.png", dpi=ex_dpi)
    else:
        ExportHelper.plot_pack_save(
            [
                v.plot_len_distribution,
                v.plot_gc_distribution,
                v.plot_nucl_positions,
                v.plot_qual_distribution,
                v.plot_qual_positions
            ], export_path + "/pl_merge.png", dpi=ex_dpi)

if arg_exp_frame:
    exp.to_csv(export_path + "/stats.csv")
