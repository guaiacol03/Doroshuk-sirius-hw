import argparse
import os
from os.path import exists

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

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
parser.add_argument('--adapters_len', type=int,
                    help='Minimal length of adapters (>=3)', default=3, required=False)
parser.add_argument('--mode', choices=['base', 'strip'], nargs='+',
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
    nucCheckRow = parsed.Sequence.apply(lambda seq: np.isin(seq, VALID_NUCLEOTIDES).all()) # check if all nucleotides are valid
    failedRows = nucCheckRow[nucCheckRow == False].index.values
    if len(failedRows) > 0:
        raise Exception("Rows " + str(failedRows) + " contain invalid nucleotides (not A,T,G,C)")
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

class ExportHelper:
    @staticmethod
    def plot_save(delegate, save_loc, **kwargs):
        plt.figure()
        fig, ax = plt.subplots()
        delegate(ax)
        plt.savefig(save_loc, **kwargs)
        plt.close()

    @staticmethod
    def plot_pack_save(delegates, save_loc, **kwargs):
        fig, axs = plt.subplots(2, 3, figsize=(16, 8))
        del_position = 0
        for i in range(axs.shape[0]):
            for j in range(axs.shape[1]):
                if del_position < len(delegates):
                    delegates[del_position](axs[i, j])
                    del_position += 1
                else:
                    fig.delaxes(axs[i, j])
        plt.tight_layout()
        plt.savefig(save_loc, **kwargs)

    @staticmethod
    def print_stat(data, worker, export_frame=False):
        if export_frame : exp = pd.DataFrame(columns=["Property", "Value"])
        print("-- General statistics --")
        print("Input file: " + args.input.split('/')[-1])
        print("Records: " + str(data.shape[0]))

        lbins = np.bincount(worker.SeqLengths)
        print("Most common length: " + str(np.argmax(lbins)))
        if export_frame : exp.loc[len(exp)] = {"Property": "Most common length", "Value": np.argmax(lbins)}

        total_nuc = np.sum(worker.NuclRowFrequency, axis=0)
        total_qua = np.sum(total_nuc)
        print("GC quantity: " + str(np.round(100 * (float(total_nuc[2] + total_nuc[3]) / total_qua), 2)) + "%")
        if export_frame : exp.loc[len(exp)] = {"Property": "GC quantity",
                                               "Value": float(total_nuc[2] + total_nuc[3]) / total_qua}

        print("Nucleotide quantities: ")
        for i, value in enumerate(VALID_NUCLEOTIDES):
            print("\t" + value.decode('utf-8') + ": " + str(np.round(float(total_nuc[i]) / total_qua, 4)))
            if export_frame : exp.loc[len(exp)] = {"Property": value.decode('utf-8') + " quantity",
                                                   "Value": float(total_nuc[i]) / total_qua}

        if export_frame: return exp

    @staticmethod
    def export_slicing_res(worker, slice_tab, final_tab):
        exp = final_tab[0].merge(final_tab[1], how='outer', on="Sequence", suffixes=("_pos", "_neg"))
        return exp

    @staticmethod
    def print_slicing_data(worker, slice_tab):
        for i, adapter in enumerate(worker.adapters):
            print("For adapter " + str(i) + " (" + "".join([nucl.decode('utf-8') for nucl in adapter]) + "): ")
            print("* In " + str(slice_tab[i][0].shape[0]) + " sequences from the start")
            tP = np.bincount(slice_tab[i][0]["Offset"].values)
            for t, v in enumerate(tP):
                if v > 0:
                    print("\t" + str(t) + " symbols in " + str(v) + " sequences")
            print("* In " + str(slice_tab[i][1].shape[0]) + " sequences from the end")
            tN = np.bincount(slice_tab[i][1]["Offset"].values)
            for t, v in enumerate(tN):
                if v > 0:
                    print("\t" + str(t) + " symbols in " + str(v) + " sequences")

arg_exp_frame = False
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
    x = ExportHelper.print_stat(parsed, v, export_frame = arg_exp_frame)
    if args.export is not None:
        if arg_exp_frame:
            exp = pd.concat((exp, x), ignore_index=True)

        ex_dpi = args.export_dpi
        if not args.export_merge:
            ExportHelper.plot_save(v.plot_len_distribution, export_path + "/pl_len_distribution.png", dpi = ex_dpi)
            ExportHelper.plot_save(v.plot_gc_distribution, export_path + "/pl_gc_distribution.png", dpi = ex_dpi)
            ExportHelper.plot_save(v.plot_nucl_positions, export_path + "/pl_nucl_positions.png", dpi = ex_dpi)
            ExportHelper.plot_save(v.plot_qual_distribution, export_path + "/pl_qual_distribution.png", dpi = ex_dpi)
            ExportHelper.plot_save(v.plot_qual_positions, export_path + "/pl_qual_positions.png", dpi = ex_dpi)
        else:
            ExportHelper.plot_pack_save(
                [
                    v.plot_len_distribution,
                    v.plot_gc_distribution,
                    v.plot_nucl_positions,
                    v.plot_qual_distribution,
                    v.plot_qual_positions
                ], export_path + "/pl_merge.png", dpi = ex_dpi)
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

    ret = stripper.check_adapters(parsed.Sequence.values, _adapters_len)
    while True:
        ret = stripper.filter_results(ret, _adapters_len)
        ExportHelper.print_slicing_data(stripper, ret)
        if args.noask:
            break

        print("\nOK or change minimal length? [y/n]")
        ans = input()
        if ans in ['y', 'Y']:
            break
        elif ans in ['n', 'N']:
            _adapters_len = int(input())
            if _adapters_len < 3:
                print('Adapters length is too small, defauling to 3')
                _adapters_len = 3
        else:
            print('Invalid input')
    if arg_exp_frame: exp.loc[len(exp)] = {"Property": "Min. inclusion", "Value": _adapters_len}

    f = stripper.slice_seq(parsed.Sequence.values, ret, return_finals=arg_exp_frame)
    if arg_exp_frame:
        parsed["Sequence"] = f[0]
        exp = pd.concat((exp, ExportHelper.export_slicing_res(stripper, ret, (f[1], f[2]))), ignore_index=True)
        print(exp)
    else:
        parsed["Sequence"] = f
    print("Slicing finished")

if arg_exp_frame:
    exp.to_csv(export_path + "/stats.csv")
