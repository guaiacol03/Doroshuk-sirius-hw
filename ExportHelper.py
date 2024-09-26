import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

VALID_NUCLEOTIDES = np.array(["A", "T", "G", "C"], dtype="S1")

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
        print("Records: " + str(data.shape[0]))
        if export_frame: exp.loc[len(exp)] = {"Property": "Total sequences", "Value": data.shape[0]}

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

    @staticmethod
    def compose_fastq(tab, save_loc):
        exp_strings = []
        for i, rec in tab.iterrows():
            exp_strings.append(rec["Data"])
            exp_strings.append("".join([d.decode('utf-8') for d in rec["Sequence"]]))
            exp_strings.append("+")
            exp_strings.append("".join([chr(v) for v in rec["Quality"]]))
        with open(save_loc, "w") as out:
            out.write("\n".join(exp_strings))
