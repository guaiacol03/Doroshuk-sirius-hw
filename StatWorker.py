import numpy as np
from matplotlib import pyplot as plt

VALID_NUCLEOTIDES = np.array(["A", "T", "G", "C"], dtype="S1")


class StatWorker:
    def __init__(self, dat):
        # sequences length distribution
        self.SeqLengths = dat.Sequence.apply(lambda n: np.shape(n)[0]).values

        # nucleotide distribution for each position
        max_len = np.max(self.SeqLengths)
        elong = dat.Sequence.apply( # elongate sequences to uniform length, so they can be made into ndarray
            lambda seq: np.concatenate((seq, np.full(max_len - len(seq), b'0')))).values
        elong_arr = np.array(list(elong))
        count_ax = np.zeros((max_len, 4), dtype=int)
        for i in range(max_len):
            x = np.unique(elong_arr[:, i], return_counts=True)
            for n_nuc, nuc in enumerate(VALID_NUCLEOTIDES): # fillers are dropped
                ct_pos = np.where(x[0] == nuc)[0]
                if len(ct_pos) > 0:
                    count_ax[i, n_nuc] = x[1][ct_pos[0]]
        self.NuclPosFrequency = count_ax

        # nucleotide distribution for each sequence
        count_row = np.zeros((max_len, 4), dtype=int)
        row_dat = dat.Sequence.apply(StatWorker._count_nucl_rows)
        row_arr = np.array(list(row_dat.values))
        self.NuclRowFrequency = row_arr

        # mean quality for sequences
        self.SeqMeanQuality = dat.Quality.apply(np.mean).values

        # quality count in each position
        qlong = dat.Quality.apply(
            lambda seq: np.concatenate((seq, np.full(max_len - len(seq), np.nan)))).values
        qlong_arr = np.array(list(qlong))
        self.PosMeanQuality = np.nanmean(qlong_arr, axis=0)
        print(1)

    @staticmethod
    def _count_nucl_rows(row):
        x = np.unique(row, return_counts=True)
        ret = np.ndarray(4, dtype=int)
        for n_nuc, nuc in enumerate(VALID_NUCLEOTIDES):
            ct_pos = np.where(x[0] == nuc)[0]
            if len(ct_pos) > 0:
                ret[n_nuc] = x[1][ct_pos[0]]
        return ret

    def plot_len_distribution(self, canvas):
        canvas.hist(self.SeqLengths, bins=20)
        canvas.set_title("Distr. of sequence lengths")
        canvas.set_xlabel("Length (bases)")
        canvas.set_ylabel("No. of sequences")

    def plot_gc_distribution(self, canvas):
        mod = np.apply_along_axis(lambda x: (float(np.sum(x[2:])) / np.sum(x)) * 100, axis=1, arr=self.NuclRowFrequency)
        ht = np.histogram(mod, bins=30)
        canvas.bar(ht[1][:-1], ht[0], align="edge", width=0.75)
        canvas.set_title("Distr. of GC content")
        canvas.set_xlabel("GC content (%)")
        canvas.set_ylabel("No. of sequences")
        print(1)

    def plot_nucl_positions(self, canvas):
        arr_p = np.apply_along_axis(lambda row: [float(k) / np.sum(row) for k in row], axis=1, arr=self.NuclPosFrequency)
        arr_limit = np.min(np.where(arr_p > 0.75)[0]) # stop plotting when distribution grows abnormal (too few reads long enough)
        for i in range(arr_p.shape[1]): # for each nucleotide
            pl_data = np.convolve(arr_p[:arr_limit, i], np.ones(10) / 10, mode='same') # smooth out with moving average
            # pl_data[:10] = arr_p[:10, i] # remove segment where moving average grows
            canvas.plot(range(0, arr_limit), pl_data * 100, label=VALID_NUCLEOTIDES[i].decode("utf-8") + "%")
        canvas.set_ylim(0, 75)
        canvas.set_title("Distr. of nucleotides")
        canvas.set_xlabel("Sequence pos.")
        canvas.set_ylabel("Base %")
        canvas.legend()

    def plot_qual_distribution(self, canvas):
        ht = np.histogram(self.SeqMeanQuality, bins=30)
        canvas.bar(ht[1][:-1], ht[0], align="edge", width=0.25)
        canvas.set_title("Distr. of sequence mean quality")
        canvas.set_xlabel("Avg. quality (phred)")
        canvas.set_ylabel("No. of sequences")
        print(1)

    def plot_qual_positions(self, canvas):
        pval_orig = self.PosMeanQuality
        pval = np.convolve(pval_orig, np.ones(10) / 10, mode='same') # smoothen graph
        pval[:10] = pval_orig[:10] # remove segment where moving average grows
        canvas.set_title("Distr. of quality")
        canvas.set_xlabel("Sequence pos.")
        canvas.set_ylabel("Avg. quality (phred)")
        canvas.plot(range(pval.shape[0]), pval)
        print(1)

