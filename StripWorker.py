import numpy as np
import pandas as pd

VALID_NUCLEOTIDES = np.array(["A", "T", "G", "C"], dtype="S1")

class StripWorker:
    adapters: []
    adapters_len: 0 # max length of an adapter

    def __init__(self, strip_strings):
        self.adapters = [] # convert adapters to byte arrays
        for i, s in enumerate(strip_strings):
            h = np.fromstring(s.strip(), dtype="S1")
            if not np.isin(h, VALID_NUCLEOTIDES).all():
                raise Exception("Adapter " + str(i) + " contains invalid nucleotides (not A,T,G,C)")
            self.adapters.insert(i, h)

        self.adapters_len = max([len(seq) for seq in self.adapters])

    def check_adapters(self, seq, dis_min):
        # prepare slices of sequences from start and end as ndarray to speed up search
        seq_slice_start = np.array([q[:self.adapters_len] for q in seq])
        seq_slice_end = np.array([q[-self.adapters_len:] for q in seq])
        if min([len(q) for q in seq]) <= self.adapters_len:
            raise Exception("Adapter sequences are longer than the shortest fasta chain - this is not supported")

        slice_pos = []
        for i, adapter in enumerate(self.adapters):
            # matrix - if there is a match between chain and an adapter offset by column id
            eq_matrix_pos = np.zeros((len(seq), self.adapters_len + 1), dtype=bool)
            eq_matrix_neg = np.zeros((len(seq), self.adapters_len + 1), dtype=bool)
            for offset in range(1, len(adapter) + 1):
                x_start = seq_slice_start[:, :offset] # which of first/last OFFSET elements of chain
                x_end = seq_slice_end[:, :offset]
                start_equal = np.equal(x_start, adapter[-offset:]) # are equal to first/last OFFSET elements of adapter
                end_equal = np.equal(x_end, adapter[:offset])
                eq_matrix_pos[:, offset] = np.all(start_equal, axis=1) # check if they are all equal
                eq_matrix_neg[:, offset] = np.all(end_equal, axis=1)
            pos_dp = pd.DataFrame(np.where(eq_matrix_pos), index=["Sequence", "Offset"]).T # get positions where they are all equal
            pos_neg = pd.DataFrame(np.where(eq_matrix_neg), index=["Sequence", "Offset"]).T
            pos_dp = pos_dp.groupby(pos_dp["Sequence"], as_index = False).max() # select the longest match
            pos_neg = pos_neg.groupby(pos_neg["Sequence"], as_index = False).max()
            slice_pos.insert(i, (pos_dp, pos_neg))
        return slice_pos

    # filter search results by minimal inclusion, copies the results
    def filter_results(self, slice_pos, dis_min):
        n_pos = []
        for i in range(len(slice_pos)):
            e_start = slice_pos[i][0]
            n_start = e_start.drop(e_start[e_start["Offset"] < dis_min].index)
            e_end = slice_pos[i][1]
            n_end = e_end.drop(e_end[e_end["Offset"] < dis_min].index)
            n_pos.insert(i, (n_start, n_end))
        return n_pos

    def slice_seq(self, data, slice_pos, return_finals=False):
        # join slice points of all adapters and select the largest
        j_table_pos = (pd.concat([p[0] for p in slice_pos], ignore_index=True)
                       .groupby("Sequence", as_index = False).max())
        j_table_neg = (pd.concat([p[1] for p in slice_pos], ignore_index=True)
                       .groupby("Sequence", as_index = False).max())

        res_data = []
        for i, seqs in enumerate(data):
            for index, values in j_table_pos.iterrows():
                t = int(values["Sequence"])
                seqs[t] = seqs[t][int(values["Offset"]):]  # slice from the start
            for name, values in j_table_neg.iterrows():
                t = int(values["Sequence"])
                seqs[t] = seqs[t][:-int(values["Offset"])]  # slice from the end
            res_data.insert(i, seqs)

        if not return_finals:
            return (res_data,)
        else:
            return res_data, j_table_pos, j_table_neg
