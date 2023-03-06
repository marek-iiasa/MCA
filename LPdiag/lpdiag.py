"""
Prototype of simple analysis of the MPS-format file
Written by Marek Makowski, ECE Program of IIASA, in March 2023
"""

# import sys		# needed for sys.exit()
import typing
import os
import numpy as np
import pandas as pd
# from datetime import datetime as dt
# from datetime import timedelta as td
# import matplotlib.pyplot as plt
# from matplotlib import cm
# from matplotlib import colors
# from matplotlib.ticker import LinearLocator    # needed for ax.set_major_locator
# import seaborn as sns
# sns.set()   # settings for seaborn plotting style


class LPdiag:
    """Process the MPS-format input file and provide its basic diagnotics.

    The diagnostics currently includes:
    - handling formal errors of the MPS file
    - basic statistics of the matrix coefficients.

    Attributes
    ----------
    rep_dir : str
        sub-directory for reports (text and plots)
    """
    def __init__(self, rep_dir):
        self.rep_dir = rep_dir    # subdirectory for reports (text, in future also plots)
        self.fname = 'undefined'    # MPS input file (to be defined, if/when rd_mps() is called)
        self.pname = 'undefined'    # problem name
        if not os.path.exists(self.rep_dir):
            os.makedirs(self.rep_dir, mode=0o755)

        # initialize the vars for processing
        self.n_lines = 0        # number of processed lines
        self.row_names = {}     # names of rows and their seq_numbers
        self.row_types = []     # types of rows
        self.gf_seq = -1        # sequence_no of the goal function (objective) row: equal = -1, if undefined
        self.col_names = {}     # names of columns and their seq_numbers
        self.mat = pd.DataFrame(columns=['row', 'col', 'val'])   # space for the LP matrix
        self.id_rows = []       # ids of rows
        self.id_cols = []       # ids of cols
        self.vals = []          # values of the matrix

    def rd_mps(self, fname):   # process the MPS file
        print(f'\nReading MPS-format file {fname}.')
        self.fname = fname
        sections = ['NAME', 'ROWS', 'COLUMNS', 'RHS', 'RANGES', 'BOUNDS', 'ENDATA']
        req_sect = [True, True, True, False, False, False, True]    # required/optional MPS sections
        row_typ = ['N', 'E', 'G', 'L']  # types of rows
        n_section = 0   # seq_no of the currently processed MPS-file section
        next_sect = 0   # seq_no of the next (to be processed) MPS-file section
        col_curr = ''   # current column (initialized to an illegal empty name)
        with open(self.fname, 'r') as reader:
            for n_line, line in enumerate(reader):
                line = line.rstrip('\n')
                # print(f'line {line}')
                if line[0] == '*' or len(line) == 0:    # skip commented and empty lines
                    continue
                words = line.split()
                n_words = len(words)
                if line[0] == ' ':  # continue the current MPS section
                    if n_section == 2:  # columns/matrix (first here because most frequently used
                        # print(f'processing line no {n_line}, n_words {n_words}: {line}')
                        assert n_words == 3 or n_words == 5, f'matrix element (line {n_line}) has {n_words} words.'
                        col_name = words[0]
                        if col_name != col_curr:    # new column
                            assert col_name not in self.col_names, f'duplicated column name: {col_name} (line {n_line})'
                            col_seq = len(self.col_names)
                            self.col_names.update({col_name: col_seq})
                            col_curr = col_name
                        row_name = words[1]
                        row_seq = self.row_names.get(row_name)
                        assert row_seq is not None, f'unknown row name {row_name} (line {n_line}).'
                        val = float(words[2])
                        assert type(val) == float, f'string  {words[2]} (line {n_line}) is not a number.'
                        # print(f' matrix element ({row_seq}, {col_seq}) = {val}')
                        # the next two lines takes far too long for large matrices; thus tmp-store in three lists
                        # df2 = pd.DataFrame({'row': row_seq, 'col': col_seq, 'val': val}, index=list(range(1)))
                        # self.mat = pd.concat([self.mat, df2], axis=0, ignore_index=True)
                        self.id_rows.append(row_seq)
                        self.id_cols.append(col_seq)
                        self.vals.append(val)
                        if n_words > 3:     # proccess second matrix element in the same MPS row
                            assert n_words == 5, f'line {n_line}) has {n_words} words, five words needed for' \
                                                 f'defining second element in the same MPS line.'
                            row_name = words[3]
                            row_seq = self.row_names.get(row_name)
                            assert row_seq is not None, f'unknown row name {row_name} (line {n_line}).'
                            val = float(words[4])
                            assert type(val) == float, f'string  {words[4]} (line {n_line}) is not a number.'
                            self.id_rows.append(row_seq)
                            self.id_cols.append(col_seq)
                            self.vals.append(val)
                    elif n_section == 1:  # rows
                        assert n_words == 2, f'row declaration (line {n_line}) has {n_words} words instead of 2.'
                        assert words[0] in row_typ, f'unknown row type {words[0]} (line {n_line}).'
                        assert words[1] not in self.row_names, f'duplicated row name: {words[1]} (line {n_line}).'
                        seq_id = len(self.row_names)
                        self.row_names.update({words[1]: seq_id})
                        self.row_types.append(words[0])
                        if words[0] == 'N' and self.gf_seq == -1:
                            self.gf_seq = seq_id
                            print(f'Row {words[1]} (seq_id {seq_id}) declared as the objective (goal function) row.')
                    elif n_section == 3:  # rhs
                        # TODO: process RHS
                        pass
                    elif n_section == 4:  # ranges
                        # TODO: process Ranges
                        pass
                    elif n_section == 5:    # bounds
                        # TODO: process Bounds
                        pass
                    elif n_section == 6:  # end data
                        raise Exception(f'Unexpected execution flow; needs to be explored.')
                    else:
                        raise Exception(f'Unknown section id {n_section}.')
                else:       # the first or a next section
                    print(f'next section found: {line} (line {n_line}).')
                    self.n_lines = n_line
                    if req_sect[next_sect]:     # required section must be defined in the sequence
                        assert words[0] == sections[next_sect], f'expect section {sections[next_sect]} found: {line}.'
                        if words[0] == sections[next_sect]:
                            n_section = next_sect
                            next_sect = n_section + 1
                            if n_section == 0:  # the only section declaration with the content
                                assert n_words > 1, f'problem name undefined: line {n_line} has {n_words} words.'
                                self.pname = words[1]  # store the problem name
                                print(f'Problem name: {self.pname}.')
                            continue
                        else:
                            raise Exception(f'Required MPS section {sections[n_section]} undefined or misplaced.')
                    else:   # optional section expected
                        if line == sections[next_sect]:     # found
                            print(f'Warning: values of section {sections[next_sect]} not stored yet.')
                            n_section = next_sect
                            next_sect = n_section + 1
                        else:       # expected section not found; process the section found
                            try:
                                n_section = sections.index(line)
                            except ValueError:
                                raise Exception(f'Unknown section id :{line} (line number = {n_line}).')
                            next_sect = n_section + 1
                            print(f'Warning: values of section {sections[n_section]} not stored yet.')
                        continue

        # check, if there was at least one N row (the first N row assumed to be the objective)
        assert self.gf_seq != -1, f'objective (goal function) row is undefined.'

        # create a df with the matrix coefficients
        self.mat = pd.DataFrame({'row': self.id_rows, 'col': self.id_cols, 'val': self.vals})
        self.mat['abs_val'] = abs(self.mat['val'])      # add column with absolute values of coeff.
        self.mat['log'] = np.log10(self.mat['abs_val']).astype(int)  # add col with int-part of the log10(coeffs)

        # Finish the MPS processing with the summary of its size
        print(f'\nFinished processing {self.n_lines} lines of the MPS file {self.fname}.')
        print(f'LP has: {len(self.row_names)} rows, {len(self.col_names)} cols, {len(self.mat)} non-zeros.')
        # TODO: add a summary of the GF row

    def stat(self, lo_tail=-7, up_tail=6):
        """Basic statistics of the matrix coefficients.

        Focus on distributions of magnitudes of non-zero coeff. represented by values of int(log10(abs(coeff))).
        Additionally, tails (low and upp) of the distributions are reported.

        Attributes
        ----------
        lo_tail: int
            Magnitude order of the low-tail (-7 denotes values < 10^(-6))
        up_tail: int
            Magnitude order of the upper-tail (6 denotes values >= 10^6)
        """

        print(f'\nDistribution of non-zeros values:\n{self.mat["abs_val"].describe()}')
        print(f'\nDistribution of log10(values):\n{self.mat["log"].describe()}')
        min_logv = self.mat["log"].min()
        max_logv = self.mat["log"].max()
        print(f'log10 values: min = {min_logv}, max = {max_logv}.')

        if lo_tail > up_tail:
            print(f'Overlapping distribution tails ({lo_tail}, {up_tail}) reset to 0.')
            lo_tail = up_tail = 0

        # low-tail of the distribution
        if lo_tail < min_logv:
            print(f'\nNo log10(values) in the requested low-tail (<= {lo_tail}) of the ditribution.')
        else:
            print(f'\nDistribution of log10(values) in the requested low-tail (<= {lo_tail}) of the ditribution.')
            print(f'{self.mat.loc[self.mat["log"] <= lo_tail].describe()}')
            for val in [*range(min_logv, lo_tail + 1)]:
                print(f'Number of log10(values) == {val}: {self.mat.loc[self.mat["log"] == val]["log"].count()}')
        # up-tail of the distribution
        if max_logv < up_tail:
            print(f'\nNo log10(values) in the requested upper-tail (>= {up_tail}) of the ditribution.')
        else:
            print(f'\nDistribution of log10(values) in the requested upp-tail (>= {up_tail}) of the ditribution.')
            print(f'{self.mat.loc[self.mat["log"] >= up_tail].describe()}')
            for val in [*range(up_tail, max_logv + 1)]:
                print(f'Number of log10(values) == {val}: {self.mat.loc[self.mat["log"] == val]["log"].count()}')

    def out_loc(self, small=True, thresh=-7, max_rec=500):
        """Locations of outlayers, i.e., elements having small/large coeff values.

        Locations of outlayers (in the term of the matrix coefficient values).
        The provided ranges of values in the corresponding row/col indicate potential of the simple scaling.

        Attributes
        ----------
        small: bool
            True/False for threshold of either small or large coefficients.
        thresh: int
            Magnitude of the threshold (in: int(log10(abs(coeff))), i.e. -7 denotes values < 10^(-6))
        max_rec: int
            Maximum number of processed coefficients
        """

        if small:
            df = self.mat.loc[self.mat['log'] <= thresh]
            print(f'\nLocations of {df["log"].count()} outlayers (coeff. with values of log10(values) <= {thresh}).')
        else:
            df = self.mat.loc[self.mat['log'] >= thresh]
            print(f'\nLocations of {df["log"].count()} outlayers (coeff. with values of log10(values) >= {thresh}).')
        df1 = df.sort_values('row')
        df1.reset_index()
        for n_rows, (indx, row) in enumerate(df1.iterrows()):
            assert n_rows < max_rec, f'To process all requested coeffs modify the safety limit assertion.'
            row_seq, row_name = self.ent_inf(row, True)     # row seq_id and name of the current coeff.
            col_seq, col_name = self.ent_inf(row, False)    # col seq_id and name of the current coeff.
            print(f'Coeff. ({row_seq}, {col_seq}): val = {row["val"]:.4e}, log(val) = {row["log"]:n}')
            df_row = df1.loc[df1['row'] == row_seq]  # df with elements in the same row
            print(f'\tRow {row_name} has {df_row["log"].count()} coeffs of magnitudes in [{df_row["log"].min()},'
                  f'{df_row["log"].max()}]')
            df_col = df1.loc[df1['col'] == col_seq]  # df with elements in the same col
            print(f'\tCol {col_name} has {df_col["log"].count()} coeffs of magnitudes in [{df_col["log"].min()},'
                  f'{df_col["log"].max()}]')
            # print(f'matrix elements in the same row:\n{df_row}')

    def ent_inf(self, row, by_row=True) -> typing.Tuple[int, str]:
        """Return info on the entity (either row or col) defining the selected matrix coefficient.

        Each row of the dataFrame contains definition (composed of the row_seq, col_seq, value, log(value))
        of one matrix coefficient.
        The function returns seq_id and name of either row or col of the currently considered coeff.

        Attributes
        ----------
        row: dataFrame row
            row/record of the df with the data of currently processed element
        by_row: bool
            True/False for returning the seq_id and name of the corresponding row/col
        """
        if by_row:
            ent_id = 'row'
            ent_seq = int(row[ent_id])
            names = [k for k, idx in self.row_names.items() if idx == ent_seq]
        else:
            ent_id = 'col'
            ent_seq = int(row[ent_id])
            names = [k for k, idx in self.col_names.items() if idx == ent_seq]
        assert len(names) == 1, f'List of {ent_id}_names for id {ent_seq} has {len(names)} elements.'
        return ent_seq, names[0]

    def plot_hist(self):
        """Plot histograms."""
