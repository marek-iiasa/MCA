"""
Prototype of simple analysis of the MPS-format file
"""
# import sys		# needed for sys.exit()
import os
import numpy as np
import pandas as pd
from datetime import datetime as dt
# from datetime import timedelta as td
# import matplotlib.pyplot as plt
# from matplotlib import cm
# from matplotlib import colors
# from matplotlib.ticker import LinearLocator    # needed for ax.set_major_locator
# import seaborn as sns
# sns.set()   # settings for seaborn plotting style


class LPdiag:
    def __init__(self, rep_dir):
        self.rep_dir = rep_dir    # subdirectory for reports
        self.fname = 'undefined'    # MPS input file
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
        self.fname = fname
        sections = ['NAME', 'ROWS', 'COLUMNS', 'RHS', 'RANGES', 'BOUNDS', 'ENDATA']
        req_sect = [True, True, True, False, False, False, True]
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
                        assert n_words == 3, f'matrix element (line {n_line}) has {n_words} words, cannot process it.'
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
                        pass
                    elif n_section == 4:  # ranges
                        pass
                    elif n_section == 5:    # bounds
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
                            print(f'Values of section {sections[next_sect]} not stored yet.')
                            n_section = next_sect
                            next_sect = n_section + 1
                        else:       # expected section not found; process the section found
                            try:
                                n_section = sections.index(line)
                            except ValueError:
                                raise Exception(f'Unknown section id :{line} (line number = {n_line}).')
                            next_sect = n_section + 1
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
        print(f'\nDistribution of non-zeros values:\n{self.mat["abs_val"].describe()}')
        print(f'\nDistribution of log10(values):\n{self.mat["log"].describe()}')

    def report(self):
        df1 = self.mat.loc[self.mat['log'] == -6]
        df2 = self.mat.loc[self.mat['log'] < -6]
        df3 = self.mat.loc[self.mat['log'] == -7]
        df4 = self.mat.loc[self.mat['log'] == -8]
        df5 = self.mat.loc[self.mat['log'] == -9]
        df6 = self.mat.loc[self.mat['log'] == -10]
        df7 = self.mat.loc[self.mat['log'] > 4]
        print(f'\nDistribution of small values of log10(values) < -6:\n{df2["log"].describe()}')
        print(f'\nNumber of small values of log10(values) == -6: {df1["log"].count()}')
        print(f'Number of small values of log10(values) == -7: {df3["log"].count()}')
        print(f'Number of small values of log10(values) == -8: {df4["log"].count()}')
        print(f'Number of small values of log10(values) == -9: {df5["log"].count()}')
        print(f'Number of small values of log10(values) == -10: {df6["log"].count()}')
        print(f'Number of large values of log10(values) > 4: {df7["log"].count()}')


if __name__ == '__main__':
    tstart = dt.now()
    # print('Started at:', str(tstart))
    wrk_dir = './'
    os.chdir(wrk_dir)

    # small MPSs, for testing the code
    # fn_mps = 'data1/aez.mps'  # 5 matrix elems in a row - processing not implemented yet
    # fn_mps = 'data1/diet.mps'
    # fn_mps = 'data1/jg_korh.mps'
    # fn_mps = 'data1/lotfi.mps'

    # trouble-makers MPSs, large (over 1G) files, not posted to gitHub
    # all large MPS files should be preferably copied to /t/fricko/for_marek/ (see the Oliver's slack post of Jan 9th)
    fn_mps = 'data1/of_led1.mps'    # posted by Oliver in /t/fricko...
    # fn_mps = 'data1/of_baselin.mps'   # second MPS from Oliver

    repdir = 'Rep/'    # subdirectory for reports
    lp = LPdiag(repdir)
    lp.rd_mps(fn_mps)
    lp.report()
    tend = dt.now()
    print('Started at: ', str(tstart))
    print('Finished at:', str(tend))
