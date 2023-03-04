"""
Prototype of simple analysis of the MPS-format file
Written by Marek Makowski, ECE Program of IIASA, in March 2023
Developed in PyCharm, with Python 3.10.4
"""

from lpdiag import *    # LPdiag class for processing and analysis of LP matrices
# import sys		# needed for sys.exit()
import os
# import numpy as np
# import pandas as pd
from datetime import datetime as dt
# from datetime import timedelta as td


if __name__ == '__main__':
    """Driver of the LP diagnostics script.
    
    Defines the working space, then constrols the flow by executing the desired functions of LPdiag class.
    """

    tstart = dt.now()
    # print('Started at:', str(tstart))
    wrk_dir = '/Users/marek/Documents/GitHub/marek_iiasa/MCA/LPdiag/'   # should be modified by each user
    os.chdir(wrk_dir)

    # small MPSs, for testing the code, posted to Data/mps_tst dir
    # aez  - agro-ecological zones, medium size; two matrix elems in a row - processing not implemented yet
    # diet - classical small LP
    # jg_korh - tiny testing problem
    # lotfi - classical medium size; two matrix elems in a row - processing not implemented yet

    # trouble-makers MPSs, large (over 1G) files, not posted to gitHub, locally in Data/mps dir
    # all large MPS files should be preferably copied to /t/fricko/for_marek/ (see the Oliver's slack post of Jan 9th)
    # of_led1    # posted by Oliver in /t/fricko... on Jan 9th, 2023
    # of_baselin   # second MPS from Oliver, posted on Jan 9th 2023

    # data_dir = 'Data/mps_tst/'
    # prob_id = 'diet'
    data_dir = 'Data/mps/'
    prob_id = 'of_led1'
    fn_mps = data_dir + prob_id
    repdir = 'Rep_tst/'    # subdirectory for reports

    lp = LPdiag(repdir)     # LPdiag ctor
    lp.rd_mps(fn_mps)       # read MPS, store the matrix in dataFrame
    lp.stat(lo_tail=-7, up_tail=5)   # statistics of the matrix coefficients, incl. distribution tails
    # lp.stat(lo_tail=0, up_tail=0)  # to get numbers of coeffs for each magnitute specify equal/overlapping tails
    lp.out_loc(small=True, thresh=-7)  # locations of outlayers
    tend = dt.now()
    print('\nStarted at: ', str(tstart))
    print('Finished at:', str(tend))
    # TODO: plots of distributions of coeffs
    # TODO: implement naive scaling
    # TODO: distributions of coeff-vals in rows/cols with small/large coefficients
