"""
Prototype of simple analysis of the MPS-format file
Written by Marek Makowski, ECE Program of IIASA, in March 2023
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
    lp.stat()               # basic statistics of the matrix coefficients
    tend = dt.now()
    print('Started at: ', str(tstart))
    print('Finished at:', str(tend))
