# -*- coding: utf-8 -*- #
# ------------------------------------------------------------------
# File Name:        gnss_read.py
# Author:           hlgou
# Version:          V 1.0
# Created:          2022-04-17
# Description:      get the data from the files
# History:
#          2022-04-17   hlgou       :Create the file
#          20xx-xx-xx   xxxxxxxx    :Do some modified
# ------------------------------------------------------------------

import os
import logging
import math
import copy
from numpy.core.records import array
from numpy.lib.function_base import append
import pandas as pd
import numpy as np


def read_enu(f_enu):
    """ 
    get enu data from enu file 
    
    > @param[in] f_enu:         the enu filename
	return: 
	< @param[out] data:         array(T,E,N,U)
    """
    try:
        with open(f_enu) as f:
            lines = f.readlines()
    except FileNotFoundError:
        logging.error(f"file not found {f_enu}")
        return
    sod, de, dn, du = [], [], [], []
    for line in lines:
        info = line.split()
        if len(info) == 4:
            if info[0] == "RMS":
                break
            sod.append(float(info[0]))
            de.append(float(info[1]))
            dn.append(float(info[2]))
            du.append(float(info[3]))
    # minus the value of the last epoch
    if math.fabs(du[-1]) > 1:
        de1 = [i - de[-1] for i in de]
        dn1 = [i - dn[-1] for i in dn]
        du1 = [i - du[-1] for i in du]
        de = de1
        dn = dn1
        du = du1
    return np.array([sod, de, dn, du]).T

