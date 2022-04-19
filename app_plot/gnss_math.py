# -*- coding: utf-8 -*- #
# ------------------------------------------------------------------
# File Name:        gnss_math.py
# Author:           hlgou
# Version:          V 1.0
# Created:          2022-04-17
# Description:      Calculation indicator function
# History:
#          2022-04-17   hlgou       :Create the file
#          20xx-xx-xx   xxxxxxxx    :Do some modified
# ------------------------------------------------------------------

import numpy as np
import math

# TODO: TTFF, Conv time, rms for enu
def get_TTFF(Data):
    """
    get the TTFF of the error serial.
    
    > @param[in] Data:          ENU error data, 
    > @param[out] TTFF:         The Time of First Fix
    return:     
    """
    time = Data[:, 0]
    EE = Data[:,1]
    NN = Data[:,2]
    UU = Data[:,3]
    
    return 0

def get_ConTime(Data):
    """
    get the Convergence Time of the error serial.
    
    > @param[in] Data:          ENU error data, 
    > @param[out] TTFF:         The Time of First Fix
    return:     
    """
    time = Data[:, 0]
    EE = Data[:,1]
    NN = Data[:,2]
    UU = Data[:,3]
    beg_time=time[0]
    ind=np.where(pow(EE,2)+pow(NN,2)<0.05, True, False)
    index=find_series(ind)
    ctime=time[index]-beg_time
    return ctime


def find_series(data,ln=10):
    """
    get the first TRUE whose next ln are all TRUE.
    
    > @param[in] data:          the series of bool value
    > @param[in] ln:            the length of we need
    > @param[out] ind:          The index which meet the require
    return:     
    """
    n=len(data)
    pos=0
    if n<ln:
        return n-1
    i=0
    while i<n:
        if data[i]:
            for j in range(i+1,i+ln+1):
                if not data[j]:
                    i=j
                    break
            return j-ln
        i=i+1
    return n-1


def get_FixRate(Data):
    
    return 0

def get_RMS(Data):
    time = Data[:, 0]
    EE = Data[:,1]
    NN = Data[:,2]
    UU = Data[:,3]

    return 0