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
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from funcs import coordinate as fcd
import mg_math as mm


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


def read_enuflt(f_flt, tag=0):
    """ 
    get xyz from flt_file which from PPPLSQ 
    
    > @param[in] f_flt:         the flt filename
    > @param[in] tag:           the result id, 0:fixed+float; 1:only fixed
	return: 
	< @param[out] time:         the time data, second of week
    < @param[out] pos:          the position data, XYZ
    < @param[out] status:       nsat, pdop, state(fixed of float) 
    < @param[out] pct:          the fixed rate
    """
    try:
        with open(f_flt) as f:
            lines = f.readlines()
    except FileNotFoundError:
        logging.error(f"file not found {f_flt}")
        return
    t, x, y, z, nsat, pdop, state = [], [], [], [], [], [], []
    float_num, fixed_num = 0, 0
    for line in lines:
        if line[0] == "#" or line[0] == "%":
            continue
        line = line[31:-1]
        info = line.split()
        value = line.split()
        if value[14] == 'Fixed':
            fixed_num = fixed_num + 1
        else:
            float_num = float_num + 1
            if tag != 0:
                continue
        t.append(float(value[0]))
        x.append(float(value[2]))
        y.append(float(value[3]))
        z.append(float(value[4]))
        nsat.append(float(value[10]))
        if value[13] == '-nan(ind)':
            pdop.append(0.0)
        else:
            pdop.append(float(value[13]))
        if value[14] == 'Fixed':
            sta = 1
        else:
            sta = 0
        state.append(sta)
    time = np.array(t).T
    pos = np.array([x, y, z]).T
    status = np.array([nsat, pdop, state]).T
    pct = fixed_num / (fixed_num + float_num)
    return (time, pos, status, pct)


def read_flt(f_flt, tag=0, beg=0, end=999999, interval=1):
    """ 
    get xyz data from flt file which from PVTFLT 
    
    > @param[in] f_flt:         the flt filename
    > @param[in] tag:           the result id, 0:fixed+float; 1:only fixed
    > @param[in] beg:           the begin time of data
    > @param[in] end:           the end time of data
    > @param[in] interval:      the interval of data
	return: 
	< @param[out] time:         the time data, second of week
    < @param[out] pos:          the position data, XYZ
    < @param[out] status:       nsat, pdop, state(fixed of float) 
    < @param[out] pct:          the fixed rate
    """
    try:
        with open(f_flt) as f:
            lines = f.readlines()
    except FileNotFoundError:
        logging.error(f"file not found {f_flt}")
        return
    t, x, y, z, vx, vy, vz, nsat, pdop, state = [], [], [], [], [], [], [], [], [], []
    float_num, fixed_num = 0, 0
    for line in lines:
        if line[0] == "#" or line[0] == "%":
            continue
        value = line.split()
        if len(value)==0:
            continue
        tmpt=float(value[0])
        if tmpt < beg or tmpt > end or tmpt % interval != 0:
            continue       
        if value[16] == 'Fixed':
            fixed_num = fixed_num + 1
        else:
            float_num = float_num + 1
            if tag != 0:
                continue
        t.append(float(value[0]))
        x.append(float(value[1]))
        y.append(float(value[2]))
        z.append(float(value[3]))
        vx.append(float(value[4]))
        vy.append(float(value[5]))
        vz.append(float(value[6]))
        nsat.append(float(value[13]))
        pdop.append(float(value[14]))
        if value[16] == 'Fixed':
            sta = 1
        else:
            sta = 0
        state.append(sta)
    time = np.array(t).T
    pos = np.array([x, y, z]).T
    vel = np.array([vx, vy, vz]).T
    status = np.array([nsat, pdop, state]).T
    pct = fixed_num / (fixed_num + float_num)
    #print(fixed_num + float_num)
    return (time, pos, status, pct)


def read_enuflt(f_flt, tag=0):
    """ 
    get xyz from flt_file which from PPPLSQ 
    
    > @param[in] f_flt:         the flt filename
    > @param[in] tag:           the result id, 0:fixed+float; 1:only fixed
	return: 
	< @param[out] time:         the time data, second of week
    < @param[out] pos:          the position data, XYZ
    < @param[out] status:       nsat, pdop, state(fixed of float) 
    < @param[out] pct:          the fixed rate
    """
    try:
        with open(f_flt) as f:
            lines = f.readlines()
    except FileNotFoundError:
        logging.error(f"file not found {f_flt}")
        return
    t, x, y, z, nsat, pdop, state = [], [], [], [], [], [], []
    float_num, fixed_num = 0, 0
    for line in lines:
        if line[0] == "#" or line[0] == "%":
            continue
        line = line[31:-1]
        info = line.split()
        value = line.split()
        if value[14] == 'Fixed':
            fixed_num = fixed_num + 1
        else:
            float_num = float_num + 1
            if tag != 0:
                continue
        t.append(float(value[0]))
        x.append(float(value[2]))
        y.append(float(value[3]))
        z.append(float(value[4]))
        nsat.append(float(value[10]))
        if value[13] == '-nan(ind)':
            pdop.append(0.0)
        else:
            pdop.append(float(value[13]))
        if value[14] == 'Fixed':
            sta = 1
        else:
            sta = 0
        state.append(sta)
    time = np.array(t).T
    pos = np.array([x, y, z]).T
    status = np.array([nsat, pdop, state]).T
    pct = fixed_num / (fixed_num + float_num)
    return (time, pos, status, pct)


def read_nl(f_nl, sys='GECR'):
    """
    get epoch upd data from epoch upd file.
    
    > @param[in] f_nl:          filename
    > @param[in] sys:           the system which you want get the data
    return:     
    < @param[out] Gupd:         GPS UPD data
    < @param[out] Eupd:         GAL UPD data
    < @param[out] Cupd:         BDS UPD data
    """
    try:
        with open(f_nl) as f:
            lines = f.readlines()
    except FileNotFoundError:
        logging.error(f"file not found {f_nl}")
        return

    if 'G' in sys:
        G_sats = list(range(1, 33))
    if 'E' in sys:
        E_sats = list(range(1, 37))
    if 'C' in sys:
        C_sats = list(range(1, 62))
    if 'R' in sys:
        R_sats = list(range(1, 26))
    Gn = len(G_sats)
    En = len(E_sats)
    Cn = len(C_sats)
    Rn = len(R_sats)
    Gupd = [[0, 0] + G_sats]
    Eupd = [[0, 0] + E_sats]
    Cupd = [[0, 0] + C_sats]
    Rupd = [[0, 0] + R_sats]
    tag = 0
    for line in lines:
        info = line.split()
        if info[0] == '%' or info[0] == 'EOF':
            continue
        if info[0] == 'EPOCH-TIME':
            if tag == 1:
                Gupd.append(Gtmp)
                Eupd.append(Etmp)
                Cupd.append(Ctmp)
                Rupd.append(Rtmp)
            Gtmp = [float(info[1]), float(info[2])] + list(np.full(Gn, np.nan))
            Etmp = [float(info[1]), float(info[2])] + list(np.full(En, np.nan))
            Ctmp = [float(info[1]), float(info[2])] + list(np.full(Cn, np.nan))
            Rtmp = [float(info[1]), float(info[2])] + list(np.full(Rn, np.nan))
            tag = 1
            continue
        if info[0][0] == 'x':
            continue
        i = int(info[0][1:3])
        if info[0][0] == 'G':
            Gtmp[i + 1] = float(info[1])
        elif info[0][0] == 'E':
            Etmp[i + 1] = float(info[1])
        elif info[0][0] == 'C':
            Ctmp[i + 1] = float(info[1])
        elif info[0][0] == 'R':
            Rtmp[i + 1] = float(info[1]) 
    return delet_nan(Gupd), delet_nan(Eupd), delet_nan(Cupd), delet_nan(Rupd)


def read_wl(f_wl, sys='GEC'):
    """
    get wl upd value from wl upd file.
    
    > @param[in] f_wl:          filename
    > @param[in] sys:           the system which you want get the data
    return:     
    < @param[out] Gupd:         GPS UPD data
    < @param[out] Eupd:         GAL UPD data
    < @param[out] Cupd:         BDS UPD data
    % @note[data]    data structure: the first row is prn, the second is the value
    """
    try:
        with open(f_wl) as f:
            lines = f.readlines()
    except FileNotFoundError:
        logging.error(f"file not found {f_wl}")
        return
    if 'G' in sys:
        G_sats = list(range(1, 33))
    if 'E' in sys:
        E_sats = list(range(1, 37))
    if 'C' in sys:
        C_sats = list(range(1, 62))
    Gn = len(G_sats)
    En = len(E_sats)
    Cn = len(C_sats)
    Gupd = [G_sats]
    Eupd = [E_sats]
    Cupd = [C_sats]
    Gtmp = list(np.full(Gn, np.nan))
    Etmp = list(np.full(En, np.nan))
    Ctmp = list(np.full(Cn, np.nan))
    for line in lines:
        info = line.split()
        if info[0] == '%' or info[0] == 'EOF' or info[0][0] == 'x':
            continue
        i = int(info[0][1:3])
        if info[0][0] == 'G':
            Gtmp[i - 1] = float(info[1])
        elif info[0][0] == 'E':
            Etmp[i - 1] = float(info[1])
        elif info[0][0] == 'C':
            Ctmp[i - 1] = float(info[1])
    Gupd.append(Gtmp)
    Eupd.append(Etmp)
    Cupd.append(Ctmp)
    return delet_wl_nan(Gupd), delet_wl_nan(Eupd), delet_wl_nan(Cupd)


def read_wl_many(path,beg,end,year=2022,sys='GEC'):
    """
    get wl-upd value from wl upd files. You need to CHANGE the file tag!!
    
    > @param[in] path:          the path of the files
    > @param[in] beg:           the begin doy of the files
    > @param[in] end:           the end doy of the files
    > @param[in] year:          the year of the files
    > @param[in] sys:           the system which you want get the data
    return:     
    < @param[out] Gupd:         GPS UPD data
    < @param[out] Eupd:         GAL UPD data
    < @param[out] Cupd:         BDS UPD data
    % @note[data]    data structure: the first row is prn, the first cloumn is the DOY
    """
    if 'G' in sys:
        G_sats = list(range(0, 33))
    if 'E' in sys:
        E_sats = list(range(0, 37))
    if 'C' in sys:
        C_sats = list(range(0, 62))
    Gn = len(G_sats)
    En = len(E_sats)
    Cn = len(C_sats)
    Gupd = [G_sats]
    Eupd = [E_sats]
    Cupd = [C_sats]
    for j in range(beg,end+1):
        f1=path+'\\upd_wl_'+str(year)+str(j).zfill(3)+'_GREC'
        try:
            with open(f1) as f:
                lines = f.readlines()
        except FileNotFoundError:
            logging.error(f"file not found {f1}")
            continue
        Gtmp = list(np.full(Gn, np.nan));Gtmp[0]=j
        Etmp = list(np.full(En, np.nan));Etmp[0]=j
        Ctmp = list(np.full(Cn, np.nan));Ctmp[0]=j
        for line in lines:
            info = line.split()
            if info[0] == '%' or info[0] == 'EOF' or info[0][0] == 'x':
                continue
            i = int(info[0][1:3])
            if info[0][0] == 'G':
                Gtmp[i] = float(info[1])
            elif info[0][0] == 'E':
                Etmp[i] = float(info[1])
            elif info[0][0] == 'C':
                Ctmp[i] = float(info[1])
        Gupd.append(Gtmp)
        Eupd.append(Etmp)
        Cupd.append(Ctmp)
    return delet_wl_nan(Gupd), delet_wl_nan(Eupd), delet_wl_nan(Cupd)


def read_ewl_many(path,beg,end,year=2022,sys='GEC'):
    """
    get ewl-upd value from ewl upd files. You need to CHANGE the file tag!!
    
    > @param[in] path:          the path of the files
    > @param[in] beg:           the begin doy of the files
    > @param[in] end:           the end doy of the files
    > @param[in] year:          the year of the files
    > @param[in] sys:           the system which you want get the data
    return:     
    < @param[out] Gupd:         GPS UPD data
    < @param[out] Eupd:         GAL UPD data
    < @param[out] Cupd:         BDS UPD data
    % @note[data]    data structure: the first row is prn, the first cloumn is the DOY
    """
    if 'G' in sys:
        G_sats = list(range(0, 33))
    if 'E' in sys:
        E_sats = list(range(0, 37))
    if 'C' in sys:
        C_sats = list(range(0, 62))
    Gn = len(G_sats)
    En = len(E_sats)
    Cn = len(C_sats)
    Gupd = [G_sats]
    Eupd = [E_sats]
    Cupd = [C_sats]
    for j in range(beg,end+1):
        f1=path+'\\upd_ewl_'+str(year)+str(j).zfill(3)+'_GEC'
        try:
            with open(f1) as f:
                lines = f.readlines()
        except FileNotFoundError:
            logging.error(f"file not found {f1}")
            continue
        Gtmp = list(np.full(Gn, np.nan));Gtmp[0]=j
        Etmp = list(np.full(En, np.nan));Etmp[0]=j
        Ctmp = list(np.full(Cn, np.nan));Ctmp[0]=j
        for line in lines:
            info = line.split()
            if info[0] == '%' or info[0] == 'EOF' or info[0][0] == 'x':
                continue
            i = int(info[0][1:3])
            if info[0][0] == 'G':
                Gtmp[i] = float(info[1])
            elif info[0][0] == 'E':
                Etmp[i] = float(info[1])
            elif info[0][0] == 'C':
                Ctmp[i] = float(info[1])
        Gupd.append(Gtmp)
        Eupd.append(Etmp)
        Cupd.append(Ctmp)
    return delet_wl_nan(Gupd), delet_wl_nan(Eupd), delet_wl_nan(Cupd)


def delet_nan(data):
    """ Util: delete the row which all the data is nan """
    data = list(map(list, zip(*data)))  # transpose
    data1 = data[:2]
    n = len(data)
    for i in range(2, n):
        tmp = data[i][1:]
        if mIsAllnan(tmp):
            continue
        data1.append(data[i])
    return list(map(list, zip(*data1)))


def delet_wl_nan(data):
    """ Util: delete the row which all the data is nan (WL) """
    data = list(map(list, zip(*data)))  # transpose
    data1 = []
    n = len(data)
    for i in range(0, n):
        tmp = data[i][1:]
        if mIsAllnan(tmp):
            continue
        data1.append(data[i])
    return list(map(list, zip(*data1)))


def mIsAllnan(data):
    """ Util: judge whether all the data is nan """
    n = len(data)
    state = False
    for i in range(n):
        if not (math.isnan(data[i])):
            break
    if i == n - 1 and math.isnan(data[i]):
        state = True
    return state


def dxyzs2ddxyzs(dxyzs,lim=0.2):
    """ 
    Util: Eliminate system error for dxyzs

    > @param[in] dxyzs:         the dxyzs series
    > @param[in] lim:           the threshold to judge whether do the minus operation
    return:     
    < @param[out] ddxyzs:       the result after eliminating system error
    """
    pos = dxyzs.T.tolist()
    dx = pos[0]
    dy = pos[1]
    dz = pos[2]
    if math.fabs(dz[-1]) > lim:
        index=-1
        ddx = [i - dx[index] for i in dx]
        ddy = [i - dy[index] for i in dy]
        ddz = [i - dz[index] for i in dz]
    dx=ddx
    dy=ddy
    dz=ddz
    return np.array([dx, dy, dz]).T


def xyzs2dxyzs_mean(xyzs):
    """ 
    Util: Convert XYZS to DXYZS, All previous epoch data minus the mean value of all data

    > @param[in] xyzs:          the xyz series
    return:     
    < @param[out] dxyzs:        the result after doing the minus operation
    """
    xyzs = xyzs.T.tolist()
    xs = xyzs[0]
    ys = xyzs[1]
    zs = xyzs[2]
    xmean, xstd = mm.rmse(xs, 1)
    ymean, ystd = mm.rmse(ys, 1)
    zmean, zstd = mm.rmse(zs, 1)
    dx = [i - xmean for i in xs]
    dy = [i - ymean for i in ys]
    dz = [i - zmean for i in zs]
    return np.array([dx, dy, dz]).T


def xyzs2dxyzs(xyzs, xyz=[0,0,0]):
    """
    pos-xyz=pos1, for static.
    
    > @param[in] xyzs:          the series of coordinate
    > @param[in] xyz:           the true value of coordinate for static
    return:     
    < @param[out] xyzs:         the crd error, xyxs = xyzs - xyz
    """
    axyz = np.array(xyz)
    n = len(xyzs)
    for i in range(0, n):
        xyzs[i] = xyzs[i] - axyz
    return xyzs


def xyzs2enus(xyzs,xyz):
    """
    Util: Convert XYZS to ENUS for static.
    
    > @param[in] xyzs:          the series of coordinate
    > @param[in] xyz:           the true value of coordinate for static, It can't be zeros!
    return:     
    < @param[out] enus:         the ENUS 
    """
    dxyzs=xyzs2dxyzs(xyzs,xyz)
    n = len(xyzs)
    ell0 = fcd.cart2ell(xyz[0],xyz[1],xyz[2])
    enus=[]
    for i in range(n):
        dx=dxyzs[i][0];dy=dxyzs[i][1];dz=dxyzs[i][2]
        enu1=fcd.dxyz2enu(ell0,[dx,dy,dz])
        enus.append(enu1)
    return np.array(enus)
