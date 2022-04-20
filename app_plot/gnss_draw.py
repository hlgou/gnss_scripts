# -*- coding: utf-8 -*- #
# ------------------------------------------------------------------
# File Name:        gnss_draw.py
# Author:           hlgou
# Version:          V 1.0
# Created:          2022-04-17
# Description:      plot function
# History:
#          2022-04-17   hlgou       :Create the file
#          20xx-xx-xx   xxxxxxxx    :Do some modified
# ------------------------------------------------------------------
import os
import math
from re import sub
from turtle import title
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from funcs import gnss_time as gt
import gnss_read as gr

import mg_math as mm

def draw_enu(Data, lgds={}, title=["ENU data"]):
    """
    plot enu.
    
    > @param[in] Data:          ENU data, 
    > @param[in] lgds:          the legends, stand for the number of ENU type
    > @param[in] title:         the titile of the figure
    return:     0
    """
    nline = Data.shape[0]  # Data: nlines x 4 x nepos
    nlgds = len(lgds)
    if nline != nlgds:
        print('ERROR::[gnss_draw::draw_enu] the dims of lgds != the dims of data')
        return
    col = 'rbgcmyk'
    cl = col[0:nline]
    cl = cl[::-1]  # the color to plot
    data = Data[0, :, :]
    time = data[:, 0]
    size = len(time)
    time = time / 3600
    shour = time[0] - 0.1  # start time(hour)
    ehour = time[-1] + 0.1  # end time(hour)
    x_ticks = np.arange(math.ceil(shour), math.floor(ehour) + 1, 1)
    # begin plot
    fig = plt.figure()
    fig.patch.set_alpha(0.5)
    font1 = {'weight': 600, 'size': 15}
    ax1 = plt.subplot(3, 1, 1)
    plt.title(title, font1)
    ax1.set_ylabel('East(m)', font1)
    ax1.set_xlim(shour, ehour)
    ax1.set_ylim(-0.1, 0.1)
    plt.xticks(x_ticks, [])
    pAll = []
    all_rms = []
    x_rms=[];y_rms=[];z_rms=[]
    for i in range(0, nline):
        p1 = ax1.scatter(time, Data[i, :, 1], s=5, color=cl[i], marker='.')
        tmp=mm.rms(Data[i, :, 1])
        x_rms.append(tmp)
        pAll.append(p1)
    plt.legend(pAll, lgds, ncol=nline, frameon=False)  # legend::a little bug
    plt.grid(axis="y")

    ax1 = plt.subplot(3, 1, 2)
    ax1.set_ylabel('North(m)', font1)
    ax1.set_ylim(-0.1, 0.1)
    ax1.set_xlim(shour, ehour)
    plt.xticks(x_ticks, [])
    for i in range(0, nline):
        ax1.scatter(time, Data[i, :, 2], s=5, color=cl[i], marker='.')
        tmp=mm.rms(Data[i, :, 2])
        y_rms.append(tmp)
    plt.grid(axis="y")

    ax1 = plt.subplot(3, 1, 3)
    ax1.set_ylabel('Up(m)', font1)
    ax1.set_ylim(-0.2, 0.2)
    ax1.set_xlim(shour, ehour)
    plt.xticks(x_ticks)
    for i in range(0, nline):
        ax1.scatter(time, Data[i, :, 3], s=5, color=cl[i], marker='.')
        tmp=mm.rms(Data[i, :, 3])
        z_rms.append(tmp)
    plt.grid(axis="y")
    plt.xlabel('GPS Time (Hour of day)', font1)
    all_rms=[x_rms,y_rms,z_rms]
    return all_rms



def plot_sat_pdop(time=[], status=[], pct=[]):
    """
    plot the Number of visible satellites, pdop and ambiguity fix rate.
    
    > @param[in] time:          The time data, seconds of week
    > @param[in] status:        The main data, contain nsat, pdop, fixed or float
    > @param[in] pct:           The fix rate
    return:     void
    """
    # RMS
    size = len(time)
    #hour=(time-time[0])/3600
    hour = time / 3600
    nsat = status[:, 0]
    pdop = status[:, 1]
    amb = status[:, 2]
    fig, ax = plt.subplots(1, 1, figsize=(12, 6), dpi=100, facecolor='w')
    ax.patch.set_alpha(0.0)
    #ax.grid(linestyle='--',linewidth=0.3, color='blue',axis='both')
    # Plot NSAT
    ax.scatter(hour, nsat, s=15, color='#0804f9', marker='s', zorder=50, label='NSAT')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])
    ax.set_xlim(hour[0] - 0.01, hour[size - 1] + 0.01)
    ax.set_ylim(4, max(nsat) + 1)
    font1 = {'weight': 600, 'size': 15}
    ax.set_xlabel('GPS Time(hour of week)', font1)
    ax.set_ylabel('Number of Satellite', font1)
    ax.tick_params(axis='both', colors='black', direction='out', labelsize=15, width=1, length=1, pad=5)
    # Plot PDOP
    ax1 = ax.twinx()
    ax1.set_ylabel('PDOP', font1)
    ax1.scatter(hour, pdop, s=15, color='#cf1b1b', marker='o', zorder=1, label='PDOP')
    ax1.tick_params(axis='both', colors='black', direction='out', labelsize=15, width=1, length=1, pad=5)
    ax1.set_ylim(0.5, max(pdop) + 1)
    # Plot AMB
    ax3 = ax.twinx()
    [x, y] = [[hour[i] for i in range(size) if amb[i] == 0], [0.99 for i in range(size) if amb[i] == 0]]
    ax3.scatter(x, y, s=50, color='red', marker='|')
    ax3.scatter([], [], s=20, color='red', marker='s', label='Float')
    ax3.set_yticks(np.arange(0, 1.1, 0.1))
    ax3.set_yticks([])
    ax2 = ax.twinx()
    [x, y] = [[hour[i] for i in range(size) if amb[i] == 1], [0.99 for i in range(size) if amb[i] == 1]]
    ax2.scatter(x, y, s=50, color='lime', marker='|')
    ax2.scatter([], [], s=20, color='lime', marker='s', label='Fixed')
    ax2.set_yticks(np.arange(0, 1.1, 0.1))
    ax2.set_yticks([])
    # Legend
    font2 = {'weight': 600, 'size': 15}
    legend = fig.legend(bbox_to_anchor=(0.585, 0.2),
                        loc='upper right',
                        prop=font2,
                        framealpha=0.0,
                        facecolor='none',
                        ncol=5,
                        numpoints=5,
                        markerscale=2,
                        handlelength=1)
    s = 'Fixing Rate:{0:.2f}%'.format(pct * 100)
    ax.text(0.745, 0.04, s, bbox=dict(facecolor='none', alpha=0.0, pad=6), fontdict=font2, transform=ax.transAxes)


def plot_flt(time=[], data=[], title='FLT'):
    """
    plot flt.
    
    > @param[in] time:          The time data, seconds of week
    > @param[in] data:          The xyz data
    > @param[in] title:         The titile of the figure
    return:     void
    """
    time = time / 3600
    shour = time[0] - 0.1
    ehour = time[-1] + 0.1
    x_ticks = np.arange(math.ceil(shour), math.floor(ehour) + 1, 1)
    # begin plot
    font1 = {'weight': 600, 'size': 15}
    fig, ax = plt.subplots(1, 1, figsize=(12, 6), dpi=100, facecolor='w')
    plt.title(title, font1)
    ax.set_ylabel('position error(m)', font1)
    #ax1.set_xlim(shour, ehour)
    yl=0.5
    ax.set_ylim(-yl,yl)
    # plt.xticks(x_ticks)
    point_size=15
    p1 = ax.scatter(time, data[:, 0], s=point_size, color='r', marker='.')
    p2 = ax.scatter(time, data[:, 1], s=point_size, color='g', marker='.')
    p3 = ax.scatter(time, data[:, 2], s=point_size, color='b', marker='.')
    lgds = ['E', 'N', 'U']
    plt.legend([p1, p2, p3], lgds, prop=font1, ncol=3, frameon=False)
    plt.grid(axis="y")
    s = mm.rms2str(data, 1)
    ax.text(0.01, 0.94, s, bbox=dict(facecolor='none', alpha=0.0, pad=6), fontdict=font1, transform=ax.transAxes)
    ax.set_xlabel('GPST(Hour of week)', font1)
    ax.tick_params(axis='both', colors='black', direction='out', labelsize=15, width=1, length=1, pad=5)
   

def draw_wl(Data, sys='G', mode='WL'):
    """
    plot sigle daily upd.
    
    > @param[in] Data:          The upd data, the first row is the prn, the second is the value
    > @param[in] sys:           The system ID
    > @param[in] mode:          The UPD mode
    return:     void
    """
    font1 = {'weight': 60, 'size': 10}
    sats = Data[0]
    nsat = len(sats)
    satlist = num2str_char(sys, sats)
    # plt.figure()
    fig, ax = plt.subplots(1, 1, dpi=100, facecolor='w')
    plt.bar(list(range(1, nsat + 1)), Data[1], color='b')
    ticks = list(range(1, nsat + 1))
    mean_val, std_val=mm.rmse(Data[1])
    plt.xticks(ticks, satlist, rotation=60, fontsize=8)
    plt.xlabel('PRN')
    plt.ylabel('UPD Value(Cycle)')
    #plt.ylim(0, 0.5)
    plt.grid(axis="y")
    s = 'Mean: {0:.3f}\n  STD: {1:.3f}'.format(mean_val,std_val)
    plt.text(0.01, 0.8, s, bbox=dict(facecolor='none', alpha=0.0, pad=6), fontdict=font1,transform=ax.transAxes)
    plt.title(mode)
    return 0


def draw_wl_line(Data,sys='G',tit='WL'):
    sats = Data[0][1:]
    result=[sats]
    rms=[]
    nsat = len(sats)
    satlist = num2str_char(sys, sats)
    nepo = len(Data) - 1
    time = []
    font1 = {'weight': 60, 'size': 10}
    for i in range(1, nepo + 1):
        doy = Data[i][0]
        time.append(doy)
    c = plt.cm.get_cmap('jet', nsat)
    # fig = plt.figure()
    ref_sat = nsat-6             # the reference satellite
    y = list(map(list, zip(*Data)))[ref_sat][1:]
    y0 = slip_wl_upd(y)
    for i in range(1, nsat + 1):
        x = time
        y = list(map(list, zip(*Data)))[i][1:]
        y = slip_wl_upd(y)
        y=[y[i]-y0[i] for i in range(len(y))]
        mean_val, std_val=mm.rmse(y)
        rms.append(std_val)
        plt.plot(x, y, 'rv-', color=c(i - 2))
    plt.legend(satlist, ncol=8, frameon=False, framealpha=False, fontsize='x-small')
    plt.xlabel('DOY', font1)
    plt.ylabel('UPD(Cycles)', font1)
    plt.ylim(-1, 1.5)
    ticks = list(range(-2, 3))
    ticks = [i * 0.5 for i in ticks]
    plt.yticks(ticks, fontsize=10)
    plt.title(S2System(sys)+' '+tit)
    result.append(rms)
    return result


def slip_wl_upd(y):
    z=y
    for i in range(len(y)):
        tmp=y[i]
        if i==0:
            if tmp>1:
                z[i]=tmp-1
            if tmp<0:
                z[i]=tmp+1
            continue
        if tmp-z[i-1]>0.51:
            z[i]=tmp-1
        if z[i-1]-tmp>0.51:
            z[i]=tmp+1  
    return z

def draw_nl(Data, sys, tit='NL'):
    """
    plot epoch upd.
    
    > @param[in] Data:          The upd data, the first row is the prn, the first two cloums is time
    > @param[in] sys:           The system ID
    > @param[in] tit:           The title of the figure
    return:     void
    """
    rms=[]
    sats = Data[0][2:]
    nsat = len(sats)
    satlist = num2str_char(sys, sats)
    mjd = Data[1][0]
    sod = Data[1][1]
    week0, how0 = gt.sod2how(mjd, sod)
    nepo = len(Data) - 1
    time = []
    font1 = {'weight': 60, 'size': 10}
    for i in range(1, nepo + 1):
        mjd = Data[i][0]
        sod = Data[i][1]
        week, how = gt.sod2how(mjd, sod)
        if week > week0:
            how = how + (week - week0) * 24 * 7
        time.append([week0, how])
    c = plt.cm.get_cmap('jet', nsat)
    # fig = plt.figure()
    for i in range(2, nsat + 2):
        x = list(map(list, zip(*time)))[1]
        y = list(map(list, zip(*Data)))[i][1:]
        mean_val, std_val=mm.rmse(y)
        rms.append(std_val)
        plt.scatter(x, y, s=5, color=c(i - 2), marker='*')
    plt.legend(satlist, ncol=8, frameon=False, framealpha=False, fontsize='x-small')
    plt.xlabel('GPS Time (Hour of week)', font1)
    plt.ylabel('UPD(Cycles)', font1)
    plt.ylim(-1, 1)
    ticks = list(range(-2, 3))
    ticks = [i * 0.5 for i in ticks]
    plt.yticks(ticks, fontsize=10)
    plt.title(S2System(sys) + ' in GPS Week: ' + '%04d' % week0)
    mean_val, std_val=mm.rmse(rms)
    return mean_val


def draw_nl_sec(Data, sys, beg, end):
    """
    plot epoch upd in a specific time period.
    
    > @param[in] Data:          The upd data, the first row is the prn, the first two cloums is time
    > @param[in] sys:           The system ID
    > @param[in] beg:           The begin time(seconds of week)
    > @param[in] end:           The end time(seconds of week)
    return:     void
    """
    mode = 'NL'
    sats = Data[0][2:]
    nsat = len(sats)
    satlist = num2str_char(sys, sats)
    mjd = Data[1][0]
    sod = Data[1][1]
    week0, how0 = gt.sod2sow(mjd, sod)
    nepo = len(Data) - 1
    time = []
    font1 = {'weight': 60, 'size': 10}
    for i in range(1, nepo + 1):
        mjd = Data[i][0]
        sod = Data[i][1]
        week, how = gt.sod2how(mjd, sod)
        if week > week0:
            how = how + (week - week0) * 24 * 7 * 3600
        if how == beg:
            sindex = i
        if how == end:
            eindex = i
        if how >= beg and how <= end:
            time.append([week0, how])
    c = plt.cm.get_cmap('jet', nsat)
    # fig = plt.figure()
    for i in range(2, nsat + 2):
        x = list(map(list, zip(*time)))[1]
        y = list(map(list, zip(*Data)))[i][sindex:eindex]
        plt.scatter(x, y, s=5, color=c(i - 2), marker='*')
    plt.legend(satlist, ncol=8, frameon=False, framealpha=False, fontsize='x-small')
    plt.xlabel('GPS Time (Hour of week)', font1)
    plt.ylabel(mode + ' UPD(Cycles)', font1)
    plt.ylim(-1, 1)
    ticks = list(range(-2, 3))
    ticks = [i * 0.5 for i in ticks]
    plt.yticks(ticks, fontsize=10)
    # plt.title(S2System(sys) + ' in GPS Week: ' + '%04d' % week0)


def draw_nlSigle(Data, sys, mode='NL'):
    """
    plot sigle upd, subFunc of draw_updAll.
    
    > @param[in] Data:          The upd data, the first row is the prn, the first two cloums is time
    > @param[in] sys:           The system ID
    > @param[in] mode:          The UPD mode
    return:     void
    """
    sats = Data[0][2:]
    nsat = len(sats)
    satlist = num2str_char(sys, sats)
    mjd = Data[1][0]
    sod = Data[1][1]
    week0, how0 = gt.sod2how(mjd, sod)
    nepo = len(Data) - 1
    time = []
    font1 = {'weight': 60, 'size': 10}
    for i in range(1, nepo + 1):
        mjd = Data[i][0]
        sod = Data[i][1]
        week, how = gt.sod2how(mjd, sod)
        if week > week0:
            how = how + (week - week0) * 24 * 7
        time.append([week0, how])
    c = plt.cm.get_cmap('jet', nsat)
    for i in range(2, nsat + 2):
        x = list(map(list, zip(*time)))[1]
        y = list(map(list, zip(*Data)))[i][1:]
        plt.scatter(x, y, s=5, color=c(i - 2), marker='*')
    plt.legend(satlist, ncol=8, frameon=False, framealpha=False, fontsize='x-small')
    #plt.xlabel('GPS Time (Hour of week)',font1)
    plt.ylabel(mode + ' UPD(Cycles)', font1)
    ticks = list(range(-2, 3))
    ticks = [i * 0.5 for i in ticks]
    plt.yticks(ticks, fontsize=10)
    plt.ylim(-1, 1)


def draw_nlSigle_sec(Data, sys, beg, end, mode='NL'):
    """
    plot sigle upd in a specific time period, subFunc of draw_updAll_sec.
    
    > @param[in] Data:          The upd data, the first row is the prn, the first two cloums is time
    > @param[in] sys:           The system ID
    > @param[in] beg:           The begin time(seconds of week)
    > @param[in] end:           The end time(seconds of week)
    > @param[in] mode:          The UPD mode
    return:     void
    """
    sats = Data[0][2:]
    nsat = len(sats)
    satlist = num2str_char(sys, sats)
    mjd = Data[1][0]
    sod = Data[1][1]
    week0, how0 = gt.sod2sow(mjd, sod)
    nepo = len(Data) - 1
    time = []
    sindex = 0
    eindex = 0
    font1 = {'weight': 60, 'size': 10}
    for i in range(1, nepo + 1):
        mjd = Data[i][0]
        sod = Data[i][1]
        week, how = gt.sod2sow(mjd, sod)
        if week > week0:
            how = how + (week - week0) * 24 * 7 * 3600
        if how == beg:
            sindex = i
        if how == end:
            eindex = i
        if how >= beg and how <= end:
            time.append([week0, how])
    c = plt.cm.get_cmap('jet', nsat)
    for i in range(2, nsat + 2):
        x = list(map(list, zip(*time)))[1]
        if beg > 3600 * 24 * 7:
            x = [i - 3600 * 24 * 7 for i in x]
        y = list(map(list, zip(*Data)))[i][sindex:eindex + 1]
        plt.scatter(x, y, s=5, color=c(i - 2), marker='*')
    plt.legend(satlist, ncol=8, frameon=False, framealpha=False, fontsize='x-small')
    #plt.xlabel('GPS Time (Hour of week)',font1)
    plt.ylabel(mode + ' UPD(Cycles)', font1)
    ticks = list(range(-2, 3))
    ticks = [i * 0.5 for i in ticks]
    plt.yticks(ticks, fontsize=10)
    plt.ylim(-1, 1)


def draw_nlAll(wl, nl, sysAll='GECR'):
    """
    plot all upd, the main plot function.
    
    > @param[in] wl:            The first epoch upd data set
    > @param[in] nl:            The second epoch upd data set
    > @param[in] sysAll:        The system ID
    return:     void
    """
    num = len(sysAll)
    for i in range(0, num):
        sys = sysAll[i]
        if sys == 'G':
            j = 0
        elif sys == 'E':
            j = 1
        elif sys == 'C':
            j = 2
        elif sys == 'R':
            j = 3
        fig = plt.figure()
        plt.subplot(2, 1, 1)
        mjd = wl[0][1][0]
        week, day = gt.mjd2gpsweek(mjd)
        plt.title(S2System(sys) + ' in GPS week: ' + '%04d' % week)
        draw_nlSigle(wl[j], sys, 'WL')
        plt.subplot(2, 1, 2)
        draw_nlSigle(nl[j], sys)
        font1 = {'weight': 60, 'size': 10}
        plt.xlabel('GPS Time (Hour of week)', font1)


def draw_nlminsAll(wl, nl, sysAll='GEC'):
    """
    plot all mins upd, the main plot function.
    
    > @param[in] wl:            The first epoch upd data set
    > @param[in] nl:            The second epoch upd data set
    > @param[in] sysAll:        The system ID
    return:     void
    """
    num = len(sysAll)
    for i in range(0, num):
        sys = sysAll[i]
        if sys == 'G':
            j = 0
        elif sys == 'E':
            j = 1
        elif sys == 'C':
            j = 2
        fig = plt.figure()
        nl1, nl2 = gr.preminusNL(wl[j], nl[j])
        wl1 = gr.minusNL(nl1, nl2)
        draw_nl(wl1, sys)
    return 0


def draw_nlmins(wl, nl, sys='G'):
    """
    plot single system mins upd.
    
    > @param[in] wl:            The first epoch upd data
    > @param[in] nl:            The second epoch upd data
    > @param[in] sys:           The system ID
    return:     void
    """
    nl1, nl2 = gr.preminusNL(wl, nl)
    wl1 = gr.minusNL(nl1, nl2)
    draw_nl(wl1, sys)
    return 0


def draw_nlAll_sec(wl, nl, sysAll, beg, end):
    """
    plot all upd in a specific time period, the main plot function.
    
    > @param[in] wl:            The first epoch upd data
    > @param[in] nl:            The second epoch upd data
    > @param[in] sysAll:        The system ID
    > @param[in] beg:           The begin time(seconds of week)
    > @param[in] end:           The end time(seconds of week)
    return:     void
    """
    num = len(sysAll)
    for i in range(0, num):
        sys = sysAll[i]
        fig = plt.figure()
        plt.subplot(2, 1, 1)
        mjd = wl[0][1][0]
        week, day = gt.mjd2gpsweek(mjd)
        plt.title(S2System(sys) + ' in GPS week: ' + '%04d' % week)
        draw_nlSigle_sec(wl[i], sys, beg, end, 'WL')
        plt.subplot(2, 1, 2)
        draw_nlSigle_sec(nl[i], sys, beg, end)
        font1 = {'weight': 60, 'size': 10}
        plt.xlabel('GPS Time (Hour of week)', font1)


def num2str_char(sys, num):
    """
    Util: convert sat_num list to sat_char list.
    
    > @param[in] sys:           The the system ID
    > @param[in] num:           The prn number-satlist
    return:     
    < @param[out] satlist:      The char-satlist 
    """
    satlist = []
    for i in num:
        satlist.append(sys + '%02d' % i)
    return satlist


def S2SYS(sys):
    """
    Util: Convert single-Sysid to 3-char sysid.
    
    > @param[in] sys:           The the system ID
    return:     
    < @param[out] SYS:          The 3-char of system 
    """
    if sys == 'G':
        return 'GPS'
    elif sys == 'E':
        return 'GAL'
    elif sys == 'C':
        return 'BDS'
    elif sys == 'R':
        return 'GLO'
    elif sys == 'Q':
        return 'QZS'
    else:
        return 'Undefined sysid'


def S2System(sys):
    """
    Util: Convert single-Sysid to system full name.
    
    > @param[in] sys:           The the system ID
    return:     
    < @param[out] system:       The full-name of system 
    """
    if sys == 'G':
        return 'GPS'
    elif sys == 'E':
        return 'Galileo'
    elif sys == 'C':
        return 'BDS'
    elif sys == 'R':
        return 'GLONASS'    
    elif sys == 'Q':
        return 'QZSS'
    else:
        return 'Undefined sysid'
