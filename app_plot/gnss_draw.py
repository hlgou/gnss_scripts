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

import math
from re import sub
from turtle import title
import numpy as np
import matplotlib.pyplot as plt

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