# ------------------------------------------------------------------
# File Name:        mg_math.py
# Author:           hlgou
# Version:          V 1.0
# Created:          2021-10-15
# Description:      some small functions about math
# History:
#          2021-10-15   hlgou       :Create the file
#          2021-xx-xx   xxxxx       :add some functions
# ------------------------------------------------------------------

import numpy as np
import math

def get_lcm(a, b):
    """ get the least common multiple of two integers """
    for i in range(min(a, b), 0, -1):
        if a % i == 0 and b % i == 0:
            return a * b // i


def get_digNum(dig):
    """ get the number of decimal places """
    s = str(dig)
    info = s.split('.')
    if len(info) < 2:
        return 0
    return len(info[1])


def get_diglcm(a, b):
    """ get the least common multiple of two decimals """
    times = 10**max(get_digNum(a), get_digNum(b))
    lcm = get_lcm(round(a * times), round(b * times))
    return lcm / times


def rms(diff=[]):
    """ get the rms of the data """
    n= len(diff)        # delete NAN
    tmp=[]
    for i in range(n):
        if not (math.isnan(diff[i])):
            tmp.append(diff[i])
    diff=tmp
    sum = 0
    for i in range(n):
        sum = sum + diff[i] * diff[i]
    return math.sqrt(sum / (n))


def rms2str(diff=[], mode=0):
    """ get the rms of the data """
    r1 = rms(diff[:, 0])
    r2 = rms(diff[:, 1])
    r3 = rms(diff[:, 2])
    if mode:
        a = 3
        r13 = r1 * a
        r23 = r2 * a
        r33 = r3 * a
        size = len(diff)
        sum1 = 0
        sum2 = 0
        sum3 = 0
        n3 = 0
        for i in range(size):
            if abs(diff[i, 0]) >= r13 or abs(diff[i, 1]) >= r23 or abs(diff[i, 2]) >= r33:
                n3 = n3 + 1
                continue
            else:
                sum1 = sum1 + diff[i, 0] * diff[i, 0]
                sum2 = sum2 + diff[i, 1] * diff[i, 1]
                sum3 = sum3 + diff[i, 2] * diff[i, 2]
        # print(n3)
        r1 = math.sqrt(sum1 / (size - n3))
        r2 = math.sqrt(sum2 / (size - n3))
        r3 = math.sqrt(sum3 / (size - n3))
    s = 'RMS:({0:.4f},{1:.4f},{2:.4f})m'.format(r1, r2, r3)
    return s


def rmse(data=[], mode=0, big_error=2):
    """
    get the rmse(std) and the mean value of the data.
    
    > @param[in] data:          the data list
    > @param[in] mode:          mode=1: Eliminate noise
    > @param[in] big_error:     the times of rms
    return: 
    < @param[out] mean_val:     the mean value
    < @param[out] std_val:      the std value
    """
    n= len(data)        # delete NAN
    tmp=[]
    for i in range(n):
        if not (math.isnan(data[i])):
            tmp.append(data[i])
    data=tmp
    mean_val = np.mean(data)
    error = [i - mean_val for i in data]
    n=len(error)
    data1=[]
    for i in range(n):              #remove big errors
        if abs(error[i])<=big_error:
            data1.append(data[i])
    data=data1

    mean_val = np.mean(data)
    error = [i - mean_val for i in data]
    std_val = np.std(error)
    while mode:
        a = 3
        data1 = []
        for i in range(len(error)):
            if abs(error[i]) < a * std_val:
                data1.append(data[i])
        mean_val = np.mean(data1)
        error = [i - mean_val for i in data1]
        std_val = np.std(error)
        if len(data1) == len(data):
            break
        data = data1
    return mean_val, std_val


def rm_pos_error(time,pos,mode=0):
    time1,pos1,x,y,z=[],[],[],[],[]
    pos1=[]
    nepo=len(time)
    [x_mean_val,x_sigma]=rmse(pos[:,0],mode)
    [y_mean_val,y_sigma]=rmse(pos[:,1],mode)
    [z_mean_val,z_sigma]=rmse(pos[:,2],mode)
    for i in range(nepo):
        if np.abs(pos[i,0]-x_mean_val)<=3*x_sigma and np.abs(pos[i,1]-y_mean_val)<=3*y_sigma and np.abs(pos[i,2]-z_mean_val)<=3*z_sigma:
            x.append(pos[i,0])
            y.append(pos[i,1])
            z.append(pos[i,2])
            time1.append(time[i])
    pos1 = np.array([x, y, z]).T
    time1=np.array(time1).T
    return time1,pos1
            