# Introduction

* Description: This file is for recording how to use the scripts.
* Author: hlgou, WHU
* Creat Date: 2022-01-17

# Sec.1 Run GREAT


# Sec.2 Draw
## 1.Position Error seriel
### a.ENU files
This is about ENU files Reading and plotting
```python
wdir=r'C:\Users\OHanlon\Desktop\plotData\xyz\enu'
pdir=r'C:\Users\OHanlon\Desktop\plotData\xyz\enu\png'
sitelist=['ABMF','ABPO']
time='2022049'
sys='GEC'
freq='2'
comb='IF'
ambfix='AR'
len = 350
lgds = ['Float','2-Fixed','3-Fixed'] 
for site in sitelist:
    f1=wdir+'\\'+ site + '_'+time +'_'+ sys+'_2_IF_F.enu'
    f2=wdir+'\\'+ site + '_'+time +'_'+ sys+'_2_IF_AR.enu'
    f3=wdir+'\\'+ site + '_'+time +'_'+ sys+'_3_IF_AR.enu'
    data1=gr.read_enu(f1)
    data2=gr.read_enu(f2)
    data3=gr.read_enu(f3)
    Data = np.array([data1[0:len,:],data2[0:len,:],data3[0:len,:]])
    tmp=data1[0:len,:]
    print(gm.get_ConTime(tmp))
    gd.draw_enu(Data,lgds,site)
    path=pdir+'\\'+site+'.jpg'
    plt.savefig(path, dpi=300)
```
### a.FLT files

# Sec.3 Utility Function
## 1.Coordinate Transformation
这部分需要的函数在 func/coordinate.py 文件中，使用前需要导入
```python
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from funcs import coordinate as fcd
```
### 1. XYZ2BLH
地心地固三维空间直角坐标系坐标 → 大地坐标系坐标
```python
xyz0=[4097216.5161,4429119.2396,-2065771.1467]
x0=xyz0[0];y0=xyz0[1];z0=xyz0[2]
ell0=fcd.cart2ell(x0,y0,z0)
print(ell0)
b0=ell0[0];l0=ell0[1];h0=ell0[2]
```
### 2. BLH2XYZ

### 3. XYZ2ENU
地心地固三维空间直角坐标系坐标 → 站心坐标系坐标
在进行这个转换的时候需要知道站心的三维空间直角坐标或者大地坐标，下面是已知空间直角坐标xyz0为例进行转换的
```python
xyz0=[4097216.5161,4429119.2396,-2065771.1467]
xyz1=[4097216.4745,4429119.2696,-2065771.1622]
dx=xyz1[0]-xyz0[0];dy=xyz1[1]-xyz0[1];dz=xyz1[2]-xyz0[2]
x0=xyz0[0];y0=xyz0[1];z0=xyz0[2]
ell0=fcd.cart2ell(x0,y0,z0)
enu1=fcd.dxyz2enu(ell0,[dx,dy,dz])      # 方法一，by jqwu
print(enu1)
```
方法二
```python
xyz0=[4097216.5161,4429119.2396,-2065771.1467]
xyz1=[4097216.4745,4429119.2696,-2065771.1622]
dx=xyz1[0]-xyz0[0];dy=xyz1[1]-xyz0[1];dz=xyz1[2]-xyz0[2]
x0=xyz0[0];y0=xyz0[1];z0=xyz0[2]
ell0=fcd.cart2ell(x0,y0,z0)
b0=ell0[0];l0=ell0[1];h0=ell0[2]
rotm=fcd.ell2topo(b0,l0,h0)
dxyz=np.matrix([dx,dy,dz])
print(dxyz*rotm[0],dxyz*rotm[1],dxyz*rotm[2])
```
这两种方法都是得到 ENU，注意，这里顺序可是没有错噢。第一种方法更为简单；当进行很多个坐标转换的时候，也许第二种方法的效率更高一些。