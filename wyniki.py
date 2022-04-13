# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from transformacje import *
geo = Transformacje(model = "wgs84")
#import sys
#sys.path.append('C:\\Users\\kinga\\projekt')



plik = "wsp_inp.txt"
# odczyt z pliku: https://docs.scipy.org/doc/numpy-1.15.1/reference/generated/numpy.genfromtxt.html
tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)

p1 = []
l1 = []
h1 = []

# zamiana XYZ na flh
for line in tablica:
    f,l,h = geo.xyz2plh(line[0], line[1], line[2])
    p1.append(f)
    l1.append(l)
    h1.append(h)
    
PLH = np.column_stack((p1,l1,h1))

x1 = []
y1 = []
z1 = []

# zamiana flh na XYZ
for line in PLH:
    X,Y,Z = geo.flh2XYZ(line[0],line[1], line[2])
    x1.append(X)
    y1.append(Y)
    z1.append(Z)

XYZ = np.column_stack((x1,y1,z1))

#przeliczenie do Gaussa-Kruggera i u1992
xgk1 = []
ygk1 = []
for line in PLH:
    xgk92, ygk92 = geo.fl2xy(line[0], line[1], 19)
    xgk1.append(xgk92)
    ygk1.append(ygk92)

XY_GK92 = np.column_stack((xgk1,ygk1))

x92_1 = []
y92_1 = []
for line in XY_GK92:
    x92, y92 = geo.u92(line[0], line[1])
    x92_1.append(x92)
    y92_1.append(y92)

WSP_92 = np.column_stack((x92_1,y92_1))

#przeliczenie do Gaussa-Kruggera i u2000

xgk2 = []
ygk2 = []
for line in PLH:
    xgk2000,ygk2000 = geo.fl2xy(line[0], line[1], 21)
    xgk2.append(xgk2000)
    ygk2.append(ygk2000)

XY_GK2000 = np.column_stack((xgk2,ygk2))

x2000_1 = []
y2000_1 = []
for line in XY_GK2000:
    x2000,y2000 = geo.u2000(line[0], line[1], 21)
    x2000_1.append(x2000)
    y2000_1.append(y2000)

WSP_2000 = np.column_stack((x2000_1,y2000_1))

#obliczenie odległoci pomiędzy kolejnymi punktami
odl1 = []
for line in XYZ:
    x1 = line[0]
    y1 = line[1]
    x2 = (line+1)[0]
    y2 = (line+1)[1]
    dist = geo.dist_xy(x1,y1,x2,y2)
    odl1.append(dist)

az1 = []
ele = []
#obliczenie kąta azymutu i elewacji
for line in PLH:
    fi = line[0]
    lam = line[1]
    h = line[2]
    for line in XYZ:
        x0 = (line+1)[0]
        y0 = (line+1)[1]
        z0 = (line+1)[2]
        azw, el = geo.azymut_elewacja(fi,lam,h,x0,y0,z0)
        az1.append(azw)
        ele.append(el)

for line in XYZ:
    x3 = line[0]
    y3 = line[1]
    z3 = line[2]
    x4 = (line+1)[0]
    y4 = (line+1)[1]
    z4 = (line+1)[2]
    neus = geo.neu(x3,y3,z3,x4,y4,z4)

tablica1 = np.column_stack((XYZ, PLH, XY_GK92, WSP_92, XY_GK2000, WSP_2000))

# zapis: https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.savetxt.html
np.savetxt("wsp_out.txt", tablica1, delimiter=',', header = 'konversja współrzednych geodezyjnych')


