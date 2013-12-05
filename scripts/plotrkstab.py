#!/usr/bin/env python

import os
import sys

import subprocess
import matplotlib.pyplot as plt

import numpy as np
import scipy as sp
from scipy import linalg as la

root_dir = os.path.dirname(os.path.realpath(sys.argv[0]))

if len(sys.argv) < 3:
	print 'method name and graph name required.'
	sys.exit(1)

method = sys.argv[1]
name = sys.argv[2]

fig = plt.figure(figsize=(6, 6), dpi=80)

xmin = -5
xmax = 5
ymin = -5
ymax = 5
nx = 100
ny = 100
dt = 1


def CallDump(method,coeff):
	return subprocess.check_output([os.path.join(root_dir,'..','pythODE++'),'-phase','runner','-method',method,coeff]).strip()

x = np.linspace(xmin, xmax, nx)
y = np.linspace(ymin, ymax, ny)
xv, yv = np.meshgrid(x, y)
z = np.zeros((nx,ny))

if method == "rodas":
	gamma = 0.25

	g = np.eye(6)*gamma
	g[1,0] = -0.3543
	g[2,0] = -0.133602505268175
	g[2,1] = -0.012897494731825
	g[3,0] =  1.526849173006459
	g[3,1] = -0.533656288750454
	g[3,2] = -1.279392884256
	g[4,0] =  6.981190951784981
	g[4,1] = -2.092930097006103
	g[4,2] = -5.870067663032724
	g[4,3] =  0.731806808253845
	g[5,0] = -2.080189494180926
	g[5,1] =  0.59576235567668
	g[5,2] =  1.701617798267255
	g[5,3] = -0.088514519835879
	g[5,4] = -0.378676139927128

	a = sp.zeros((6,6))
	a[1,0] = 0.386
	a[2,0] = 0.146074707525418
	a[2,1] = 0.063925292474582
	a[3,0] = -0.330811503667722
	a[3,1] = 0.711151025168282
	a[3,2] = 0.24966047849944
	a[4,0] = -4.552557186318003
	a[4,1] = 1.710181363241322
	a[4,2] = 4.014347332103150
	a[4,3] = -0.171971509026469
	a[5,0] = 2.428633765466978
	a[5,1] = -0.382748733764781
	a[5,2] = -1.855720330929574
	a[5,3] = 0.559835299227375
	a[5,4] = 0.25

	b = sp.zeros(6)
	b[0] = 0.348444271286054
	b[1] = 0.213013621911897
	b[2] = -0.154102532662319
	b[3] = 0.471320779391497
	b[4] = -0.128676139927129
	b[5] = 0.25
	
	A  = a + g
	At = A
	bt = b

elif method == "Radau5":
	A = sp.array([ [(88.-7.*np.sqrt(6.))/360.,     (296.-169.*np.sqrt(6.))/1800., (-2.+3.*np.sqrt(6.))/225.],
				   [(296.+169.*np.sqrt(6.))/1800., (88.+7.*np.sqrt(6.))/360.,     (-2.-3.*np.sqrt(6.))/225.],
				   [(16.-np.sqrt(6.))/36.,            (16.+np.sqrt(6.))/36.,             1./9] ])
	b = sp.array([(16.-np.sqrt(6.))/36., (16.+np.sqrt(6.))/36., 1./9])
	At = A
	bt = b
	
else:
	A = np.array(map(lambda x: map(float, x.strip().split(' ')), CallDump(method,'A').split('\n')))
	b = np.array(map(float, CallDump(method,'b').split(' ')))
	At = np.array(map(lambda x: map(float, x.strip().split(' ')), CallDump(method,'A2').split('\n')))
	bt = np.array(map(float, CallDump(method,'b2').split(' ')))

s = len(b)

for i in range(nx):
	for j in range(ny):
		lam = xv[i,j]
		mu = yv[i,j]*1j
		z[i,j] = abs(1+sum(dt*(lam*b + mu*bt).dot(la.inv(np.eye(s,s) - lam*dt*A - mu*dt*At))))

plt.contourf(xv, yv, z, [0, 1], colors=['#999999'])
plt.contour(xv, yv, z, [1], colors='k', linewidths=1)

plt.plot([xmin, xmax],[0, 0],'k')
plt.plot([0, 0],[ymin, ymax],'k')

plt.axis([xmin, xmax, ymin, ymax])
plt.xlabel('$(\Delta t) \lambda$')
plt.ylabel('$(\Delta t) \mu$')
plt.title(name)
#plt.title(name + ' Region of Absolute Stability')

#plt.show()
plt.savefig(name+"stab.pdf",bbox_inches='tight')

