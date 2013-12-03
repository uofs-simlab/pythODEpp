#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys

def cheb(j,z):
	if j == 0:
		return 1.
	if j == 1:
		return 1.*z
	return 2*z*cheb(j-1,z) - cheb(j-2,z)

def chebp(j,z):
	if j == 0:
		return 0.
	if j == 1:
		return 1.
	return 2*cheb(j-1,z) + 2*z*chebp(j-1,z) - chebp(j-2,z)

def chebpp(j,z):
	if j == 0:
		return 0.
	if j == 1:
		return 0.
	return 4*chebp(j-1,z) + 2*z*chebpp(j-1,z) - chebpp(j-2,z)

epsilon = 2./13

def rkc1_stab(z,s,w0,w1):
	br = 1. / cheb(s,w0)
	w1 = cheb(s,w0)/chebp(s,w0)
	return br*cheb(s, w0 + w1*z)

def rkc2_stab(z,s,w0,w1):
	br = chebpp(s,w0) / chebp(s,w0)**2
	ar = 1 - br*cheb(s,w0)
	return ar + br*cheb(s, w0 + w1*z)

def prkc_stab(z,s,w0,w1):
	la = z.real
	mu = 1j*z.imag
	if s == 2:
		cmm1 = w1*chebpp(2,w0)/chebp(2,w0)
	else:
		cmm1 = w1*chebpp(s-1,w0)/chebp(s-1,w0)
	
	v = 1.
	a3 = 0

	a0 = 1./2
	a1 = -1./2 + v*(3-4*v)
	a2 = 2*v*(2*v-1) - a3
	a4 = (1-3*v)/(6*v)
	a5 = (1+3*v*(1-2*v)+4*cmm1*v*(3*v-2)) / (6*cmm1*v*(2*v-1))
	a6 = (3*v*(2*v-1) - 1) / (6*cmm1*v*(2*v-1))
	a7 = 1./(6*v*(2*v-1))

	stab = rkc2_stab(la,s,w0,w1)*(1 + (a0+a7)*mu + a0*a7*mu**2)
	stab += rkc2_stab(la,s-1,w0,w1)*(a6*mu + (a0*a6+a3*a7)*mu**2 + a0*a3*a7*mu**3)
	stab += (a4+a5)*mu + (a0*a5 + (a1+a2)*a7)*mu**2 + a0*a2*a7*mu**3
	return stab

def irkc_stab(z,s,w0,w1):
	br = chebpp(s,w0) / chebp(s,w0)**2
	ar = 1 - br*cheb(s,w0)
	return ar + br*cheb(s, w0 + w1*z/(1-w1/w0*(z.imag*1j)))

for method in sys.argv[1:]:
	for r in range(2,8):
		fig = plt.figure(figsize=(5, 5), dpi=80)
	
		xmin = -6*r + 8
		xmax = 2
		nx = 200
		ny = 200
		dt = 1
	
		if method == "prkc":
			ymin = -4
			ymax = 4
			stabf = prkc_stab
			title = str(r) + '-stage PRKC'
			filename = 'prkc-' + str(r) + 'stage-'
		elif method == "irkc":
			ymin = -20
			ymax = 20
			stabf = irkc_stab
			title = str(r) + '-stage IRKC'
			filename = 'irkc-' + str(r) + 'stage-'
		elif method == "rkc1":
			xmin = -18*r + 26
			ymin = -8
			ymax = 8
			stabf = rkc1_stab
			title = str(r) + '-stage RKC1'
			filename = 'rkc1-' + str(r) + 'stage-'
		elif method == "rkc2":
			ymin = -6
			ymax = 6
			stabf = rkc2_stab
			title = str(r) + '-stage RKC2'
			filename = 'rkc2-' + str(r) + 'stage-'
		else:
			print "invalid rkc type."
			sys.exit(1)
		
		x, y = np.meshgrid(np.linspace(xmin,xmax,nx),np.linspace(ymin,ymax,ny))
		D = x + 1j*y
		w0 = 1 + epsilon/r**2	
		w1 = chebp(r,w0)/chebpp(r,w0)
		
		f  = np.zeros(D.shape)
		for i in range(D.shape[0]):
			for j in range(D.shape[1]):
				f[i,j] = abs(stabf(D[i,j],r,w0,w1))
	
		plt.contourf(x, y, f, [0, 1], colors=['#999999'])
		plt.contour(x, y, f, [1], colors='k', linewidths=1)
	
		plt.plot([xmin, xmax],[0, 0],'k')
		plt.plot([0, 0],[ymin, ymax],'k')
	
		plt.axis([xmin, xmax, ymin, ymax])
		plt.xlabel('$(\Delta t) \lambda$')
		plt.ylabel('$(\Delta t) \mu$')
		#plt.title('Region of Absolute Stability for ' + title)
		plt.text(0.05, 0.05, 'm='+str(r), style='italic', bbox={'facecolor': 'white', 'pad':10}, transform=plt.axes().transAxes)
	
		name = filename+"stab.pdf"
		print "writing", name
		plt.savefig(name,bbox_inches='tight')
	
