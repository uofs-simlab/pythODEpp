#!/usr/bin/env python
import sys
import os
import shutil
import subprocess

basepath = os.path.join(os.getenv('HOME'), 'tmp', 'eigenvalue-plots')
runpath = os.path.join(basepath,'runs')
plotpath = os.path.join(basepath,'plots')

if os.path.exists(basepath):
	shutil.rmtree(basepath)
os.makedirs(plotpath)

def GetEigenvalues(plotname, ivp, title, dt, params=''):
	os.makedirs(runpath)
	print 'Running',ivp,'for',plotname
	subprocess.call('../pythODE++ -phase runner -method Radau5 -ivp ' + ivp + ' -solver EmbeddedSolver ' + params + ' -path ' + runpath, shell=True)
	plotfile = os.path.join(plotpath,plotname)
	print '  setting title.'
	subprocess.call('echo set title \\\''+title+'\\\' > ' + plotfile, shell=True)
	print '  computing eigenvalues.'
	subprocess.call('../pythODE++ -phase eigenvalues -ivp ' + ivp + ' -path ' + runpath + ' -dt ' + dt + ' ' + params + ' | python >> ' + plotfile, shell=True)
	print '  writing pdf.'
	subprocess.call('gnuplot < ' + plotfile + ' > ' + plotfile + '.pdf', shell=True)
	print '  complete.'
	shutil.rmtree(runpath)

GetEigenvalues('crete1d-sink.txt', 'ConcreteRewetting', 'Sink Boundary Condition', '13.999', '-sink\ bc 1')
GetEigenvalues('crete1d-insulated.txt', 'ConcreteRewetting', 'Insulated Boundary Condition', '13.999', '-sink\ bc 0')
GetEigenvalues('adv1d.txt', 'AdvectionDiffusion1D', 'Linear Advection-Diffusion', '1.')
GetEigenvalues('heat1d.txt', 'HeatTransfer', 'Heat Transfer', '50.', '-NX 70 -NY 10')
GetEigenvalues('cusp1d.txt', 'CUSP', 'CUSP Model', '0.5')
GetEigenvalues('bruss1d.txt', 'Brusselator1D', '1D Brusselator', '3.', '-N 100')
GetEigenvalues('bruss2d.txt', 'Brusselator2D', '2D Brusselator', '3.')
GetEigenvalues('angio1d-0.txt', 'Angiogenesis1D', 'd=1', '0.34', '-d 1 -N 200 -print\ time 1')
GetEigenvalues('angio1d-1.txt', 'Angiogenesis1D', 'd=0.001', '0.1', '-d 0.001 -atol 1e-11 -rtol 1e-11 -newton\ tol 1e-13 -print\ time 1 -sparse 1 -jacobian Forward -N 200')

for U in (0,0.75,0.99):
	for r in ('FKPP','Ignition','Fisher'):
		GetEigenvalues('combustion-U'+str(U)+'-'+r+'.txt', 'CombustionARD', r, '15.', '-reaction '+r+' -U0 '+str(U))

