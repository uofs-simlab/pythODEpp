#!/usr/bin/env python
import sys
import os
import shutil
import subprocess
import numpy as np
import numpy.linalg as la

verbose = False

def GetSolution(args):
	if os.path.exists(runpath):
		shutil.rmtree(runpath)
	os.makedirs(runpath)

	command = ['../pythode++', '-phase', 'runner', '-solver', 'ConstantSolver', '-path', runpath] + args
	if verbose:
		print ' '.join(command)
	p = subprocess.Popen(command)
	p.wait()
	
	p = subprocess.Popen(['../pythode++', '-phase', 'dumpsolution', '-path', runpath], stdout=subprocess.PIPE)
	out,_ = p.communicate()

	if os.path.exists(runpath):
		shutil.rmtree(runpath)

	return np.array([float(x) for x in out.split("\n")[1].strip().split(" ") ])

methods = ['RK4','ForwardEuler','DIRKCF1','DIRKCF2','DIRKCF3','ERKCF2']
ivp = 'AdvectionDiffusion1D'
dtstart = 1e-2

solparams = ['-N','10','-fdorder','2','-sparse','1']
runpath = os.path.join(os.getenv('HOME'), 'tmp', 'convergence-test')
solexact = GetSolution(['-min write time', '20', '-method', 'RK4', '-ivp', ivp, '-dt', '1e-5'] + solparams)

for m in methods:
	errorList = []
	for dt in [dtstart*2**-i for i in range(10)]:
		solution = GetSolution(['-min write time', '20',
							    '-method', m,
							    '-ivp', ivp,
							    '-dt', str(dt)] + solparams)
		error = la.norm(solution - solexact)
		errorList.append(error)

	convergence = np.array(errorList[:-1]) / np.array(errorList[1:])
	print m, ' '.join([ '%7.4f' % c for c in convergence ])

