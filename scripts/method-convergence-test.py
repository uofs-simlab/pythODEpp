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
	
	p = subprocess.Popen(['../pythode++', '-phase','dumpsolution','-path', runpath], stdout=subprocess.PIPE)
	out,_ = p.communicate()

	if os.path.exists(runpath):
		shutil.rmtree(runpath)

	return np.array([float(x) for x in out.split("\n")[1].strip().split(" ") ])

erk = ['ForwardEuler', 'Heun2', 'Runge2', 'Heun3', 'Kutta3', 'Runge3', 'RK4', 'RK38', 'RKF45', 'Merson43', 'Zonneveld43', 'DOPR54', 'Verner65', 'FEHL78']
ark = ['ARK3','ARK4','ARK5']
rkc = ['RKC1','RKC2','PRKC']
irk = ['BackwardEuler','RODAS','Radau5']
exprk = ['DIRKCF1','ERKCF2']
methods = erk + ark + rkc + irk + exprk

ivp = 'NonstiffE5'
solexact = [1.411797390542629e+01,  2.400000000000002e+00]

runpath = os.path.join(os.getenv('HOME'), 'tmp', 'convergence-test')
dtstart = 2e0

for m in methods:
	errorList = []
	for dt in [dtstart*2**-i for i in range(10)]:
		solution = GetSolution(['-min write time', '20',
								'-flipexp', '1',
							    '-method', m,
								'-sparse', '1',
								'-jacobian','Forward',
							    '-jacobian splitting', '1',
							    '-ivp', ivp,
							    '-dt', str(dt)])
		error = la.norm(solution - np.array(solexact))
		errorList.append(error)

	convergence = np.array(errorList[:-1]) / np.array(errorList[1:])
	print m, ' '.join([ '%7.4f' % c for c in convergence ])

