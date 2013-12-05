#!/usr/bin/env python
import sys
import os
import shutil
import subprocess
import numpy as np
import numpy.linalg as la

ivp = 'Angiogenesis1D'
runpath = os.path.join(os.getenv('HOME'), 'tmp', 'convergence-test')
jactype = 'Analytic'
N = '200'

def GetSolution(args):
	if os.path.exists(runpath):
		shutil.rmtree(runpath)
	os.makedirs(runpath)
	
	command = ['../pythODE++', '-phase', 'runner', '-path', runpath] + args
	p = subprocess.Popen(command)
	p.wait()
	
	p = subprocess.Popen(['../pythODE++', '-phase', 'dumpsolution', '-path', runpath], stdout=subprocess.PIPE)
	out,_ = p.communicate()

	if os.path.exists(runpath):
		shutil.rmtree(runpath)

	return np.array([float(x) for x in out.split("\n")[1].strip().split(" ") ])

exactsol = GetSolution(['-min write time', '1.',
						'-solver', 'EmbeddedSolver',
					    '-method', 'Radau5',
					    #'-method', 'ARK5',
					    '-jacobian splitting', '1',
					    #'-print time', '1',
					    '-atol', '1e-12',
					    '-rtol', '1e-12',
					    '-newton tol', '1e-14',
						'-sparse', '1',
						'-N', N,
						'-jacobian', jactype,
						'-ivp', ivp])

tols = { '0': ('1e-5','1e-6','1e-7','1e-8'),
		 '1': ('1e-5','1e-6','1e-7','1e-8','1e-9','1e-10','1e-11','1e-12') }
for js in ('0','1'):
	print 'vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv'
	for tol in tols[js]:
		solution = GetSolution(['-min write time', '20',
								'-solver', 'EmbeddedSolver',
							    '-method', 'ARK5',
							    '-atol', tol,
							    '-rtol', tol,
								'-newton tol', str(float(tol)*1e-2),
								'-sparse', '1',
								'-print stats', '1',
								'-N', N,
								'-jacobian', jactype,
								'-jacobian scaling', '1.1',
								'-jacobian splitting', js,
								'-ivp', ivp ])
		error = la.norm(solution - exactsol)
		print 'JS = '+js+', Tol = '+tol+', Error = ', error
	print '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'


