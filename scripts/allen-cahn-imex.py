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
	
	command = ['../pythode++', '-phase', 'runner', '-sparse', '1',
			   '-jacobian', 'Analytic', '-min write time', '0.02',
			   '-path', runpath, '-ivp', 'AllenCahn', '-N', '150' ] + args
	if verbose:
		print ' '.join(command)
	p = subprocess.Popen(command)
	p.wait()
	
	p = subprocess.Popen(['../pythode++', '-phase','dumpsolution','-path', runpath], stdout=subprocess.PIPE)
	out,_ = p.communicate()

	if os.path.exists(runpath):
		shutil.rmtree(runpath)

	return np.array([float(x) for x in out.split("\n")[1].strip().split(" ") ])

methods = [('ForwardEuler',  0.0002, '0'),
		   ('BackwardEuler', 0.0015, '0'),
		   ('ARK1',          0.0002, '0'),
		   ('ARK1',          0.0002, '1')]

runpath = os.path.join(os.getenv('HOME'), 'tmp', 'allen-cahn-test')
solexact = GetSolution(['-method', 'Radau5',
						'-solver', 'EmbeddedSolver',
						'-atol', '1e-10',
						'-rtol', '1e-10',
						'-newton tol', '1e-12'])

for m in sorted(methods):
	print m[0] if m[2] == '0' else m[0]+' (JS)'
	print 'dt=%7.5f' % m[1]
	solution = GetSolution(['-solver', 'ConstantSolver',
						    '-method', m[0],
							'-print stats', '1',
							'-dt', str(m[1]),
						    '-jacobian splitting', m[2]])
	error = la.norm(solution - solexact) #/la.norm(solexact)
	print 'Error: %.10f' % error
	print '---------------'

