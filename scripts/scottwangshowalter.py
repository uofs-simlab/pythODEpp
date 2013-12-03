import os
import itertools

simname = 'ScottWangShowalter'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

method_solvers = { #'ForwardEuler':  ['StepDoublingSolver'],
#				   'BackwardEuler': ['StepDoublingSolver'],
#				   'RK4':           ['StepDoublingSolver'],
#				   'RODAS':         ['EmbeddedSolver'],
#				   'Radau5':        ['EmbeddedSolver'],
				   'ARK3':          ['EmbeddedSolver'],
				   'ARK4':          ['EmbeddedSolver'],
				   'ARK5':          ['EmbeddedSolver'],
#				   'RKC2':          ['EmbeddedSolver'],
#				   'PRKC':          ['EmbeddedSolver'],
#				   'IRKC':          ['EmbeddedSolver'],
#				   'RKF45':         ['EmbeddedSolver'],
#				   'Merson43':      ['EmbeddedSolver'],
#				   'Verner65':      ['EmbeddedSolver'],
				   'DOPR54':        ['EmbeddedSolver'] }

methods = [ m for m in method_solvers ]
ark_methods = [ 'ARK3', 'ARK4', 'ARK5' ]
ivps = ['ScottWangShowalter']
refsol = 'DOPR54'

tolerances = [(10**-t,10**-t) for t in range(4,11)]
unknowns = (100,200,300)
problems = ('Problem1', 'Problem2', 'Problem3')
problemTimes = { 'Problem1': (0.011, 0.033, 0.059),
				 'Problem2': (0.008, 0.012, 0.018, 0.040),
				 'Problem3': (0.008, 0.020, 0.039, 0.057) }
 
def GenerateRunList():
	tolArgs = [ {"dt": 0.01, "atol": t[0], "rtol": t[1], "newton sol": t[0]*1e-2} for t in tolerances ]

	runlist = []
	for N, m, ivp, p in itertools.product(unknowns, methods, ivps, problems):
		for ms in method_solvers[m]:
			for tol in tolArgs:
				argdict = {"ivp": ivp,
						   "method": m,
						   "solver": ms,
						   "N": N,
						   "problem": p,
						   "sparse": 1,
						   "jacobian": "Analytic",
						   "max steps": 1000000,
						   "min write time": 1e-3 }
				argdict.update(tol)
				runlist.append(argdict.copy())
					
				if m in ark_methods:
					argdict.update({"jacobian splitting": 1})
					runlist.append(argdict.copy())
	
	timingList = []
	timingGroup = 0
	for r in runlist:
		for t in range(1):
			r.update({"timing group": timingGroup})
			timingList.append(r.copy())
		timingGroup += 1
		
	return timingList

def SolutionLegendName(runinfo):
	return runinfo['method name'] + ", tol=" + str(runinfo['rtol'])

def AccuracyLegendName(runinfo):
	common = runinfo[0][1]['method name']
	if runinfo[0][1]['method'] in ark_methods:
		common += '(' + ('Jacobian Splitting' if 'jacobian splitting' in runinfo[0][1] and int(runinfo[0][1]['jacobian splitting']) else 'Physical Splitting') + ')'
	return common

def AccuracySymbol(runinfo):
	if 'jacobian splitting' in runinfo[0][1] and runinfo[0][1]['jacobian splitting'] == '1':
		return 7
	if runinfo[0][1]['method'] in ('RKC1','RKC2'):
		return 6
	return 8

def AccuracyColor(runinfo):
	method = runinfo[0][1]['method']
	if method == 'ARK3':
		return (1.0,0.0,1.0)
	if method == 'ARK4':
		return (0.7,0.0,0.0)
	if method == 'ARK5':
		return (0.0,0.0,1.0)
	if method == 'RKC1':
		return (0.5,1.0,0.5)
	if method == 'RKC2':
		return (0.7,0.7,0.7)
	return (0,0,0)

def GenerateAnalysisPasses():
	passes = []
	for N, p, ivp in itertools.product(unknowns, problems, ivps):
		referenceMatch = {'ivp':      ivp,
						  'N':        N,
						  'problem':  p,
						  'method':   refsol,
						  'atol':     tolerances[-1][0] }
		for t in problemTimes[p]:
			for q in (0,1):
				passes.append({ 'mode': 'Solution2D',
								'match': referenceMatch,
								'filename': ivp.lower() + '-n=%d-%s-t=%.3f-%s' % (N,p,t,'temperature' if q else 'concentration'),
								'title': '%s, N=%d, %s, t=%.3f' % ('Temperature' if q else 'Concentration',N,p,t),
								'xlabel': 'X',
								'ylabel': 'Y',
								'xsize': 6,
								'ysize': 4,
								'xdim': N,
								'ydim': N,
								'cbmin': 0,
								'cbmax': 60 if q else 3,
								'solution offset': q,
								'solution stride': 2,
								'solution time': t })

		for matchUpdate in [ ('all',{}), ('ark',{'method': ['ARK3','ARK4','ARK5']}) ]:
			accuracy = { 'mode': 'Accuracy',
						 'legend': AccuracyLegendName,
						 'reference run': referenceMatch,
						 'match': {'ivp': ivp, 'N': N, 'problem': p },
						 'group': ['method','solver','jacobian splitting'] }
			accuracy['match'].update(matchUpdate[1])

			if matchUpdate[0] == 'ark':
				accuracy['symbol'] = AccuracySymbol
				accuracy['color'] = AccuracyColor
				accuracy['xsize'] = 7.
				accuracy['ysize'] = 3.

			metrics = {"time": "CPU Time (ms)", "steps": "Steps"}		
			for m in metrics:
				title ="%d Unknowns, %s" % (N,p)
				filename = ivp.lower() + '-%s-%d-%s-%s' % (m, N, p, matchUpdate[0])
				accuracy.update({ 'title': title, 'filename': filename, 'comparison': m, "xlabel": "Accuracy", "ylabel": metrics[m]})
				passes.append(accuracy.copy())
	
	return passes

