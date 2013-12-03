import os
import itertools

simname = 'allen-cahn'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

method_solvers = { 'Radau5':      ['EmbeddedSolver'],
				   'RKC1':        ['EmbeddedSolver'],
				   'RKC2':        ['EmbeddedSolver'],
				   'ARK3':        ['EmbeddedSolver'],
				   'ARK4':        ['EmbeddedSolver'],
				   'ARK5':        ['EmbeddedSolver'] }

methods = [ m for m in method_solvers ]
ark_methods = [ 'ARK3', 'ARK4', 'ARK5' ]
ivps = ['AllenCahn']
reference_method = 'Radau5'

tolerances = [(10**-t,10**-t) for t in range(4,9)]
unknowns = (100,200)
lambdas = (0.1,1.)
gammas = (0.1,1.)
 
def GenerateRunList():
	tolArgs = [ {"dt": 0.01, "atol": t[0], "rtol": t[1], "newton tol": t[0]*1e-2 } for t in tolerances ]

	runlist = []
	for N, l, g, m, ivp in itertools.product(unknowns, lambdas, gammas, methods, ivps):
		for ms in method_solvers[m]:
			for tol in tolArgs:
				argdict = {"ivp": ivp,
						   "method": m,
						   "solver": ms,
						   "N": N,
						   "lambda": l,
						   "gamma": g,
						   "sparse": 1,
						   "jacobian": "Autodiff",
						   "max steps": 1000000,
						   "min write time": 1. }
				argdict.update(tol)
				runlist.append(argdict.copy())
					
				if m in ark_methods:
					argdict.update({"jacobian splitting": 1})
					runlist.append(argdict.copy())
	
	timingList = []
	timingGroup = 0
	for r in runlist:
		for t in range(10):
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
	for N, l, g, ivp in itertools.product(unknowns, lambdas, gammas, ivps):
		referenceMatch = {'ivp':      ivp,
						  'N':        N,
						  'lambda':   l,
						  'gamma':    g,
						  'method':   reference_method,
						  'atol':     tolerances[-1][0] }
		for t in [10*i for i in range(11)]:
			passes.append({ 'mode': 'Solution2D',
							'match': referenceMatch,
							'filename': ivp.lower() + '-%d-%.2f' % (N,t),
							'title': '%dx%d, lambda=%.1f, gamma=%.1f, t=%.3f' % (N,N,l,g,t),
							'xlabel': 'X',
							'ylabel': 'Y',
							'xdim': N,
							'ydim': N,
							'xsize': 5,
							'ysize': 5,
							'cbmin': 0,
							'cbmax': 1,
							'solution time': t })

		for matchUpdate in [ ('all',{}), ('ark',{'method': ['ARK3','ARK4','ARK5']}), ('arkrkc',{'method': ['ARK3','ARK4','ARK5','RKC1','RKC2']}) ]:
			if matchUpdate[0] == 'ark':
				continue

			accuracy = { 'mode': 'Accuracy',
						 'plottxt': 1,
						 'legend': AccuracyLegendName,
						 'reference run': referenceMatch,
						 'match': {'ivp': ivp, 'N': N, 'lambda': l, 'gamma': g},
						 'group': ['method','solver','jacobian splitting'] }
			accuracy['match'].update(matchUpdate[1])

			if matchUpdate[0] in ('ark','arkrkc'):
				accuracy['symbol'] = AccuracySymbol
				accuracy['color'] = AccuracyColor
				accuracy['xsize'] = 7.
				accuracy['ysize'] = 3.

			metrics = {"time": "CPU Time (ms)", "steps": "Steps"}		
			for m in metrics:
				title ="%dx%d Unknowns, lambda=%.1f, gamma=%.1f" % (N,N,l,g)
				filename = ivp.lower() + '-%s-%d-%.1f-%.1f-%s' % (m, N, l, g, matchUpdate[0])
				accuracy.update({ 'title': title, 'filename': filename, 'comparison': m, "xlabel": "Accuracy", "ylabel": metrics[m]})
				passes.append(accuracy.copy())
	
	return passes

