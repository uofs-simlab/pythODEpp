import os
import itertools

simname = 'heat-transfer'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

method_solvers = { 'DOPR54':        ['EmbeddedSolver'],
				   'RKC1':        ['EmbeddedSolver'],
				   'RKC2':        ['EmbeddedSolver'],
				   'ARK3':        ['EmbeddedSolver'],
				   'ARK4':        ['EmbeddedSolver'],
				   'ARK5':        ['EmbeddedSolver'] }

methods = [ m for m in method_solvers ]
ark_methods = [ 'ARK3', 'ARK4', 'ARK5' ]
ivps = ['HeatTransfer']
reference_method = 'DOPR54'

tolerances = [(10**-t,10**-t) for t in range(4,9)]
diffusivities = (1e-7,)

unknowns = ((70,10),(140,20),(210,30),(280,40))
#unknowns = ((70,10),(350,50),(700,100))
 
def GenerateRunList():
	tolArgs = [ {"dt": 0.01, "atol": t[0], "rtol": t[1], "newton tol": t[0]*1e-2 } for t in tolerances ]

	runlist = []
	for N, d, m, ivp in itertools.product(unknowns, diffusivities, methods, ivps):
		for ms in method_solvers[m]:
			for tol in tolArgs:
				argdict = {"ivp": ivp,
						   "method": m,
						   "solver": ms,
						   "d":  d,
						   "NX": N[0],
						   "NY": N[1],
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
	for N, d, ivp in itertools.product(unknowns, diffusivities, ivps):
		referenceMatch = {'ivp':      ivp,
						  'd':        d,
						  'NX':       N[0],
						  'NY':       N[1],
						  'method':   reference_method,
						  'atol':     tolerances[-1][0] }
		for t in [5*i for i in range(11)]:
			passes.append({ 'mode': 'Solution2D',
							'match': referenceMatch,
							'filename': ivp.lower() + '-%dx%d-%.2f' % (N[0],N[1],t),
							'title': 'NX=%d, NY=%d, t=%.3f' % (N[0],N[1],t),
							'xlabel': 'X',
							'ylabel': 'Y',
							'xdim': N[0],
							'ydim': N[1],
							'xsize': 10,
							'ysize': 3,
							'cbmin': 20 + 273.15,
							'cbmax': 40 + 273.15,
							'solution time': t })

		for matchUpdate in [ ('all',{}), ('ark',{'method': ['ARK3','ARK4','ARK5']}), ('arkrkc',{'method': ['ARK3','ARK4','ARK5','RKC1','RKC2']}) ]:
			if matchUpdate[0] == 'ark':
				continue

			accuracy = { 'mode': 'Accuracy',
						 'plottxt': 1,
						 'legend': AccuracyLegendName,
						 'reference run': referenceMatch,
						 'match': {'ivp': ivp, 'NX': N[0], 'NY': N[1], 'd': d },
						 'group': ['method','solver','jacobian splitting'] }
			accuracy['match'].update(matchUpdate[1])

			if matchUpdate[0] in ('ark','arkrkc'):
				accuracy['symbol'] = AccuracySymbol
				accuracy['color'] = AccuracyColor
				accuracy['xsize'] = 7.
				accuracy['ysize'] = 3.

			metrics = {"time": "CPU Time (ms)", "steps": "Steps"}		
			for m in metrics:
				title ="%dx%d Unknowns" % (N[0],N[1])
				filename = ivp.lower() + '-%s-%dx%d-%s' % (m, N[0], N[1], matchUpdate[0])
				accuracy.update({ 'title': title, 'filename': filename, 'comparison': m, "xlabel": "Accuracy", "ylabel": metrics[m]})
				passes.append(accuracy.copy())
	
	return passes

