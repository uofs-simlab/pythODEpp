import os
import itertools

simname = 'combustion'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

method_solvers = { #'ForwardEuler':  ['StepDoublingSolver'],
#				   'BackwardEuler': ['StepDoublingSolver'],
#				   'RK4':           ['StepDoublingSolver'],
#				   'RODAS':         ['EmbeddedSolver'],
				   'Radau5':        ['EmbeddedSolver'],
				   'ARK3':          ['EmbeddedSolver'],
				   'ARK4':          ['EmbeddedSolver'],
				   'ARK5':          ['EmbeddedSolver'],
				   'RKC1':          ['EmbeddedSolver'],
				   'RKC2':          ['EmbeddedSolver'] }
#				   'PRKC':          ['EmbeddedSolver'],
#				   'IRKC':          ['EmbeddedSolver'],
#				   'RKF45':         ['EmbeddedSolver'],
#				   'Merson43':      ['EmbeddedSolver'],
#				   'Verner65':      ['EmbeddedSolver'],
#				   'DOPR54':        ['EmbeddedSolver'] }

#methods = [ m for m in method_solvers ]
methods = [ 'Radau5' ]
ark_methods = [ 'ARK3', 'ARK4', 'ARK5' ]
ivps = ['CombustionARD']
uvalues = ( 0., 0.75, 0.99 )
reactions = ( 'FKPP', 'Ignition', 'Fisher' )

tolerances = [(10**-t,10**-t) for t in range(4,9)]
unknowns = (1600,)
 
def GenerateRunList():
	tolArgs = [ {"dt": 0.01, "atol": t[0], "rtol": t[1], "newton tol": t[0]*1e-2 } for t in tolerances ]

	runlist = []
	for N, m, ivp, U0, r in itertools.product(unknowns, methods, ivps, uvalues, reactions):
		for ms in method_solvers[m]:
			for tol in tolArgs:
				argdict = {"ivp": ivp,
						   "method": m,
						   "solver": ms,
						   "N": N,
						   "U0": U0,
						   "reaction": r,
						   "sparse": 1,
						   "jacobian": "Analytic",
						   "max steps": 1000000,
						   "min write time": 0.1 }
				argdict.update(tol)
				runlist.append(argdict.copy())
					
				if m in ark_methods:
					argdict.update({"jacobian splitting": 1})
					runlist.append(argdict.copy())
	
	timingList = []
	timingGroup = 0
	for r in runlist:
		for t in range(5):
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
	for N, ivp, U0, r in itertools.product(unknowns, ivps, uvalues, reactions):
		referenceMatch = {'ivp':      ivp,
						  'N':        N,
						  'U0':       U0,
						  'reaction': r,
						  'method':   'Radau5',
						  'atol':     tolerances[-1][0] }

		passes.append({ 'mode': 'Solution1D',
			'match': referenceMatch,
			'title': "Solution Plot, N=%d, U0=%.2f, %s" % (N,U0,r),
			'filename': "solutionPlot-n-%d--U0-%.2f-%s" % (N,U0,r),
			'xmin': 10.,
			'xmax': 50.,
			'xlabel': 'N',
			'ylabel': 'c(x,t)',
			'xsize': 6,
			'ysize': 4,
			'plottxt': 1,
			'solution times': [ 0., 15.] })

		for matchUpdate in [ ('all',{}), ('ark',{'method': ['ARK3','ARK4','ARK5','RKC1','RKC2']}) ]:
			accuracy = { 'mode': 'Accuracy',
						 'plottxt': 1,
						 'legend': AccuracyLegendName,
						 'reference run': referenceMatch,
						 'match': {'ivp': ivp, 'N': N, 'U0': U0, 'reaction': r },
						 'group': ['method','solver','jacobian splitting'] }
			accuracy['match'].update(matchUpdate[1])

			if matchUpdate[0] == 'ark':
				accuracy['symbol'] = AccuracySymbol
				accuracy['color'] = AccuracyColor
				accuracy['xsize'] = 7.
				accuracy['ysize'] = 3.
			else: # temporary, so we don't produce as many plots
				continue

			metrics = {"time": "CPU Time (ms)", "steps": "Steps"}		
			for m in metrics:
				title ="%s, N=%d, U_0=%.2f" % (r,N,U0)
				filename = ivp.lower() + '-%s-%d-%.2f-%s-%s' % (m, N, U0, r, matchUpdate[0])
				accuracy.update({ 'title': title, 'filename': filename, 'comparison': m, "xlabel": "Accuracy", "ylabel": metrics[m]})
				passes.append(accuracy.copy())
	
	return passes

