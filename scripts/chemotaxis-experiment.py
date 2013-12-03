import os
import itertools

simname = 'chemotaxis-experiment'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

method_solvers = { #'ForwardEuler':  ['StepDoublingSolver'],
#				   'BackwardEuler': ['StepDoublingSolver'],
#				   'RK4':           ['StepDoublingSolver'],
#				   'RODAS':         ['EmbeddedSolver'],
				   'Radau5':        ['EmbeddedSolver'],
				   'ARK3':          ['EmbeddedSolver'],
				   'ARK4':          ['EmbeddedSolver'],
				   'ARK5':          ['EmbeddedSolver'],
#				   'RKC1':          ['EmbeddedSolver'],
				   'RKC2':          ['EmbeddedSolver'] }
#				   'PRKC':          ['EmbeddedSolver'],
#				   'IRKC':          ['EmbeddedSolver'],
#				   'RKF45':         ['EmbeddedSolver'],
#				   'Merson43':      ['EmbeddedSolver'],
#				   'Verner65':      ['EmbeddedSolver'],
#				   'DOPR54':        ['EmbeddedSolver'] }

methods = [ m for m in method_solvers ]
ark_methods = [ 'ARK3', 'ARK4', 'ARK5' ]
ivps = ['Angiogenesis1D']

tolerances = [(t,t) for t in [ 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12 ]]
unknowns = (500,1000)
deltas = (1.,1e-3)
 
def GenerateRunList():
	tolArgs = [ {"dt": 0.01, "atol": t[0], "rtol": t[1], "newton tol": t[0]*1e-2} for t in tolerances ]

	runlist = []
	for N, d, m, ivp, tol in itertools.product(unknowns, deltas, methods, ivps, tolArgs):
		for ms in method_solvers[m]:
			argdict = {"ivp": ivp,
					   "method": m,
					   "solver": ms,
					   "N": N, "d": d,
					   "sparse": 1,
					   "jacobian": "Analytic",
					   "max steps": 100000,
					   "min write time": 0.01 }
			argdict.update(tol)

			if m in ark_methods:
				# Skip high tolerances for physical splitting
				if argdict["atol"] >= 1e-8:
					argdict.update({"jacobian splitting": 0})
					runlist.append(argdict.copy())

				argdict.update({"jacobian splitting": 1})
				runlist.append(argdict.copy())
			else: 
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
	return (0,0,0)

def GenerateAnalysisPasses():
	passes = []
	for N, d, ivp in itertools.product(unknowns, deltas, ivps):
		for offset in (0,1):
			solution1d = { 'mode': 'Solution1D',
						   'match': {'ivp': ivp, 'N': N, 'd': d, 'method': 'Radau5', 'atol': tolerances[-1][0] },
						   'title': ivp + ' Solution (N=%d)' % N,
						   'filename': ivp.lower() + '-%d-%.2e-%s' % (N,d,'c' if offset else 'p'),
						   'xlabel': 'N',
						   'ylabel': 'Value',
						   'plottxt': 1,
						   'solution times': [ 0, 0.1, 0.3, 0.5, 0.7 ],
						   'solution stride': 2,
						   'solution offset': offset }
			passes.append(solution1d)

		for matchUpdate in [ ('all',{}), ('ark',{'method': ['ARK3','ARK4','ARK5']}) ]:
			accuracy = { 'mode': 'Accuracy',
						 'plottxt': 1,
						 'legend': AccuracyLegendName,
						 'reference run': {'ivp':    ivp,
										   'N':      N,
										   'd':      d,
										   'method': 'Radau5',
										   'atol':   tolerances[-1][0] },
						 'match': {'ivp': ivp, 'N': N, 'd': d },
						 'group': ['method','solver','jacobian splitting'] }
			accuracy['match'].update(matchUpdate[1])

			if matchUpdate[0] == 'ark':
				accuracy['symbol'] = AccuracySymbol
				accuracy['color'] = AccuracyColor
				accuracy['xsize'] = 7.
				accuracy['ysize'] = 3.

			metrics = {"time": "CPU Time (ms)", "steps": "Steps"}		
			for m in metrics:
				title = ivp + ' -- ' + metrics[m] + ' vs. Accuracy' + \
					(" (N=%d, d=%.1e)" % (N,d))
				filename = ivp.lower() + '-%s-%d-%.1e-%s' % (m, N, d, matchUpdate[0])
				accuracy.update({ 'title': title, 'filename': filename, 'comparison': m, "xlabel": "Accuracy", "ylabel": metrics[m]})
				passes.append(accuracy.copy())
	
	return passes

