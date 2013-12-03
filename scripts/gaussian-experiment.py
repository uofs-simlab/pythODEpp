import os
import itertools

simname = 'gaussian-experiment'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

method_solvers = { #'ForwardEuler':  ['StepDoublingSolver'],
#				   'BackwardEuler': ['StepDoublingSolver'],
#				   'RK4':           ['StepDoublingSolver'],
				   'RODAS':         ['EmbeddedSolver'],
				   'Radau5':        ['EmbeddedSolver'],
#				   'ARK3':          ['EmbeddedSolver'],
#				   'ARK4':          ['EmbeddedSolver'],
#				   'ARK5':          ['EmbeddedSolver'],
#				   'RKC2':          ['EmbeddedSolver'],
#				   'PRKC':          ['EmbeddedSolver'],
#				   'IRKC':          ['EmbeddedSolver'],
#				   'RKF45':         ['EmbeddedSolver'],
#				   'Merson43':      ['EmbeddedSolver'],
#				   'Verner65':      ['EmbeddedSolver'],
				   'DOPR54':        ['EmbeddedSolver'] }

methods = [ m for m in method_solvers ]
ark_methods = [ 'ARK3', 'ARK4', 'ARK5' ]
ivps = ['Gaussian']

tolerances = [(t,t) for t in [1e-4, 1e-5, 1e-6, 1e-7, 1e-8 ]]
unknowns = (100,)
 
def GenerateRunList():
	tolArgs = [ {"dt": 0.01, "atol": t[0], "rtol": t[1]} for t in tolerances ]

	runlist = []
	for N, m, ivp in itertools.product(unknowns, methods, ivps):
		for ms in method_solvers[m]:
			for tol in tolArgs:
				argdict = {"ivp": ivp,
						   "method": m,
						   "solver": ms,
						   "N": N,
						   "sparse": 1,
						   "jacobian": "Forward",
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

def GenerateAnalysisPasses():
	passes = []
	for N, ivp in itertools.product(unknowns, ivps):
		passes.append({ 'mode': 'Solution1D',
						'match': {'ivp': ivp, 'N': N, 'method': 'RODAS', 'atol':   tolerances[-1][0] },
						'title': "Gaussian Plot",
						'filename': "gaussian",
						'xlabel': 'N',
						'ylabel': 'Value',
						'xsize': 6,
						'ysize': 4,
						'solution times': [ t*1e0 for t in range(5)] })

		for matchUpdate in [ ('all',{}), ('ark',{'method': ['ARK3','ARK4','ARK5']}) ]:
			accuracy = { 'mode': 'Accuracy',
						 'legend': AccuracyLegendName,
						 'reference run': {'ivp':    ivp,
										   'N':      N,
										   'method': 'RODAS',
										   'atol':   tolerances[-1][0] },
						 'match': {'ivp': ivp, 'N': N },
						 'group': ['method','solver','jacobian splitting'] }
			accuracy['match'].update(matchUpdate[1])

			if matchUpdate[0] == 'ark':
				accuracy['symbol'] = AccuracySymbol
				accuracy['xsize'] = 7.
				accuracy['ysize'] = 3.

			metrics = {"time": "CPU Time (ms)", "steps": "Steps"}		
			for m in metrics:
				title ="%d Unknowns" % N
				filename = ivp.lower() + '-%s-%d-%s' % (m, N, matchUpdate[0])
				accuracy.update({ 'title': title, 'filename': filename, 'comparison': m, "xlabel": "Accuracy", "ylabel": metrics[m]})
				passes.append(accuracy.copy())
	
	return passes

