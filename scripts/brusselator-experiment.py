import os
import itertools

simname = 'brusselator-experiment'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

method_solvers = { #'ForwardEuler':  ['StepDoublingSolver'],
#				   'BackwardEuler': ['StepDoublingSolver'],
#				   'RK4':           ['StepDoublingSolver'],
#				   'RODAS':         ['EmbeddedSolver'],
#				   'Radau5':        ['EmbeddedSolver'],
				   'ARK3':          ['EmbeddedSolver'],
				   'ARK4':          ['EmbeddedSolver'],
				   'ARK5':          ['EmbeddedSolver'],
				   'RKC1':          ['EmbeddedSolver'],
				   'RKC2':          ['EmbeddedSolver'],
#				   'PRKC':          ['EmbeddedSolver'],
#				   'IRKC':          ['EmbeddedSolver'],
#				   'RKF45':         ['EmbeddedSolver'],
#				   'Merson43':      ['EmbeddedSolver'],
#				   'Verner65':      ['EmbeddedSolver'],
				   'Radau5':        ['EmbeddedSolver'] }

methods = [ m for m in method_solvers ]
ark_methods = [ 'ARK3', 'ARK4', 'ARK5' ]

ivp_unknowns = {'Brusselator1D': (500, 1000),
				'Brusselator2D': (50, 75, 100) }
#ivps = [ i for i in ivp_unknowns ]
#ivps = ['Brusselator1D']
ivps = ['Brusselator2D']


tolerances = [(t,t) for t in [1e-4, 1e-5, 1e-6, 1e-7, 1e-8 ]]
alphas = (2e-2,2e-3,)
 
def GenerateRunList():
	tolArgs = [ {"dt": 0.01, "atol": t[0], "rtol": t[1], "newton tol": t[0]*1e-2} for t in tolerances ]

	runlist = []
	for m, ivp, alpha in itertools.product(methods, ivps, alphas):
		for ms, N in itertools.product(method_solvers[m], ivp_unknowns[ivp]):
			for tol in stepArgs if ms == 'ConstantSolver' else tolArgs:
				argdict = {"ivp": ivp,
						   "method": m,
						   "solver": ms,
						   "N": N,
						   "sparse": 1,
						   "alpha": alpha,
						   "jacobian": "Analytic",
						   #"step control": "Predictive",
						   "max steps": 10000000,
						   "min write time": 0.1 }
				argdict.update(tol)
				runlist.append(argdict.copy())
					
				if m in ark_methods:
					argdict.update({"jacobian splitting": 1})
					runlist.append(argdict.copy())
	
	timingList = []
	timingGroup = 0
	for r in runlist:
		for t in range(3):
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
	for ivp, alpha in itertools.product(ivps, alphas):
		for N in ivp_unknowns[ivp]:
			for matchUpdate in [ ('all',{}), ('ark',{'method': ['ARK3','ARK4','ARK5']}), ('arkrkc',{'method': ['ARK3','ARK4','ARK5','RKC1','RKC2']}) ]:
				if matchUpdate[0] in ('all','ark'):
					continue

				accuracy = { 'mode': 'Accuracy',
							 'legend': AccuracyLegendName,
							 'plottxt': 1,
							 'reference run': {'ivp':    ivp,
											   'N':      N,
											   'method': 'Radau5',
											   'alpha':  alpha,
											   'atol':   tolerances[-1][0] },
							 'match': {'ivp': ivp, 'N': N, 'alpha': alpha },
							 'group': ['method','solver','jacobian splitting'] }
				accuracy['match'].update(matchUpdate[1])
	
				if matchUpdate[0] in ('ark','arkrkc'):
					accuracy['symbol'] = AccuracySymbol
					accuracy['color'] = AccuracyColor
					accuracy['xsize'] = 7.
					accuracy['ysize'] = 3.

				metrics = {"time": "CPU Time (ms)", "steps": "Steps"}		
				for m in metrics:
					if ivp == 'Brusselator1D':
						title = "%d Unknowns, Alpha=%g" % (N,alpha)
					else:
						title = "%dx%d Unknowns, Alpha=%g" % (N,N,alpha)

					filename = ivp.lower() + '-%s-%d-%g-%s' % (m, N, alpha, matchUpdate[0])
					accuracy.update({ 'title': title, 'filename': filename, 'comparison': m, "xlabel": "Accuracy", "ylabel": metrics[m]})
					passes.append(accuracy.copy())
	
	return passes

