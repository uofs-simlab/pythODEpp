import os
import itertools

simname = 'advection-experiment'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

ivpname = 'AdvectionDiffusion1D'
method_solvers = { #'ForwardEuler':  ['StepDoublingSolver'],
				   #'BackwardEuler': ['StepDoublingSolver'],
				   #'RK4':           ['StepDoublingSolver'],
				   #'Celledoni1':    ['StepDoublingSolver'],
				   #'Celledoni2':    ['StepDoublingSolver'],
				   #'AddExpERK2':    ['StepDoublingSolver'],
				   'RODAS':         ['EmbeddedSolver'],
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

methods = [ m for m in method_solvers ]

ark_methods = [ 'ARK3', 'ARK4', 'ARK5' ]

tolerances = [(t,t) for t in [1e-4, 1e-5, 1e-6, 1e-7, 1e-8]]
unknowns     = (500,)
advList      = (1e1,1e2)
diffList     = (1e1,1e2)
fdOrderList  = (2,) # Order 4 not yet implemented

def GenerateRunList():
	tolArgs = [ {"dt": 0.01, "atol": t[0], "rtol": t[1]} for t in tolerances ]

	runlist = []
	for N, adv, diff, fdo, m, tol in itertools.product(unknowns, advList, diffList, fdOrderList, methods, tolArgs):
		for ms in method_solvers[m]:
			argdict = {"ivp": ivpname,
					   "method": m,
					   "solver": ms,
					   "N": N,
					   "sparse": 1,
					   "fdorder": fdo,
					   "adv": adv,
					   "diff": diff,
					   "max steps": 10000000,
					   "jacobian": "Analytic",
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
	if runinfo['solver'] == 'ConstantSolver':
		return runinfo['method name'] + ", h=" + str(runinfo['dt'])
	else:
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
	for N, adv, diff, fdo in itertools.product(unknowns, advList, diffList, fdOrderList):
		for matchUpdate in [ ('all',{}), ('ark',{'method': ['ARK3','ARK4','ARK5']}),
							('arkrkc',{'method': ['ARK3','ARK4','ARK5','RKC1','RKC2']}) ]:	
			if matchUpdate[0] in ('ark','all'):
				continue
			accuracy = { 'mode': 'Accuracy',
						 'plottxt': 1,
						 'legend': AccuracyLegendName,
						 'reference run': {'ivp':     ivpname,
										   'N':       N,
										   'fdorder': fdo,
										   'adv':     adv,
										   'diff':    diff,
										   'method':  'RODAS',
										   'atol':    tolerances[-1][0] },
						 'match': {'ivp': ivpname, 'N': N, 'adv': adv, 'diff': diff },
						 'group': ['method','solver','jacobian splitting'] }
			accuracy['match'].update(matchUpdate[1])

			if matchUpdate[0] != 'all':
				accuracy['symbol'] = AccuracySymbol
				accuracy['color'] = AccuracyColor
				accuracy['xsize'] = 7.
				accuracy['ysize'] = 3.

			metrics = {"time": "CPU Time (ms)", "steps": "Steps"}		
			for m in metrics:
				title = "%d Unknowns, Adv=%g, Diff=%g, FD Order=%d" % (N, adv, diff, fdo)
				filename = 'advection-diffusion1d-%s-n%d-fd%d-a%.2f-d%.2f-%s' % (m, N, fdo, adv, diff, matchUpdate[0])
				accuracy.update({ 'title': title, 'filename': filename, 'comparison': m, "xlabel": "Accuracy", "ylabel": metrics[m]})
				passes.append(accuracy.copy())
	
	return passes

