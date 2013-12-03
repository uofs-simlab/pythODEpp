import os
import itertools

simname = 'concrete-experiment'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

ivpname = 'ConcreteRewetting'
method_solvers = { #'ForwardEuler':  ['StepDoublingSolver'],
#				   'BackwardEuler': ['StepDoublingSolver'],
#				   'RK4':           ['StepDoublingSolver'],
#				   'RODAS':         ['EmbeddedSolver'],
#				   'Radau5':        ['EmbeddedSolver'],
				   'ARK3':          ['EmbeddedSolver'],
				   'ARK4':          ['EmbeddedSolver'],
				   'ARK5':          ['EmbeddedSolver'],
				   'RKC1':          ['EmbeddedSolver'],
				   'RKC2':          ['EmbeddedSolver'] }
#				   'PRKC':          ['EmbeddedSolver'],
#				   'RKF45':         ['EmbeddedSolver'],
#				   'Merson43':      ['EmbeddedSolver'],
#				   'Verner65':      ['EmbeddedSolver'],
#				   'DOPR54':        ['EmbeddedSolver'] }

methods = [ m for m in method_solvers ]
ark_methods = [ 'ARK3', 'ARK4', 'ARK5' ]

stepsizes = [1e-1, 1e-2]
tolerances = [(t,t) for t in [1e-4, 1e-5, 1e-6, 1e-7 ]]
unknowns = (100,200)
isos = (0,)
 
def GenerateRunList():
	stepArgs = [ {"dt": s} for s in stepsizes ]
	tolArgs = [ {"dt": 0.01, "atol": t[0], "rtol": t[1], "newton tol": t[1]*1e-2} for t in tolerances ]

	runlist = []
	# concrete rewetting problem
	for N, iso, sbc in itertools.product(unknowns, isos, (0,1)):
		runlist.append({"ivp": 'ConcreteRewetting',
						#"step control": 'Predictive',
						"method": 'Radau5',
						"solver": 'EmbeddedSolver',
						"sparse": '1',
						"atol": 1e-9,
						"rtol": 1e-9,
						"newton tol": 1e-11,
						"N": N,
						"sink bc": sbc,
						"isopropanol": iso,
						"jacobian": "Forward",
						"min write time": 0.1})
		for m in methods:
			for ms in method_solvers[m]:
				for tol in stepArgs if ms == 'ConstantSolver' else tolArgs:
					argdict = {"ivp": 'ConcreteRewetting',
							  # "step control": 'Predictive',
							   "method": m,
							   "solver": ms,
							   "N": N,
							   "sparse": 0,
							   "sink bc": sbc,
							   "isopropanol": iso,
							   "jacobian": "Forward",
							   "min write time": 0.1 }
					argdict.update(tol)
					
					if m in ark_methods:
						for js in (0,1):
							argdict.update({"jacobian splitting": js})
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
	if runinfo['solver'] == 'ConstantSolver':
		return runinfo['method name'] + ", h=" + str(runinfo['dt'])
	else:
		return runinfo['method name'] + ", tol=" + str(runinfo['rtol'])

def AccuracyLegendName(runinfo):
	common = runinfo[0][1]['method name']
	if runinfo[0][1]['method'] in ark_methods:
		common += '(' + ('Jacobian Splitting' if int(runinfo[0][1]['jacobian splitting']) else 'Physical Splitting') + ')'
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
	for N, iso, sbc in itertools.product(unknowns, isos, (0,1)):
		for matchUpdate in [ ('all',{}), ('ark',{'method': ['ARK3','ARK4','ARK5']}),
										 ('arkrkc',{'method': ['ARK3','ARK4','ARK5','RKC1','RKC2']}) ]:
			#if matchUpdate[0] == 'ark':
			#	continue

			accuracy = { 'mode': 'Accuracy',
						 'legend': AccuracyLegendName,
						 'plottxt': 1,
						 'reference run': {'ivp': ivpname, "N": N, "sink bc": sbc, "isopropanol": iso, 'method': 'Radau5' },
						 'match': {'ivp': ivpname, "N": N, "sink bc": sbc, "isopropanol": iso},
						 'group': ['method','solver','jacobian splitting'] }
			accuracy['match'].update(matchUpdate[1])

			if matchUpdate[0] != 'all':
				accuracy['symbol'] = AccuracySymbol
				accuracy['color'] = AccuracyColor
				accuracy['xsize'] = 7.
				accuracy['ysize'] = 3.

			metrics = {"time": "CPU Time (ms)", "steps": "Steps"}
			for m in metrics:
				title = "%d unknowns, %s boundary, %s" % (N, "Sink" if sbc else "Insulated", "Isopropanol" if iso else "Water")
				filename = 'concrete-%s-%d-%s-%s-%s' % (m, N, "sink" if sbc else "insulated", "iso" if iso else "water", matchUpdate[0])
				accuracy.update({ 'title': title, 'filename': filename, 'comparison': m, 'xlabel': "Accuracy", 'ylabel': metrics[m]})
				passes.append(accuracy.copy())
	
	return passes

