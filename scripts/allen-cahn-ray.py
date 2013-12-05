import os
import itertools

simname = 'allen-cahn-ray'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

method_solvers = { 'ForwardEuler': ['StepDoublingSolver'],
				   'BackwardEuler': ['StepDoublingSolver'],
				   'ARK1': ['StepDoublingSolver'],
				   'ARK3': ['StepDoublingSolver'],
				   'Radau5': ['EmbeddedSolver'] }
methods = [ m for m in method_solvers ]
ark_methods = [ 'ARK1', 'ARK3', 'ARK4', 'ARK5' ]
ivps = ['AllenCahn']
reference_method = 'Radau5'

tolerances = [(10**-t,10**-t) for t in range(4,10)]
unknowns = (50,200)
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
						   "jacobian": "Analytic",
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

def GenerateAnalysisPasses():
	passes = []
	for N, l, g, ivp in itertools.product(unknowns, lambdas, gammas, ivps):
		referenceMatch = {'ivp':      ivp,
						  'N':        N,
						  'lambda':   l,
						  'gamma':    g,
						  'method':   reference_method,
						  'atol':     tolerances[-1][0] }

		accuracy = { 'mode': 'Accuracy',
					 'plottxt': 1,
					 'legend': AccuracyLegendName,
					 'reference run': referenceMatch,
					 'match': {'ivp': ivp, 'N': N, 'lambda': l, 'gamma': g},
					 'group': ['method','solver','jacobian splitting'] }

		metrics = {"time": "CPU Time (ms)", "steps": "Steps"}		
		for m in metrics:
			title ="%dx%d Unknowns, lambda=%.1f, gamma=%.1f" % (N,N,l,g)
			filename = ivp.lower() + '-%s-%d-%.1f-%.1f' % (m, N, l, g)
			accuracy.update({ 'title': title, 'filename': filename, 'comparison': m, "xlabel": "Accuracy", "ylabel": metrics[m]})
			passes.append(accuracy.copy())
	
	return passes

