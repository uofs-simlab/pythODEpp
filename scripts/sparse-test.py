import os
import itertools

simname = 'sparsity-test'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

method_solvers = { 'RODAS': ['EmbeddedSolver'],
				   'ARK3':  ['EmbeddedSolver'],
				   'ARK4':  ['EmbeddedSolver'],
				   'ARK5':  ['EmbeddedSolver'] }

methods = [ m for m in method_solvers ]
ark_methods = [ 'ARK3', 'ARK4', 'ARK5' ]
ivps = ['Brusselator1D']

stepsizes = [1e-2, 1e-3]
tolerances = [(t,t) for t in [1e-4, 1e-5, 1e-6, 1e-7, 1e-8 ]]
unknowns = (100,)
alphas = (2e-1,)
sparsity = (0,1)
 
def GenerateRunList():
	stepArgs = [ {"dt": s} for s in stepsizes ]
	tolArgs = [ {"dt": 0.01, "atol": t[0], "rtol": t[1]} for t in tolerances ]

	runlist = []
	for N, m, ivp, alpha, sparse in itertools.product(unknowns, methods, ivps, alphas, sparsity):
		for ms in method_solvers[m]:
			for tol in stepArgs if ms == 'ConstantSolver' else tolArgs:
				argdict = {"ivp": ivp,
						   "method": m,
						   "solver": ms,
						   "N": N,
						   "alpha": alpha,
						   "sparse": sparse,
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

def AccuracyLegendName(runinfo):	
	common = runinfo[0][1]['method name'] + ' ' + ('Sparse' if int(runinfo[0][1]['sparse']) else 'Dense')
	if runinfo[0][1]['method'] in ark_methods:
		common += ' (' + ('Jacobian Splitting' if 'jacobian splitting' in runinfo[0][1] and int(runinfo[0][1]['jacobian splitting']) else 'Physical Splitting') + ')'
	return common

def GenerateAnalysisPasses():
	passes = []
	for N, ivp, alpha in itertools.product(unknowns, ivps, alphas):
		accuracy = { 'mode': 'Accuracy',
					 'legend': AccuracyLegendName,
					 'reference run': {'ivp':    ivp,
									   'N':      N,
									   'method': 'RODAS',
									   'alpha':  alpha,
									   'sparse': 1,
									   'atol':   tolerances[-1][0] },
					 'xsize': 7,
					 'ysize': 3,
					 'match': {'ivp': ivp, 'N': N, 'alpha': alpha },
					 'group': ['method','solver','jacobian splitting','sparse'] }

		metrics = {"time": "CPU Time (ms)", "steps": "Steps"}		
		for m in metrics:
			title = ivp + ' -- ' + metrics[m] + ' vs. Accuracy' + \
				(" (N=%d)" % (N,))
			filename = ivp.lower() + '-%s-%d-%.2e' % (m, N, alpha)
			accuracy.update({ 'title': title, 'filename': filename, 'comparison': m, "xlabel": "Accuracy", "ylabel": metrics[m]})
			passes.append(accuracy.copy())
	
	return passes

