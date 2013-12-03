import os

simname = 'cell-model'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

ivps = ['CellModel']
method = 'ForwardEuler'
solver = 'ConstantSolver'

def GenerateRunList():
	argdict = {"ivp": "CellModel",
			   "method":     method,
			   "solver":     solver,
			   "dt":         5e-4,
			   "tf":         5e-2,
			   "rtol":       1e-5,
			   "atol":       1e-5,
			   "max steps":  1000000}
	return [argdict]

def SolutionLegendName(runinfo):
	return runinfo['method name']

def AccuracyLegendName(runinfo):
	info = runinfo[0][1]
	return info['method name'] + '(' + info['solver name'] + ')'

def GenerateAnalysisPasses():
	passes = []
	for ivp in ivps:
		passes.append({ 'mode': 'Solutions',
						'title': ivp + ' Solutions',
						'xlabel': 'Time (s)',
						'ylabel': 'Solution',
						'legend': SolutionLegendName,
						'solnames': ['[Na+]', '[K+]', '[Cl-]' ],
						'match': {'ivp': ivp, 'method': method} })
	return passes

