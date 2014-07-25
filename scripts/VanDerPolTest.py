import os

simname = 'VanDerPol'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

methods = {	'Radau5' : 1,
			'RKC1' : 2,
			'RKC2' : 3,
			'DOPR54' : 4,
			'ARK3' : 5,
			'ARK4' : 6,
			'ARK5' : 7,
			'BS54' : 8 }
ivp = 'VanDerPol'
tol = [(10**-t, 10**-t) for t in range(5, 11)]
tf = 0.55139
timingRuns = 25
epsilon = 0.01
y0 = 2
z0 = -0.6654321

reference = {	'ivp' : ivp,
				'epsilon' : epsilon,
				'y0' : y0,
				'z0' : z0,
				'tf' : tf,
				'method' : 'Radau5',
				'solver' : 'EmbeddedSolver',
				'min write time' : tf / 100,
				'atol' : tol[-1][0] / 1000,
				'rtol' : tol[-1][1] / 1000	}

def GenerateRunList():
	runlist = []
	for m in methods:
		for t in range(len(tol)):
			for x in range(timingRuns):
				ne = {	'ivp' : ivp,
						'epsilon' : epsilon,
						'y0' : y0,
						'z0' : z0,
						'tf' : tf,
						'method' : m,
						'solver' : 'EmbeddedSolver',
						'atol' : tol[t][0],
						'rtol' : tol[t][1],
						'timing group' : methods[m] * 100 + t	}
				runlist.append(ne.copy())
	runlist.append(reference)
	return runlist

def LegendName(runinfo):
	return runinfo[0][1]['method name']

def GetColour(runinfo):
	method = runinfo[0][1]['method']
	if method == 'Radau5':
		return (0.2, 0.8, 0.6)
	if method == 'RKC1':
		return (0.2, 0.2, 0.4)
	if method == 'RKC2':
		return (0.4, 0.8, 1.0)
	if method == 'DOPR54':
		return (0.4, 0.8, 0.0)
	if method == 'ARK3':
		return (0.2, 0.0, 0.2)
	if method == 'ARK4':
		return (0.8, 0.6, 0.2)
	if method == 'ARK5':
		return (0.0, 0.0, 0.6)
	if method == 'BS54':
		return (0.6, 0.4, 0.8)

def GetSymbol(runinfo):
	method = runinfo[0][1]['method']
	if method == 'Radau5':
		return 5
	if method == 'RKC1':
		return 7
	if method == 'RKC2':
		return 8
	if method == 'DOPR54':
		return 9
	if method == 'ARK3':
		return 6
	if method == 'ARK4':
		return 4
	if method == 'ARK5':
		return 10
	if method == 'BS54':
		return 3

def GenerateAnalysisPasses():
	passes = []
	accuracy = {	'mode' : 'Accuracy',
					'comparison' : 'steps',
					'title' : 'Method Accuracy vs Steps',
					'symbol' : GetSymbol,
					'filename' : 'acc_vs_steps',
					'color' : GetColour,
					'xlabel' : 'Error Norm',
					'xsize' : 10.0,
					'ysize' : 5.0,
					'ylabel' : 'Steps',
					'legend' : LegendName,
					'plottxt' : 1,
					'reference run' : {	'ivp' : ivp,
										'method' : reference['method'],
										'atol' : reference['atol']	},
					'match' : {	'atol' : [t[0] for t in tol] },
					'group' : [ 'method' ] }
	time = {	'mode' : 'Accuracy',
				'comparison' : 'time',
				'title' : 'Method Accuracy vs Time',
				'symbol' : GetSymbol,
				'filename' : 'acc_vs_time',
				'color' : GetColour,
				'xlabel' : 'Error Norm',
				'xsize' : 10.0,
				'ysize' : 5.0,
				'ylabel' : 'Time (ms)',
				'legend' : LegendName,
				'plottxt' : 1,
				'reference run' : {	'ivp' : ivp,
									'method' : reference['method'],
									'atol' : reference['atol']	},
				'match' : {	'atol' : [t[0] for t in tol] },
				'group' : [ 'method' ] }
	passes.append(accuracy)
	passes.append(time)
	return passes

