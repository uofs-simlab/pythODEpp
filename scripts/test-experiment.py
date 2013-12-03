import os

simname = 'test-equations'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

ivps    = ['NonstiffA1','NonstiffA2','NonstiffA3','NonstiffA4','NonstiffA5',
           'NonstiffB1','NonstiffB2','NonstiffB3','NonstiffB4','NonstiffB5',
           'NonstiffC1','NonstiffC2','NonstiffC3','NonstiffC4','NonstiffC5',
           'NonstiffD1','NonstiffD2','NonstiffD3','NonstiffD4','NonstiffD5',
           'NonstiffE1','NonstiffE2','NonstiffE3','NonstiffE4','NonstiffE5',
           'NonstiffF1','NonstiffF2','NonstiffF3','NonstiffF4','NonstiffF5']

method_solvers = { 'ForwardEuler':  ['ConstantSolver','StepDoublingSolver'],
				   'BackwardEuler': ['ConstantSolver','StepDoublingSolver'],
				   'RK4':           ['ConstantSolver','StepDoublingSolver'],
				   'RODAS':         ['EmbeddedSolver','StepDoublingSolver'],
				   'Radau5':        ['EmbeddedSolver'],
				   'ARK3':          ['EmbeddedSolver','StepDoublingSolver'],
				   'ARK4':          ['EmbeddedSolver','StepDoublingSolver'],
				   'ARK5':          ['EmbeddedSolver','StepDoublingSolver'],
				   'RKC1':          ['EmbeddedSolver'],
				   'RKC2':          ['EmbeddedSolver'],
				   'PRKC':          ['EmbeddedSolver'],
				   'RKF45':         ['EmbeddedSolver','StepDoublingSolver'],
				   'DOPR54':        ['EmbeddedSolver','StepDoublingSolver'] }
methods = [ m for m in method_solvers ]

reference_solutions = {'NonstiffA1': {'t': 20., 'y': [ 2.061153353012535e-09] },
					   'NonstiffA2': {'t': 20., 'y': [ 2.182178902359887e-01] },
					   'NonstiffA3': {'t': 20., 'y': [ 2.491650271850414e+00] },
					   'NonstiffA4': {'t': 20., 'y': [ 1.773016648131483e+01] },
					   'NonstiffA5': {'t': 20., 'y': [-7.887826688964196e-01] },
					   'NonstiffB1': {'t': 20., 'y': [ 6.761876008576667e-01, 1.860816099640036e-01] },
					   'NonstiffB2': {'t': 20., 'y': [ 1.000000001030576e+00, 1.000000000000000e+00,
													   9.999999989694235e-01] },
					   'NonstiffB3': {'t': 20., 'y': [ 2.061153488557776e-09, 5.257228022048349e-02,
													   9.474277177183630e-01] },
					   'NonstiffB4': {'t': 20., 'y': [ 9.826950928005993e-01, 2.198447081694832e+00,
													   9.129452507276399e-01] },
					   'NonstiffB5': {'t': 20., 'y': [-9.396570798729192e-01, -3.421177754000779e-01,
													   7.414126596199957e-01] },
					   'NonstiffC1': {'t': 20., 'y': [ 2.061153622240064e-09, 4.122307244619555e-08,
                                                       4.122307244716968e-07, 2.748204829855288e-06,
                                                       1.374102414941961e-05, 5.496409659803266e-05,
                                                       1.832136553274552e-04, 5.234675866508716e-04,
                                                       1.308668966628220e-03,
                                                       9.979127409508656e-01] },
					   'NonstiffC2': {'t': 20., 'y': [ 2.061153577984930e-09, 2.061153573736588e-09,
                                                       2.061153569488245e-09, 2.061153565239902e-09,
                                                       2.061153560991560e-09, 2.061153556743217e-09,
                                                       2.061153552494874e-09, 2.061153548246532e-09,
                                                       2.061153543998189e-09, 9.999999814496180e-01] },
                       'NonstiffC3': {'t': 20., 'y': [ 2.948119211022058e-03, 5.635380154844266e-03,
                                                       7.829072515926013e-03, 9.348257908594937e-03,
                                                       1.007943610301970e-02, 9.982674171429909e-03,
                                                       9.088693332766085e-03, 7.489115195185912e-03,
                                                       5.322964130953349e-03, 2.762434379029886e-03] },
                       'NonstiffC4': {'t': 20., 'y': [ 3.124111453721466e-03, 6.015416842150318e-03,
                                                       8.470021834842650e-03, 1.033682931733337e-02,
                                                       1.153249572873923e-02, 1.204549525737964e-02,
                                                       1.192957068015293e-02, 1.128883207111195e-02,
                                                       1.025804501391024e-02, 8.982017581934167e-03,
                                                       7.597500902492453e-03, 6.219920556824985e-03,
                                                       4.935916341009131e-03, 3.801432544256119e-03,
                                                       2.844213677587894e-03, 2.069123394222672e-03,
                                                       1.464687282843915e-03, 1.009545263941126e-03,
                                                       6.779354330227017e-04, 4.437815269118510e-04,
                                                       2.833264542938954e-04, 1.765005798796805e-04,
                                                       1.073342592697238e-04, 6.374497601777217e-05,
                                                       3.698645309704183e-05, 2.097466832643746e-05,
                                                       1.162956710412555e-05, 6.306710405783322e-06,
                                                       3.346286430868515e-06, 1.737760074184334e-06,
                                                       8.835366904275847e-07, 4.399520411127637e-07,
                                                       2.146181897152360e-07, 1.025981211654928e-07,
                                                       4.807864068784215e-08, 2.209175152474847e-08,
                                                       9.956251263138180e-09, 4.402193653748924e-09,
                                                       1.910149382204028e-09, 8.135892921473050e-10,
                                                       3.402477118549235e-10, 1.397485617545782e-10,
                                                       5.638575303049199e-11, 2.235459707956947e-11,
                                                       8.710498036398032e-12, 3.336554275346643e-12,
                                                       1.256679567784939e-12, 4.654359053128788e-13,
                                                       1.693559145599857e-13, 5.996593816663054e-14,
                                                       1.891330702629865e-14] },
                       'NonstiffC5': {'t': 20., 'y': [ -4.792730224323733e+00, -2.420550725448973e+00,
                                                       -9.212509306014886e-01, -4.217310404035213e+00,
                                                        7.356202947498970e+00,  3.223785985421212e+00,
                                                        4.035559443262270e+00,  1.719865528670555e+01,
                                                        7.478910794233703e+00, -2.998759326324844e+01,
                                                       -4.107310937550929e+00, -9.277008321754407e-01,
                                                       -2.442125302518482e+01,  2.381459045746554e+01,
                                                        1.492096306951359e+01,  3.499208963063806e-01,
                                                       -5.748487687912825e-01, -2.551694020879149e-01,
                                                       -5.237040978903326e-01, -2.493000463579661e-01,
                                                       -8.045341642044464e-02, -3.875289237334110e-01,
                                                        5.648603288767891e-02,  3.023606472143342e-02,
                                                        4.133856546712445e-02, -2.862393029841379e-01,
                                                       -1.183032405136207e-01, -1.511986457359206e-01,
                                                       -2.460068894318766e-01, -3.189687411323877e-02] },
                       'NonstiffD1': {'t': 20., 'y': [  2.198835352008397e-01,  9.427076846341813e-01,
													   -9.787659841058176e-01,  3.287977990962036e-01] },
                       'NonstiffD2': {'t': 20., 'y': [ -1.777027357140412e-01,  9.467784719905892e-01,
													   -1.030294163192969e+00,  1.211074890053952e-01] },
                       'NonstiffD3': {'t': 20., 'y': [ -5.780432953035361e-01,  8.633840009194193e-01,
													   -9.595083730380727e-01, -6.504915126712089e-02] },
                       'NonstiffD4': {'t': 20., 'y': [ -9.538990293416394e-01,  6.907409024219432e-01,
													   -8.212674270877433e-01, -1.539574259125825e-01] },
                       'NonstiffD5': {'t': 20., 'y': [ -1.295266250987574e+00,  4.003938963792321e-01,
													   -6.775390924707566e-01, -1.270838154278686e-01] },
                       'NonstiffE1': {'t': 20., 'y': [  1.456723600728308e-01, -9.883500195574063e-02] },
                       'NonstiffE2': {'t': 20., 'y': [  2.008149762174948e+00, -4.250887527320057e-02] },
                       'NonstiffE3': {'t': 20., 'y': [ -1.004178858647128e-01,  2.411400132095954e-01] },
                       'NonstiffE4': {'t': 20., 'y': [  3.395091444646555e+01,  2.767822659672869e-01] },
                       'NonstiffE5': {'t': 20., 'y': [  1.411797390542629e+01,  2.400000000000002e+00] },
                       'NonstiffF1': {'t': 20., 'y': [ -1.294460621213470e1,   -2.208575158908672e-15 ] },
                       'NonstiffF2': {'t': 20., 'y': [ 70.03731057008607e0 ] },
                       'NonstiffF3': {'t': 20., 'y': [ -3.726957553088175e-1,  -6.230137949234190e-1] },
                       'NonstiffF4': {'t': 20., 'y': [  9.815017249707434e-11] },
                       'NonstiffF5': {'t': 20., 'y': [  1.0 ] }}

stepsizes = [1e-1, 1e-2, 1e-3, 1e-4 ]
tolerances = [(t,t) for t in [1e-4, 1e-5, 1e-6, 1e-7, 1e-8]]

def GenerateRunList():
	stepArgs = [ {"dt": s} for s in stepsizes ]
	tolArgs = [ {"dt": 0.01, "atol": t[0], "rtol": t[1], "newton tol": t[1]*1e-2} for t in tolerances ]

	runlist = []

	for ivp in ivps:
		for m in methods:
			for ms in method_solvers[m]:
				for tol in stepArgs if ms == 'ConstantSolver' else tolArgs:
						argdict = {"ivp": ivp,
								   "method": m,
								   "solver": ms,
								   "jacobian splitting": 1,
								   "jacobian": 'Forward',
								   "max steps": 500000,
								   "min write time": 0.1 }
						argdict.update(tol)
						runlist.append(argdict)
	
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
		tol = ", h=" + str(runinfo['dt'])
	else:
		tol = ", tol=" + str(runinfo['rtol'])

	return runinfo['method name'] + " - " + runinfo['solver'] + tol

def AccuracyLegendName(runinfo):
	info = runinfo[0][1]
	return info['method name'] + '(' + info['solver name'] + ')'

def GenerateAnalysisPasses():
	passes = []
	for ivp in ivps:
		passes.append({ 'mode': 'Solutions',
						'title': ivp + ' Solutions',
						'filename': ivp + '-solutions',
						'xlabel': 'Time (s)',
						'ylabel': 'Solution',
						'legend': SolutionLegendName,
						'match': {'ivp': ivp, 'method': methods, 'atol': tolerances[-1][0] } })
		passes.append({ 'mode': 'Accuracy',
						'title': ivp + ' CPU Time vs. Accuracy',
						'filename': ivp + '-cputime',
						'xlabel': 'Accuracy',
						'ylabel': 'CPU Time (ms)',
						'legend': AccuracyLegendName,
						'reference solution': reference_solutions[ivp],
						'match': {'ivp': ivp},
						'comparison': 'time',
						'group': ['method','solver','jacobian'] })
		passes.append({ 'mode': 'Accuracy',
						'title': ivp + ' Steps vs. Accuracy',
						'filename': ivp + '-steps',
						'xlabel': 'Accuracy',
						'ylabel': 'Steps',
						'legend': AccuracyLegendName,
						'reference solution': reference_solutions[ivp],
						'match': {'ivp': ivp},
						'comparison': 'steps',
						'group': ['method','solver','jacobian'] })
	return passes

