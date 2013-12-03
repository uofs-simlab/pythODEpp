import os
import numpy as np

simname = 'thesis-plots'
simpath = os.path.join(os.getenv('HOME'), 'tmp', simname)

def GenerateRunList():
	runs = []
	runs.append({"ivp": 'AdvectionDiffusion1D',
				 "method": 'RK4',
				 "solver": 'ConstantSolver',
			     "dt": 1e-6,
			     "N": 100, "adv": 1e2, "diff": 1e1,
				 "min write time": 1e-4,
			     "jacobian": "Forward" })
	runs.append({"ivp": 'Brusselator1D',
				 "method": 'RK4',
				 "solver": 'ConstantSolver',
			     "dt": 1e-4,
			     "N": 40,
				 "min write time": 1e-1,
			     "jacobian": "Forward" })	
	runs.append({"ivp": 'Brusselator2D',
				 "method": 'Radau5',
				 "solver": 'EmbeddedSolver',
				 "jacobian": 'Analytic',
				 "print time": 1,
				 "sparse": 1,
				 "atol": 1e-7,
				 "rtol": 1e-7,
			     "N": 30,
				 "min write time": 1e-2 })
	for iso in (0,1):
		for bc in (0,1):
			runs.append({"ivp": 'ConcreteRewetting',
						 "method": 'RKC2',
						 "solver": 'EmbeddedSolver',
					     "N": 100,
						 "sink bc": bc,
						 "isopropanol": iso,
						 "min write time": 1e-1 })
	return runs

def SolutionLegendName(runinfo):
	return runinfo['method name']

def GenerateAnalysisPasses():
	passes = []
	passes.append({ 'mode': 'Solution1D',
		  'match': {'ivp': 'AdvectionDiffusion1D' },
		  'title': "Advection 1D Solution",
		  'filename': "advectiondiffusion1d",
		  'plottxt': 1,
		  'xlabel': 'N',
		  'ylabel': 'Value',
		  'xsize': 6,
		  'ysize': 4,
		  'solution times': [ t*2e-3 for t in range(6)] })
	for iso in (0,1):
		for bc in (0,1):
			passes.append({ 'mode': 'Solution1D',
				  'match': {'ivp': 'ConcreteRewetting', 'sink bc': bc, 'isopropanol': iso },
				  'title': "%s Boundary - %s Solution" % ("Sink" if bc else "Insulated", "Isopropanol" if iso else "Water"),
				  'filename': "concrete-bc%d-iso%d" % (bc,iso),
				  'plottxt': 1,
				  'xlabel': 'Height (cm)',
				  'ylabel': 'C_3S Concentration (g/cm^3)',
				  'xmin': 0,
				  'xmax': 7,
				  'xsize': 6,
				  'ysize': 4,
				  'solution count': 70,
				  'solution times': [ t for t in np.linspace(0,28,10) ] })

	for offset in (0,1):
		passes.append({ 'mode': 'Solutions',
			  'match': {'ivp': 'Brusselator1D' },
			  'title': "Brusselator 1D solution of %s" % ("v" if offset else "u"),
			  'filename': "brusselator1d-%s" % ("v" if offset else "u"),
			  'plottxt': 1,
			  'xlabel': 'Time',
			  'ylabel': 'x',
			  'zlabel': ("v" if offset else "u"),
			  'ymin': 0,
			  'ymax': 1,
			  'ytics': 0.2,
			  'surface': 1,
			  'xsize': 5,
			  'ysize': 3.5,
			  'solution stride': 2,
			  'solution offset': offset })
		for t in range(12):
			passes.append({ 'mode': 'Solution2D',
				  'match': {'ivp': 'Brusselator2D' },
				  'title': "t=%d" % int(t),
				  'filename': "bruss2d-%s-t%.1f" % ("v" if offset else "u",t),
				  'plottxt': 1,
				  'plotpdf': 0,
				  'xlabel': 'x',
				  'ylabel': 'y',
				  'xsize': 6,
				  'ysize': 4,
				  'xdim': 30,
				  'ydim': 30,
				  'cbmin': 0,
				  'cbmax': 4,
				  'solution time': t,
				  'solution stride': 2,
				  'solution offset': offset })
	return passes

