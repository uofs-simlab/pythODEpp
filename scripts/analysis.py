#!/usr/bin/env python
import os, sys, subprocess, importlib

if __name__ != '__main__':
	sys.exit(0)

if len(sys.argv) != 3:
	print "usage: ivpanalysis modulename runpath"
	sys.exit(1)

root_dir = os.path.dirname(os.path.realpath(sys.argv[0]))

# Gather runs from module
try:
	expmodule = importlib.import_module(sys.argv[1])
except ImportError:
	print "module", sys.argv[1], "not found."
	sys.exit(1)

if not os.path.exists(expmodule.simpath):
	print expmodule.simpath, "does not exist."
	sys.exit(1)

if sys.argv[2] == "newest":
	filelist = [ f for f in os.listdir(expmodule.simpath) if f[0] != '.' ]
	dirs = map(lambda x: os.path.join(expmodule.simpath,x), filelist)
	simpath = max(dirs, key=lambda x: os.path.getmtime(x))
else:
	simpath = os.path.join(expmodule.simpath,sys.argv[2])

if not os.path.exists(simpath):
	print simpath, "does not exist."
	sys.exit(1)

# Make the plot path
plotpath = os.path.join(simpath, "Plots")
if not os.path.exists(plotpath):
	os.makedirs(plotpath)

def CheckMatch(paramMatch,paramDict):
	for p in paramMatch:
		if p not in paramDict:
			return False

		matchStr = paramMatch[p]
		if type(matchStr) not in (list, tuple):
			matchStr = [matchStr]

		if not any(str(m) == paramDict[p] for m in matchStr):
			return False
	return True

def FindAllMatchingRuns(paramMatch):
	matches = {}
	untimed_matches = []

	searchDirs = map(lambda x: os.path.join(simpath,x),os.listdir(simpath))
	for path in searchDirs:
		# Fill the parameters dictionary
		paramDict = {}
		infopath = os.path.join(path,'.runinfo')
		if not os.path.exists(infopath):
			continue
		infofile = open(infopath)
		for line in infofile:
			line = line.rstrip('\n')
			key, value = line.split(':')
			paramDict[key]=value
		infofile.close()

		# Skip if it doesn't match
		if not CheckMatch(paramMatch, paramDict):
			continue

		if 'timing group' not in paramDict:
			untimed_matches.append((path,paramDict))
		else:
			tg = paramDict['timing group']
			if tg in matches:
				if float(paramDict['time']) < float(matches[tg][1]['time']):
					matches[tg][1]['time'] = paramDict['time']
			else:
				matches[tg] = (path,paramDict)
	return untimed_matches + [ matches[m] for m in matches ]

def BuildMatchGroups(matches,matchNames):
	groups = { }
	for m in matches:
		matchKey = '-'.join([ m[1][mn] if mn in m[1] else '' for mn in matchNames])
		if matchKey not in groups:
			groups[matchKey] = []
		groups[matchKey].append(m)
	return [ groups[k] for k in groups]

def AddArgs(argname, arglist, ap):
	if argname in ap:
		arglist += ['-'+argname, str(ap[argname])]

def RunAnalysisPass(ap, program):
	matches = FindAllMatchingRuns(ap['match'])
	plotfile = os.path.join(plotpath,ap['filename']+'.pdf')

	if ap['mode'] == 'Solutions':
		arglist = ['-mode', 'Solutions',
				   '-title', ap['title'],
				   '-plotfile', plotfile,
				   '-linetype', 'lines',
				   '-number', str(len(matches)) ]

		if 'surface' in ap: arglist += ['-surface', str(ap['surface']) ]

		AddArgs('xsize', arglist, ap)
		AddArgs('ysize', arglist, ap)
		AddArgs('xlabel', arglist, ap)
		AddArgs('ylabel', arglist, ap)
		AddArgs('zlabel', arglist, ap)
		AddArgs('ymin', arglist, ap)
		AddArgs('ymax', arglist, ap)
		AddArgs('ytics', arglist, ap)
		AddArgs('plotpdf', arglist, ap)
		AddArgs('plottxt', arglist, ap)
		AddArgs('cuspplot', arglist, ap)

		if 'solution offset' in ap:
			arglist += ['-solution offset', str(ap['solution offset']) ]
		if 'solution stride' in ap:
			arglist += ['-solution stride', str(ap['solution stride']) ]
		if 'solution count' in ap:
			arglist += ['-solution count', str(ap['solution count']) ]

		if 'solnames' in ap:
			arglist += ['-solnames',str(len(ap['solnames']))]
	
		p = subprocess.Popen(program+arglist, stdin=subprocess.PIPE)

		if 'solnames' in ap:
			for solname in ap['solnames']:
				p.stdin.write(solname + "\n")

		for path, params in matches:
			p.stdin.write(path + "\n")
			if 'legend' in ap:
				p.stdin.write(ap['legend'](params) + "\n")
			else:
				p.stdin.write("Untitled\n")
			
		p.stdin.close()
		p.wait()
	
	elif ap['mode'] == 'Solution1D':
		if len(matches) == 0:
			print "error: no solution found."
			sys.exit(1)
		elif len(matches) > 1:
			print "error:", len(matches), "multiple matching solutions found."
			sys.exit(1)
			
		arglist = ['-mode', 'Solution1D',
				   '-title', ap['title'],
				   '-xlabel', ap['xlabel'],
				   '-ylabel', ap['ylabel'],
				   '-plotfile', plotfile,
				   '-solution times', str(len(ap['solution times'])),
				   '-linetype', 'lines' ]
		AddArgs('xmin', arglist, ap)
		AddArgs('xmax', arglist, ap)
		AddArgs('xsize', arglist, ap)
		AddArgs('ysize', arglist, ap)
		AddArgs('solution offset', arglist, ap)
		AddArgs('solution stride', arglist, ap)
		AddArgs('solution count', arglist, ap)
		AddArgs('plotpdf', arglist, ap)
		AddArgs('plottxt', arglist, ap)

		p = subprocess.Popen(program+arglist, stdin=subprocess.PIPE)
		p.stdin.write(matches[0][0]+"\n")
		for t in ap['solution times']:
			p.stdin.write(str(t)+"\n")
		p.stdin.close()
		p.wait()

	elif ap['mode'] == 'Solution2D':
		if len(matches) == 0:
			print "error: no solution found."
			sys.exit(1)
		elif len(matches) > 1:
			print "error:", len(matches), "multiple matching solutions found."
			sys.exit(1)
			
		arglist = ['-mode', 'Solution2D',
				   '-title', ap['title'],
				   '-xlabel', ap['xlabel'],
				   '-ylabel', ap['ylabel'],
				   '-plotfile', plotfile,
				   '-colormap', '1',
				   '-solution time', str(ap['solution time'])]
		AddArgs('cbmin', arglist, ap)
		AddArgs('cbmax', arglist, ap)
		AddArgs('xdim', arglist, ap)
		AddArgs('ydim', arglist, ap)
		AddArgs('xsize', arglist, ap)
		AddArgs('ysize', arglist, ap)
		AddArgs('solution offset', arglist, ap)
		AddArgs('solution stride', arglist, ap)
		AddArgs('solution count', arglist, ap)
		AddArgs('plotpdf', arglist, ap)
		AddArgs('plottxt', arglist, ap)

		p = subprocess.Popen(program+arglist, stdin=subprocess.PIPE)
		p.stdin.write(matches[0][0]+"\n")
		p.stdin.close()
		p.wait()
	
	elif ap['mode'] == 'Accuracy':
		groups = BuildMatchGroups(matches,ap['group'])

		arglist = ['-mode', 'Accuracy',
				   '-title', ap['title'],
				   '-xlabel', ap['xlabel'],
				   '-ylabel', ap['ylabel'],
				   '-filename', ap['filename'],
				   '-logscale', "1",
				   '-plotfile', plotfile,
				   '-groups', str(len(groups)) ]
		AddArgs('xsize', arglist, ap)
		AddArgs('ysize', arglist, ap)
		AddArgs('plotpdf', arglist, ap)
		AddArgs('plottxt', arglist, ap)

		if 'symbol' in ap:
			arglist += ['-symbol', '1']

		if 'color' in ap:
			arglist += ['-color', '1']

		# Make sure there is one and only one type of reference
		if ('reference run' in ap) == ('reference solution' in ap):
			print 'error: analysis pass must have one and only one of reference run and reference solution.'
			sys.exit(1)

		# Find reference solution
		if 'reference run' in ap:
			reference = FindAllMatchingRuns(ap['reference run'])
			if len(reference) == 0:
				print 'error: no reference solution found.'
				sys.exit(1)
			if len(reference) != 1:
				print 'error:', len(reference), 'multiple matching reference solutions found.'
				sys.exit(1)
			arglist += ['-reference', reference[0][0]]

		p = subprocess.Popen(program+arglist, stdin=subprocess.PIPE)

		# Write reference solution
		if 'reference solution' in ap:
			refsln = ap['reference solution']
		
			p.stdin.write(str(len(refsln['y'])) + "\n")
			p.stdin.write(str(refsln['t']) + "\n")
			for y in refsln['y']:
				p.stdin.write(str(y) + "\n")

		# Sort groups based on their legends (inefficient but who cares!)
		groups.sort(lambda g1, g2: cmp(ap['legend'](g1), ap['legend'](g2)))

		for group in groups:
			p.stdin.write(str(len(group)) + "\n")
			p.stdin.write(ap['legend'](group) + "\n")
			
			if 'symbol' in ap:
				p.stdin.write(str(ap['symbol'](group)) + "\n")
			if 'color' in ap:
				color = ap['color'](group)
				p.stdin.write("%f %f %f\n"%(color[0], color[1], color[2]))

			for path, params in group:
				p.stdin.write(path + "\n")
				p.stdin.write(params[ap['comparison']] + "\n")
		p.stdin.close()
		p.wait()

	else:
		print "Mode", mode, "is not implemented"


passes = expmodule.GenerateAnalysisPasses()
for i in range(len(passes)):
	RunAnalysisPass(passes[i],["../pythonde++","-phase","gnuplot1d"])

