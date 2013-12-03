#!/usr/bin/env python
import os, sys, threading, datetime, time, \
	random, string, socket, importlib, json

from mpi4py import MPI
comm = MPI.COMM_WORLD

if __name__ != '__main__':
	sys.exit(0)

# -----------------------------------------------------------------------------
#
# KillMPI()
#
# Send kill messages to each of the MPI processes via sockets
#
# -----------------------------------------------------------------------------
def KillMPI():
	global hosts
	print "Sending kill messages to MPI processes"
	for h in hosts:
		hosts[h]['socket'].send(json.dumps("kill"))
		hosts[h]['socket'].close()

# -----------------------------------------------------------------------------
#
# SafeSendRecv(s,msg)
#
# Performs a send and receive to one of the MPI processes. This function
# involves a global lock so that the input thread and the dispatching thread
# cannot garble messages.
#
# -----------------------------------------------------------------------------
sendRecvLock = threading.Lock()
messageRemainder = ""

def RecvUntilZero(s):
	global messageRemainder

	while True:
		msg = messageRemainder + s.recv(1024)
		for i in range(len(msg)):
			if msg[i] == "\0":
				messageRemainder = msg[i+1:]
				msg = msg[:i]
				return msg
		messageRemainder = msg

def SafeSendRecv(s, msg):
	global sendRecvLock
	sendRecvLock.acquire()
	s.send(json.dumps(msg))
	msg = RecvUntilZero(s)
	ret = json.loads(msg)
	sendRecvLock.release()
	return ret

# -----------------------------------------------------------------------------
#
# InputThread
#
# Reads stdin and sends messages to the MPI processes so the user can
# track simulation progress.
#
# -----------------------------------------------------------------------------
class InputThread(threading.Thread):
	def __init__(self):
		threading.Thread.__init__(self)
		self.commands = {
			'help':      [ self.cmdHelp,      'Displays a list of commands and their descrptions'],
			'status':    [ self.cmdStatus,    'Displays what problems are currently being solved'],
			'details':   [ self.cmdStatusDet, 'Displays the details of what problems are currently being solved'],
			'progress':  [ self.cmdProgress,  'Displays the overall simulation progress'] }
						
	def run(self):
		while True:
			line = sys.stdin.readline().rstrip('\n')
			if line not in self.commands:
				line = 'help'
			self.commands[line][0]()

	def cmdHelp(self):
		print 'List of commands'
		print '-------------------'
		for c in self.commands:
			print '%-20s %s' % (c, self.commands[c][1])
		print ''

	def inputArgs(self,q):
		retlist = []
		for e in q:
			retlist.append("'-%s'" % e)
			retlist.append("'%s'" % q[e])
		return ' '.join(retlist)

	def cmdStatusBody(self,details):
		global listLock, queuedlist
		completedRuns = []
		for h in hosts:
			ret = SafeSendRecv(hosts[h]['socket'], "completed")
			ret = ret.split("|") if ret else []
			if len(ret) == 1 and ret[0] == '':
				continue
			completedRuns += ret

		listLock.acquire()
		for q in queuedlist:
			if q['path'] not in completedRuns:
				print (self.inputArgs(q) if details else q['path'])
		listLock.release()
	
	def cmdStatus(self):
		self.cmdStatusBody(0)
	
	def cmdStatusDet(self):
		self.cmdStatusBody(1)

	def cmdProgress(self):
		print str(len(queuedlist)) + "/" + str(problemcount) + " queued for running."

# -----------------------------------------------------------------------------
#
# RandomString
#
# Generates a random string of the specified length, primarily used to
# generate run ids and paths.
#
# -----------------------------------------------------------------------------
def RandomString(length):
	return ''.join([random.choice(string.uppercase + string.digits) for i in range(length)])

# -----------------------------------------------------------------------------
#
# DateAndTimeString
#
# Generates a string containing the date and time, primarily used to
# generate run aths.
#
# -----------------------------------------------------------------------------
def DateAndTimeString():
	dateandtime = datetime.datetime.now()
	timestamp = "%04d-%02d-%02d--%02d-%02d-%02d" % (dateandtime.year,
													dateandtime.month,
													dateandtime.day,
													dateandtime.hour,
													dateandtime.minute,
													dateandtime.second)
	return timestamp

def MakePath(path):
	if not os.path.exists(path):
		os.makedirs(path)

# -----------------------------------------------------------------------------
#
#     * * * * * * * * * * Procedural code starts here * * * * * * * * * *
#

# This program needs at least 2 processes
if comm.size == 1:
	print "More than 1 process is required."
	sys.exit()

if len(sys.argv) != 2:
	if not comm.rank: print "usage: ivprunner modulename"
	sys.exit()

# Set up communication sockets and pass control to workers
if comm.rank != 0:
	s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	s.bind(('',0))
	port = (socket.gethostname(), s.getsockname()[1])
	comm.gather(port,root=0)
	s.listen(1)
	conn, addr = s.accept()
	s.close()

	from mpiworker import MPIWorker
	MPIWorker(conn)

ports = comm.gather(0,root=0)

hosts = { h: {'status': 'free',
			  'hostname': ports[h][0],
			  'port': ports[h][1],
			  'socket': socket.socket(socket.AF_INET, socket.SOCK_STREAM)}
			  for h in range(1,comm.size) }
for h in hosts:
	hosts[h]['socket'].connect((hosts[h]['hostname'],hosts[h]['port']))

# Gather runs from module
try:
	expmodule = importlib.import_module(sys.argv[1])
except ImportError:
	print "module", sys.argv[1], "not found."
	KillMPI()
	sys.exit(1)

simname = expmodule.simname
simpath = expmodule.simpath

# Create run path if it does not exist
simpath = os.path.join(simpath, '--'.join([simname, DateAndTimeString(), RandomString(8)]))
MakePath(simpath)

# Lists of runs that need to be queued and those that have already been queued
listLock = threading.Lock()
runlist = expmodule.GenerateRunList()
queuedlist = []

# Add path for each run
for r in runlist:
	r["simpath"] = simpath
	r['path'] = os.path.join(simpath,"--".join([r['ivp'],r['method'],r['solver'],RandomString(8)]))

problemcount = len(runlist)
print "Running", problemcount, "problems with", len(hosts) ,("processes." if len(hosts) > 1 else "process.")

# Start up the input thread
it = InputThread()
it.setDaemon(True)
it.start()

def CheckFreeHosts():
	global hosts
	allfree = True
	for h in hosts:
		hosts[h]['status'] = SafeSendRecv(hosts[h]['socket'], "status")
		if hosts[h]['status'] == 'busy':
			allfree = False
	return allfree

def GetFreeHost():
	global hosts
	while True:
		CheckFreeHosts()
		for h in hosts:
			if hosts[h]['status'] == 'free':
				return h
		time.sleep(0.5)

while True:
	if not len(runlist):
		break

	h = GetFreeHost()

	listLock.acquire()
	r = runlist.pop()
	queuedlist.append(r)
	listLock.release()

	MakePath(r['path'])
	result = SafeSendRecv(hosts[h]['socket'], r)

while not CheckFreeHosts():
	time.sleep(0.5)

KillMPI()

