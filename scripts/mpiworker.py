import os, sys, json, threading, Queue, subprocess
from time import sleep

root_dir = os.path.dirname(os.path.realpath(sys.argv[0]))

busylock = threading.Lock()
completedRuns = []
workerThreadBusy = False

class WorkerThread(threading.Thread):
	def __init__(self,queue):
		threading.Thread.__init__(self)
		self.queue = queue

	def run(self):
		global workerThreadBusy
		while True:
			argdict = self.queue.get()
			arglist = []
			for a in argdict:
				arglist.append('-'+a)
				arglist.append(str(argdict[a]))
			os.chdir(argdict['path'])
			result = subprocess.call([os.path.join(root_dir,'..','pythODE++'),'-phase','runner']+arglist)
			if result != 0:
				print "run failed:", argdict['path']
				print ' '.join(map(lambda x: "'"+x+"'",arglist))
			self.queue.task_done()
			busylock.acquire()
			workerThreadBusy = False
			completedRuns.append(argdict['path'])
			busylock.release()

def MPIWorker(s):
	global workerThreadBusy
	workerQueue = Queue.Queue()
	wt = WorkerThread(workerQueue)
	wt.setDaemon(True)
	wt.start()

	while True:
		message = json.loads(s.recv(4096))
		if message == "kill":
			s.close()
			sys.exit(0)
		elif message == "status":
			busylock.acquire()
			s.send(json.dumps("busy" if workerThreadBusy else "free") + "\0")
			busylock.release()
		elif message == "completed":
			busylock.acquire()
			s.send(json.dumps("|".join(completedRuns)) + "\0")
			busylock.release()
		else:
			workerQueue.put(message)
			busylock.acquire()
			workerThreadBusy = True
			s.send(json.dumps("queued") + "\0")
			busylock.release()

