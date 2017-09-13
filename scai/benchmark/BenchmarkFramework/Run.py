import os
import sys
import subprocess
import signal
import commands
import shlex
import time
import select
import platform

from Exceptions import Interrupt

def benchmark( envar, param, verbose ):

	dir = os.path.realpath(os.path.dirname(sys.argv[0]))

	executable = dir + "/RunBenchmark"
 
        # In build directory we will find RunBenchmark in src directory

	if not os.path.exists( executable ):
		executable = os.path.normpath(dir + "/src/RunBenchmark")

	os_cmd = str( envar )+ " " + executable + " " + str( param )
	
	#TODO: Add verbosity flag to print os_cmd in verbose mode

	if verbose:
		print("Issued command: %s"%os_cmd)

	pro = 0
	if platform.system() != 'Windows':
		pro = subprocess.Popen( os_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,close_fds=True )
	else:
		pro = subprocess.Popen( os_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,close_fds=False )
	
	ostream = pro.stdout
	estream = pro.stderr
	istream = sys.stdin

	fd1 = ostream.fileno( )
	fd2 = estream.fileno( )
	fd3 = istream.fileno( )

	# Create a poll to ask for any input of the three streams

	poll = select.poll()
	poll.register(fd1, select.POLLIN)
	poll.register(fd2, select.POLLIN)
	poll.register(fd3, select.POLLIN)

	ostarted = False
	eerror   = False
	ewarning = False
	map      = { 'ERROR': False,'WARNING':False,'MESSAGE': '','MINEX': False }
	i        = 0
	num_omes = 0
	num_emes = 0

	# We stop if output and error stream both have a hangup
	# otherwise we might miss some lines and get bad output

	stop1 = False
	stop2 = False

	while not stop1 or not stop2:

	   events = poll.poll(0)   # just ask for an event, do not wait for it

	   # in case of no event sleep some time to avoid CPU overhead due to poll

	   if len(events) == 0: time.sleep(0.3)

	   for event in events:

		   (fd, ev) = event

		   if ev == select.POLLHUP and pro.poll() != None:

			  if fd == fd1: stop1 = True
			  if fd == fd2: stop2 = True

		   elif fd == fd3:

			  iline = istream.readline( )

			  # ^F used to interrupt current benchmark and skip to next

			  if iline.startswith("\x06"):
				 if pro.poll( ) is None:
					os.kill( pro.pid,signal.SIGTERM )
					raise Interrupt( )

		   elif fd == fd1:
			  oline = ostream.readline( ).rstrip( )
			  # parse outputstream
			  if oline == '%_BENCHMARK_FRAMEWORK_!MESSAGE_START_!':
				  ostarted = True
				  num_omes = num_omes + 1
			  elif oline == '%_BENCHMARK_FRAMEWORK_!MESSAGE_END_!':
				  ostarted = False
				  num_omes = num_omes + 1
			  elif ostarted:
				  values = oline.split( '%,' )
				  map['NAME']        = values[0]
				  map['ISET']        = values[1]
				  map['GROUP']       = values[2]
				  map['THREADS']     = values[3]
				  map['VTS']         = values[4]
				  map['FLOPS']       = values[5]
				  map['BANDWIDTH']   = values[6]
				  map['SETUP']       = values[7]
				  map['MINEX']       = values[8]
				  map['MAXEX']       = values[9]
				  map['EXECUTION']   = values[10]
				  map['TEARDOWN']    = values[11]
				  map['TOTEX']       = values[12]
				  map['EXECUTED']    = values[13]
				  map['REPETITIONS'] = values[14]
			  elif oline != '':
				  print oline
		   elif fd == fd2:
			  # parse error stream
			  eline = estream.readline( ).rstrip( )
			  if eline == '%_BENCHMARK_FRAMEWORK_!ERROR_START_!':
					eerror = True
					num_emes = num_emes + 1
			  elif eline == '%_BENCHMARK_FRAMEWORK_!ERROR_END_!':
					eerror = False
					num_emes = num_emes + 1
			  elif eline == '%_BENCHMARK_FRAMEWORK_!WARNING_START_!':
					ewarning = True
					num_emes = num_emes + 1
			  elif eline == '%_BENCHMARK_FRAMEWORK_!WARNING_END_!':
					ewarning = False
					num_emes = num_emes + 1
			  elif eerror:
					map['ERROR']   = True
					map['MESSAGE'] += i*'\n' + eline
					i = 1
			  elif ewarning:
					map['WARNING'] = True
					map['MESSAGE'] += i*'\n' + eline
					i = 1
			  elif eline != '':
					print eline
 
	ostream.close( )
	estream.close( )

	if not map['MINEX']:
		map['WARNING'] = True
		map['MESSAGE'] = "Uncaught Exception.\n"

	#TODO: Add verbosity flag to print this only in verbose mode
	if not map['WARNING'] and not map['ERROR']:
		print "  Minimum Execution Time = ", map['MINEX']
	return map

def getInputSets( verbose ):
	dir = os.path.realpath(os.path.dirname(sys.argv[0]))
	
	executable = dir + "/RegisteredInputSets"
 
        # In build directory we will find RegisteredInputSets in src directory

	if not os.path.exists(executable):
		executable = os.path.normpath(dir + "/src/RegisteredInputSets")

        os_cmd = executable  # no arguments

	if verbose:
            print("Issued command: %s"%os_cmd)

	pro = 0
	if platform.system() != 'Windows':
		pro = subprocess.Popen( os_cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,close_fds=True )
	else:
		pro = subprocess.Popen( os_cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,close_fds=False )

	ostream = pro.stdout
	estream = pro.stderr

	started = False
	line    = ostream.readline( )
	map     = {}
	iSets   = []

	# parse outputstream
	while line != '':
		line = line.rstrip( )
		if line == '%_BENCHMARK_FRAMEWORK_!MESSAGE_START_!':
			started = True
		elif line == '%_BENCHMARK_FRAMEWORK_!MESSAGE_END_!':
			started = False
		elif started:
			values = line.split( '%,' )
			for v in values:
				iSets.append( v )
		else:
			print line
		line = ostream.readline( )
	ostream.close( )

	map['ISETS']   = iSets

	error   = False
	warning = False

	line = estream.readline( )
	map['MESSAGE'] = ''
	map['ERROR']   = False
	map['WARNING'] = False

	i = 0

	# parse errorstream
	while line != '':
		line = line.rstrip( )
		if line == '%_BENCHMARK_FRAMEWORK_!ERROR_START_!':
			error   = True
		elif line == '%_BENCHMARK_FRAMEWORK_!ERROR_END_!':
			erro    = False
		elif line == '%_BENCHMARK_FRAMEWORK_!WARNING_START_!':
			warning = True
		elif line == '%_BENCHMARK_FRAMEWORK_!WARNING_END_!':
			warning = False
		elif error:
			map['ERROR']    = True
			map['MESSAGE'] += i*'\n' + line
			i = 1
		elif warning:
			map['WARNING']  = True
			map['MESSAGE'] += i*'\n' + line
			i = 1
		else:
			print line
		line = estream.readline( )
	estream.close( )

	return map

def getBenchmarks( verbose ):

	dir = os.path.realpath(os.path.dirname(sys.argv[0]))

	executable = dir + "/RegisteredBenchmarks"
 
        # In build directory we will find RunBenchmark in src directory

	if not os.path.exists(executable):
		executable = os.path.normpath(dir + "/src/RegisteredBenchmarks")

        os_cmd = executable  # no arguments

	if verbose:
		print( "Issued command = %s"%os_cmd )
		
	pro = 0
	if platform.system() != 'Windows':
		pro = subprocess.Popen( os_cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,close_fds=True )
	else:
		pro = subprocess.Popen( os_cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,close_fds=False )

	ostream = pro.stdout
	estream = pro.stderr

	started    = False
	line       = ostream.readline( )
	map        = {}
	benchmarks = {}

	# parse outputstream
	while line != '':
		line = line.rstrip( )
		if line == '%_BENCHMARK_FRAMEWORK_!MESSAGE_START_!':
			started = True
		elif line == '%_BENCHMARK_FRAMEWORK_!MESSAGE_END_!':
			started = False
		elif started:
			values = line.split( '%,' )
			for v in values:
				key = v.split( '%_:_%' )[0]
				val = v.split( '%_:_%' )[1]
				benchmarks[key] = val
		else:
			print line
		line = ostream.readline( )
	ostream.close( )

	map['BENCHMARKS']   = benchmarks

	error   = False
	warning = False

	line = estream.readline( )
	map['MESSAGE'] = ''
	map['ERROR']   = False
	map['WARNING'] = False

	i = 0

	# parse errorstream
	while line != '':
		line = line.rstrip( )
		if line == '%_BENCHMARK_FRAMEWORK_!ERROR_START_!':
			error   = True
		elif line == '%_BENCHMARK_FRAMEWORK_!ERROR_END_!':
			erro    = False
		elif line == '%_BENCHMARK_FRAMEWORK_!WARNING_START_!':
			warning = True
		elif line == '%_BENCHMARK_FRAMEWORK_!WARNING_END_!':
			warning = False
		elif error:
			map['ERROR']    = True
			map['MESSAGE'] += i*'\n' + line
			i = 1
		elif warning:
			map['WARNING']  = True
			map['MESSAGE'] += i*'\n' + line
			i = 1
		else:
			print line
		line = estream.readline( )

	estream.close( )

	return map
