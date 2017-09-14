import os
import sys
import subprocess
import signal
import fcntl
import commands
import popen2
import shlex
import time
import select

from Exceptions import Interrupt

def benchmark( envar, param ):
	dir = os.path.realpath(os.path.dirname(sys.argv[0]))
	os_cmd = str( envar )+ " " + dir + "/src/BenchmarkRunner " + str( param )

	pro = subprocess.Popen( os_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,close_fds=True )

	ostream = pro.stdout
	estream = pro.stderr
	istream = sys.stdin

        poll = select.poll()

	fd1 = ostream.fileno( )
	fd2 = estream.fileno( )
        fd3 = istream.fileno( )

        poll.register(fd1, select.POLLIN)
        poll.register(fd2, select.POLLIN)

        ostarted = False
        eerror   = False
        ewarning = False
        map      = { 'ERROR': False,'WARNING':False,'MESSAGE': '' }
        i        = 0
        num_omes = 0
        num_emes = 0

        stop = False

        while not stop:

           events = poll.poll(0)

           if len(events) == 0:
              time.sleep(0.3)

           for event in events:

               fd = event[0]
               ev = event[1]
 
               if ev == select.POLLHUP and pro.poll() != None:
                  print 'Finished'
                  stop = True
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
                  eline = estream.readline( ).rstrip( )
                  print 'eline = ', eline
 
	# parse outputstream, while process is running.
	while pro.poll( ) is None:
		# sleep is necessary because otherwise the polling needs to much
		# cpu resources and distrubes the benchmarks
		time.sleep(1.0)
		try:
			iline = istream.readline( ).rstrip( )
			thrown = False
		except:
			thrown = True
		if not thrown:
			if iline=='':
				if pro.poll( ) is None:
					os.kill( pro.pid,signal.SIGTERM )
				raise Interrupt( )
		try:
			oline = ostream.readline( ).rstrip( )
			thrown = False
		except:
			thrown = True
			# do nothing
		if not thrown:
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
		try:
			eline = estream.readline( ).rstrip( )
			thrown = False
		except:
			thrown = True
			# do nothing
		if not thrown:
			# parse error stream
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

	# in case, process ended, but not all output has been read, yet
	if num_omes < 2 and num_emes < 2:
		for eline in pro.stderr.readlines( ):
			eline = eline.rstrip( )
			if eline == '%_BENCHMARK_FRAMEWORK_!ERROR_START_!':
                                eerror = True
                        elif eline == '%_BENCHMARK_FRAMEWORK_!ERROR_END_!':
                                eerror = False
                        elif eline == '%_BENCHMARK_FRAMEWORK_!WARNING_START_!':
                                ewarning = True
                        elif eline == '%_BENCHMARK_FRAMEWORK_!WARNING_END_!':
                                ewarning = False
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
		for oline in pro.stdout.readlines( ):
			oline = oline.rstrip( )
			if oline == '%_BENCHMARK_FRAMEWORK_!MESSAGE_START_!':
                                ostarted = True
                        elif oline == '%_BENCHMARK_FRAMEWORK_!MESSAGE_END_!':
                                ostarted = False
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
	ostream.close( )
	estream.close( )

	return map


def getInputSets( ):
	dir = os.path.realpath(os.path.dirname(sys.argv[0]))
	pro = subprocess.Popen( dir + '/src/RegisteredInputSets',stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,close_fds=True )

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

def getBenchmarks( ):
	dir = os.path.realpath(os.path.dirname(sys.argv[0]))
        pro = subprocess.Popen( dir + '/src/RegisteredBenchmarks',stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,close_fds=True )

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
	
