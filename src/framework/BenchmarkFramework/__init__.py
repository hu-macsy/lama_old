# @file __init__.py
#
# @license
# Copyright (c) 2011
# Fraunhofer Institute for Algorithms and Scientific Computing SCAI
# for Fraunhofer-Gesellschaft
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# @endlicense
#
# @brief Module initialization for Benchmark Framework
# @author robin
# @date 06.04.2011
# $Id$

import Options
import Configure
import Column
import Run
import BenchmarkOut
from Exceptions import Interrupt

import texttable

try:
    import cPickle as pickle
except:
    import pickle
import pprint

import sys
import os
import socket
import glob
import getpass
import time
import datetime
import math
import tempfile

MAP_NAME_LIST = ['EXECUTION','MINEX','MAXEX']

class colors:
	OKPINK  = '\033[95m'
	OKGREEN = '\033[32m'
	MESSAGE = '\033[94m'
	WARNING = '\033[33m'
	FAIL    = '\033[31m'
	END     = '\033[0m'

	def __init__( self ):
		if os.name != 'posix':
			colors.OKPINK  = ''
			colors.OKGREEN = ''
			colors.MESSAGE = ''
			colors.WARNING = ''
			colors.FAIL    = ''
			colors.END     = ''

def getColumnMap( ):
	columnMap = {}
	columnMap['n']     = Column.ColumnIdName( )
	columnMap['e']     = Column.ColumnExecutionTime( Column.ColumnTime.TimeUnit.SECONDS )
	columnMap['em']    = Column.ColumnExecutionTime( Column.ColumnTime.TimeUnit.MILLISECONDS )
	columnMap['emax']  = Column.ColumnMaxExecutionTime( Column.ColumnTime.TimeUnit.SECONDS )
        columnMap['emaxm'] = Column.ColumnMaxExecutionTime( Column.ColumnTime.TimeUnit.MILLISECONDS )
        columnMap['emin']  = Column.ColumnMinExecutionTime( Column.ColumnTime.TimeUnit.SECONDS )
        columnMap['eminm'] = Column.ColumnMinExecutionTime( Column.ColumnTime.TimeUnit.MILLISECONDS )
	columnMap['s']     = Column.ColumnSetupTime( Column.ColumnTime.TimeUnit.SECONDS )
        columnMap['sm']    = Column.ColumnSetupTime( Column.ColumnTime.TimeUnit.MILLISECONDS )
	columnMap['t']     = Column.ColumnTearDownTime( Column.ColumnTime.TimeUnit.SECONDS )
        columnMap['tm']    = Column.ColumnTearDownTime( Column.ColumnTime.TimeUnit.MILLISECONDS )
	columnMap['b']     = Column.ColumnBandwidth( Column.ColumnBandwidth.ByteUnit.GIGA_BYTE )
        columnMap['bm']    = Column.ColumnBandwidth( Column.ColumnBandwidth.ByteUnit.MEGA_BYTE )
        columnMap['bt']    = Column.ColumnBandwidth( Column.ColumnBandwidth.ByteUnit.TERA_BYTE )
        columnMap['bp']    = Column.ColumnBandwidth( Column.ColumnBandwidth.ByteUnit.PETA_BYTE )
        columnMap['be']    = Column.ColumnBandwidth( Column.ColumnBandwidth.ByteUnit.EXA_BYTE )
	columnMap['f']     = Column.ColumnFlops( Column.ColumnFlops.FlopUnit.GIGA_FLOP )
        columnMap['fm']    = Column.ColumnFlops( Column.ColumnFlops.FlopUnit.MEGA_FLOP )
        columnMap['ft']    = Column.ColumnFlops( Column.ColumnFlops.FlopUnit.TERA_FLOP )
        columnMap['fp']    = Column.ColumnFlops( Column.ColumnFlops.FlopUnit.PETA_FLOP )
        columnMap['fe']    = Column.ColumnFlops( Column.ColumnFlops.FlopUnit.EXA_FLOP )
	columnMap['nt']    = Column.ColumnNumThreads( )
	columnMap['nir']   = Column.ColumnInnerRepetition( )
	columnMap['nor']   = Column.ColumnOuterRepetition( )
	columnMap['ns']    = Column.ColumnExecutedTime( Column.ColumnTime.TimeUnit.SECONDS )
        columnMap['nms']   = Column.ColumnExecutedTime( Column.ColumnTime.TimeUnit.MILLISECONDS )
	columnMap['lc']    = Column.ColumnConfLaunch( )
	return columnMap

def getBenchmarkOut( columns ):
	columnMap  = getColumnMap( )
	columnList = []

	for col in columns.split( ',' ):
		try:
			columnList.append( columnMap[col] )
		except KeyError:
			print colors.FAIL + "Unknown argument for column: '%s'. Type --longhelp to get the valid arguments of column and their meaning.%s" % ( col,colors.END )
			sys.exit( 1 )

	return BenchmarkOut.BenchmarkOut( columnList )

def printBeru( benchmarkMap,isetList,stream,filename ):

	print >> stream, '#'
	print >> stream, '# %s' % filename
	print >> stream, '#' 
	print >> stream, '# Created on %s' % datetime.datetime.fromtimestamp( time.mktime( datetime.datetime.now( ).timetuple( ) ) )
	print >> stream, '#         by %s' % getpass.getuser( )
	print >> stream, '#\n'

	print >> stream, '# Comment in, to add an description of your file, which will be printed out, when using it.'
	print >> stream, '# DESCRIPTION: |'
	print >> stream, '#     description of %s by %s\n' % ( filename,getpass.getuser( ) )

	print >> stream, '# Comment in, to add definitions, for not typing benchmark or inputset twice.'
        print >> stream, '# DEFINITIONS: |'
        print >> stream, '#     &alias Name\n'

	print >> stream, '# The minimal number of repetition of execute( ) within the benchmark.'
	print >> stream, 'REP_MIN: 1\n'
	print >> stream, '# The minimal time of repetition of execute( ) within the benchmark in seconds.'
	print >> stream, 'TIME_MIN: 0\n'
	print >> stream, '# The minimal number of repetition of the benchmark, to avoid aberrations.'
	print >> stream, 'REP_CTRL: 3\n'

	print >> stream, '# Registered Benchmarks.'
	print >> stream, '# Benchmarks needing Arguments, are marked the following:'
	print >> stream, '#     Benchmark (<args>)'
	print >> stream, 'BENCHMARKS:'
	for key,val in benchmarkMap['BENCHMARKS'].iteritems( ):
		if val == '<parametered>':
			key += ' (<args>)'
		if ':' in key or '-' in key or '\'' in key:
			key = '"%s"' % key
		print >> stream, '    - %s' % key
	print >> stream, ''

	print >> stream, '# Registered InputSets.'
	print >> stream, 'INPUT_SETS:'
	for iset in isetList:
		print >> stream, '    - %s' % iset
	print >> stream, ''

	print >> stream, '# Configuration of environment. Replace ID with the name of your configuration'
        print >> stream, '# launch and the empty spaces with the command you wish to have executed, right'
        print >> stream, '# before the execution of the benchmark.'
	print >> stream, 'LAUNCH_CONF:'
	print >> stream, '    default: " "\n'

	print >> stream, '# Put commandline arguments, here.'
	print >> stream, '# CMD: <command-line-arguments>\n'

	stream.close( )

def __min( list,key ):
	values = []
	for map in list:
		values.append( float( map[key] ) )
	return min( values )


def __aberration( list,key,threshold ):
	if len( list )==0:
		return True
	if len( list )==1:
		return False
	min     = __min( list,key )
	repeat  = True
	deleted = False

	while repeat:
		aberration = 0
		key_list   = []
		val_list   = []
		for map in list:
			key_list.append( ( float( map[key] )-min )**2 )
			val_list.append( float( map[key] ) )
		aberration = math.sqrt( 1.0/( len( key_list )-1 ) * sum( key_list ) )
		aberration = aberration / min

		if aberration > threshold:
			deleted   = True
			deletable = []
			maximum   = max( val_list )
			for map in list:
				if float( map[key] ) == maximum:
					deletable.append( map )
			for map in deletable:
				list.remove( map )
			repeat = ( len( list )-1 )
		else:
			repeat = False
	return deleted

# @description: returns, whether an aberration exists in given list.
def _aberration( list,quantity,threshold ):
	if len( list ) < quantity:
		return True

	abb       = False
	global MAP_NAME_LIST

	for key in MAP_NAME_LIST:
		if len( list ) < quantity:
			abb = True
			break
		if __aberration( list,key,threshold ):
			abb = True
	return abb


def _minimum_result( list ):

	if len( list ) < 2:
		if len( list ) == 0:
			return {}
		else:
			return list[0]

	map = list[0]

	minex = []
	for benchmap in list:
		minex.append( float( benchmap['MINEX'] ) )

	minexecutiontime = min( minex )

	for benchmap in list:
		if float( benchmap['MINEX'] ) == minexecutiontime:	
			map = benchmap

	return map

def __get_write_stream( filename ):
        while os.path.exists( filename[0] ):
                answered = False

	        while not answered:
	                a = raw_input( '\n' + filename[0] + " already exists. Overwrite? (y|n) " )
	                if a in ['y','n']:
	                        answered = True
	                else:
	                        print "Please type y or n!"

	        if a == 'n':
	                filename[0] = raw_input( "\nType new filename or q[uit] or b[reak]: " )
	                if filename[0] in ['q','quit']:
	                        sys.exit( 0 )
	                elif filename[0] in ['b','break']:
	                        return None
		else:
			break 
	try:
		stream = open( filename[0],'w' )
		return stream
	except IOError, e:
		print colors.FAIL,e,colors.END
		return None

def __print_results( benchOut,cmdOpt,stream=sys.stdout,csv=False ):
	if 'decimals' in cmdOpt:
                benchOut.setPrecision( int( cmdOpt['decimals'] ) )
        if 'language' in cmdOpt and not csv:
		try:
	        	benchOut.setLanguage( cmdOpt['language'] )
		except ValueError,e:
			print "%sERROR: %s%s" % ( colors.FAIL,e,colors.END )
			sys.exit( 1 )
	dir = os.path.realpath(os.getcwd( ))
        data_stream = open( dir + '/.log/result_map.map','r' )
        data_string = ''
        line        = data_stream.readline( )

        while line:
                data_string += line
                line         = data_stream.readline( )

        try:
                result_list  = pickle.loads( data_string )
        except EOFError:
                print colors.FAIL + "Benchmark results cannot be loaded." + colors.END
        data_stream.close( )

        minimum_list = []
        for results in result_list:
                minimum_map = _minimum_result( results )
                if minimum_map:
                        minimum_list.append( minimum_map )
        benchOut.print_results( minimum_list,stream,csv )

def main( ):
	c = colors( )
	# parse arguments
	cmdOpt = Options.get_opt( sys.argv[1:] )

        if cmdOpt['verbose']:
           print('cmdOpt = %s'%cmdOpt)

	# parsing file
	if cmdOpt['file']:
		fileMap = Configure.config( cmdOpt['file'] )
		fileOpt = {}

		# initialization of arguments, having no default value.
		output         = ''
		beruout        = ''
		console        = False
		run            = False
		errorreport    = ''
		listbenchmarks = False
		listinputsets  = False

		if fileMap['CMD']:
			fileOpt = Options.get_opt( fileMap['CMD'].split( ) )

			output         = fileOpt['output']
                        beruout        = fileOpt['beru-out']
			console        = fileOpt['console']
			run            = fileOpt['run']
			errorreport    = fileOpt['error-report']
			listbenchmarks = fileOpt['list-benchmarks']
			listinputsets  = fileOpt['list-input-sets']

		if 'columns' in cmdOpt:
			benchOut = getBenchmarkOut( cmdOpt['columns']  )
		elif 'columns' in fileOpt:
			benchOut = getBenchmarkOut( fileOpt['columns'] )
		else:
			benchOut = BenchmarkOut.BenchmarkOut( )		
		if cmdOpt['output']:
                        output = cmdOpt['output']
		if cmdOpt['beru-out']:
                        beruout = cmdOpt['beru-out']
		if 'global-repetition' in cmdOpt:
			globalrepetition = int( cmdOpt['global-repetition'] )
		elif 'global-repetition' in fileOpt:
			globalrepetition = int( fileOpt['global-repetition'] )
		else:
			globalrepetition = 1
		if 'avoid-aberration' in cmdOpt:
			avoidaberration = int( cmdOpt['avoid-aberration'] )
		elif 'avoid-aberration' in fileOpt:
                        avoidaberration = int( fileOpt['avoid-aberration'] )
		else:
			avoidaberration = 'unset'
		if 'maxrepetition' in cmdOpt:
			maxrepetition = int( cmdOpt['maxrepetition'] )
			if maxrepetition <= 0:
				print colors.FAIL + "The number of maximum iterations over benchmarks to avoid aberrations needs to be greater than zero.",
				print colors.END
				sys.exit( 1 )
		elif 'maxrepetition' in fileOpt:
			maxrepetition = int( fileOpt['maxrepetition'] )
                        if maxrepetition <= 0:
                                print colors.FAIL + "The number of maximum iterations over benchmarks to avoid aberrations needs to be greater than zero.",
				print colors.END
                                sys.exit( 1 )
		else:
			maxrepetition = -1
		if 'itemlist' in cmdOpt:
			itemlist = cmdOpt['itemlist']
		elif 'itemlist' in fileOpt:
			itemlist = fileOpt['itemlist']
		else:
			itemlist = ''
		if itemlist:
			global MAP_NAME_LIST
			ilists = itemlist.split( ';' )
			if len( ilists ) > 2:
				print colors.FAIL + itemlist + ": Too many lists (",len( ilists),") expected at maximum 2." + colors.END
				sys.exit( 1 )
			elif len( ilists ) < 1:
				print colors.FAIL + itemlist + ": Too little lists (",len( ilists),") expected at minimum 1." + colors.END
				sys.exit( 1 )
			for istring in ilists:
				if istring.startswith( '+' ):
					istring = istring.replace( '+','' )
					items = istring.split( ',' )
					for val in items:
						if val == 'S':
							MAP_NAME_LIST.append( 'SETUP' )
						elif val == 'T':
							MAP_NAME_LIST.append( 'TEARDOWN' )
						elif val not in ['S','T']:
							print colors.FAIL + "Invalid item: " + val + colors.END
							sys.exit( 1 )
				elif istring.startswith( '-' ):
					istring = istring.replace( '-','' )
					items = istring.split( ',' )
					for val in items:
						try:
							if val == 'E':
								MAP_NAME_LIST.remove( 'EXECUTION' )
							elif val == 'MIN':
								MAP_NAME_LIST.remove( 'MINEX' )
							elif val == 'MAX':
								MAP_NAME_LIST.remove( 'MAXEX' )
							else:
								print colors.FAIL + "Argument cannot be removed from list: " + val + colors.END
								sys.exit( 1 )
						except ValueError, e:
							print colors.FAIL,e,' (' + val + ')' + colors.END
							sys.exit( 1 )
				else:
					print colors.FAIL + istring + ': List must start either with + or -.' + colors.END
					sys.exit( 1 )
		if 'min-time' in cmdOpt:
			mintime = cmdOpt['min-time']
		elif 'min-time' in fileOpt:
			mintime = fileOpt['min-time']
		else:
			mintime = 0
		if 'threshold' in cmdOpt:
			threshold = float( cmdOpt['threshold'] )
		elif 'threshold' in fileOpt:
			threshold = float( fileOpt['threshold'] )
		else:
			threshold = 0.05
		if not ( 0<=threshold<=1 ):
			print colors.FAIL + "Unexpected value for threshold: %s. Expected threshold to be in range of [0,1].%s" % ( threshold,colors.END )
			sys.exit(1)
		if 'language' in cmdOpt:
			language = True
			try:
				benchOut.setLanguage( cmdOpt['language'] )
			except ValueError,e:
                	        print "%sERROR: %s%s" % ( colors.FAIL,e,colors.END )
                        	sys.exit( 1 )
	
		elif 'language' in fileOpt:
			language = True
			try:
				benchOut.setLanguage( fileOpt['language'] )
			except ValueError,e:
                        	print "%sERROR: %s%s" % ( colors.FAIL,e,colors.END )
	                        sys.exit( 1 )
		else:
			language = False
		if 'decimals' in cmdOpt:
			benchOut.setPrecision( int( cmdOpt['decimals'] ) )
		elif 'decimals' in fileOpt:
			benchOut.setPrecision( int( fileOpt['decimals'] ) )
		if 'path' in cmdOpt:
			path = cmdOpt['path']
		elif 'path' in fileOpt:
			path = fileOpt['path']
		else:
			path = './input/'
		if cmdOpt['error-report']:
			errorreport = cmdOpt['error-report']

                console        = cmdOpt['console'] or console
                run            = cmdOpt['run'] or run
                listbenchmarks = cmdOpt['list-benchmarks'] or listbenchmarks
                listinputsets  = cmdOpt['list-input-sets'] or listinputsets
	# no file is set.
	else:
		listbenchmarks   = cmdOpt['list-benchmarks']
                listinputsets    = cmdOpt['list-input-sets']
		beruout          = cmdOpt['beru-out']
		if 'columns' in cmdOpt:
			columns  = cmdOpt['columns']
		else:
			columns  = ''

		csv              = bool( cmdOpt['output'] ) 
		console          = cmdOpt['console']
		language         = 'language' in cmdOpt
		run              = False

		if not listbenchmarks and not listinputsets and not beruout and not columns and not 'decimals' in cmdOpt and not language and not csv:
			# Wrong parameters: Print help and exit.
			Options.get_opt( [] )

		if columns:
			benchOut    = getBenchmarkOut( columns )
			if csv:
				output = cmdOpt['output']
				print colors.OKPINK + "Creating csv-file %s.%s" % ( output,colors.END )
				output = [output]
	                        stream = __get_write_stream( output )
				output = output[0]
				if stream:
		                        if not language:
		                                benchOut.setLanguage( 'de' )
			                __print_results( benchOut,cmdOpt,stream,csv ) 
			                print >> stream,"Created on;%s" % datetime.datetime.fromtimestamp( time.mktime( datetime.datetime.now( ).timetuple( ) ) )
			                print >> stream,"By;Fraunhofer SCAI (c) 2011 Benchmark-Framework"
			                print >> stream,"Author;Robin Rehrmann"
			                stream.close( )
			                print colors.OKGREEN + "Created %s.%s" % ( output,colors.END )
		                if not console and not beruout:
		                        sys.exit( 0 )
		                elif not language:
		                        benchOut.setLanguage( 'en' )
			__print_results( benchOut,cmdOpt )
		elif 'decimals' in cmdOpt or language or csv:
			benchOut    = BenchmarkOut.BenchmarkOut( )
			if csv:
                                output = cmdOpt['output']
                                print colors.OKPINK + "Creating csv-file %s.%s" % ( output,colors.END )
				output = [output]
                                stream = __get_write_stream( output )
				output = output[0]
				if stream:
	                                if not language:
	                                        benchOut.setLanguage( 'de' )
	                                __print_results( benchOut,cmdOpt,stream,csv )
	                                print >> stream,"Created on;%s" % datetime.datetime.fromtimestamp( time.mktime( datetime.datetime.now( ).timetuple( ) ) )
	                                print >> stream,"By;Fraunhofer SCAI (c) 2011 Benchmark-Framework"
	                                print >> stream,"Author;Robin Rehrmann"
	                                stream.close( )
	                                print colors.OKGREEN + "Created %s.%s" % ( output,colors.END )
                                if not console and not beruout:
                                        sys.exit( 0 )
                                elif not language:
                                        benchOut.setLanguage( 'en' )
			if console or not csv:
				__print_results( benchOut,cmdOpt )
			

	#                                                             #
	# Print out list of benchmarks and/or inputsets, if required. #
	#                                                             #
	if listinputsets:
                isetMap = Run.getInputSets( cmdOpt['verbose'] )
                if isetMap['ERROR']:
                        print colors.FAIL + isetMap['MESSAGE'] + colors.END
                        sys.exit( 1 )
                isetList = isetMap['ISETS']
	if listbenchmarks:
		benchmarkMap = Run.getBenchmarks( cmdOpt['verbose'] )
		if benchmarkMap['ERROR']:
			print colors.FAIL + benchmarkMap['MESSAGE'] + colors.END
                        sys.exit( 1 )

	if listbenchmarks and listinputsets:
		table = texttable.Texttable( 0 )
		table.set_deco( texttable.Texttable.BORDER | texttable.Texttable.HEADER | texttable.Texttable.VLINES )
		table.header( ['Benchmark ID','Benchmark Name','   ','Registered InputSets'] )
		i = 0
		for key,value in benchmarkMap['BENCHMARKS'].iteritems( ):
			row = [ key,value,'' ]
			if i < len( isetList ):
				row.append( isetList[i] )
			else:
				row.append( '' )
			i += 1
			table.add_row( row )
		while i < len( isetList ):
			table.add_row( [ '','','',isetList[i] ] )
			i += 1
		print table.draw( )
		if not run or not beruout:
			sys.exit( 0 )
	elif listbenchmarks:
		table = texttable.Texttable( 0 )
                table.set_deco( texttable.Texttable.BORDER | texttable.Texttable.HEADER | texttable.Texttable.VLINES )
		table.header( ['Benchmark ID','Benchmark Name'] )
		for key,value in benchmarkMap['BENCHMARKS'].iteritems( ):
                        table.add_row( [ key,value, ] )
		print table.draw( )
		if not run or not beruout:
			sys.exit( 0 )
	elif listinputsets:
		table = texttable.Texttable( 0 )
                table.set_deco( texttable.Texttable.BORDER | texttable.Texttable.HEADER | texttable.Texttable.VLINES )
		table.header( ['Registered InputSets'] )
		for iset in isetList:
			table.add_row( [iset] )
		print table.draw( )
		if not run or not beruout:
			sys.exit( 0 )

	#                                                                   #
	# Print out all benchmarks and inputsets to beru-file, if required. #
	#                                                                   #

	if beruout:
		if not vars( ).has_key( 'isetList' ):
			isetMap = Run.getInputSets( cmdOpt['verbose'] )
	                if isetMap['ERROR']:
	                        print colors.FAIL + isetMap['MESSAGE'] + colors.END
	                        sys.exit( 1 )
	                isetList = isetMap['ISETS']
		if not vars( ).has_key( 'benchmarkMap' ):
			benchmarkMap = Run.getBenchmarks( 0 )
	                if benchmarkMap['ERROR']:
	                        print colors.FAIL + benchmarkMap['MESSAGE'] + colors.END
	                        sys.exit( 1 )

		beruout = [beruout]
		stream = __get_write_stream( beruout )
		beruout = beruout[0]

		if stream:
			print colors.OKPINK + "Creating file %s.%s" % ( beruout,colors.END )
			try:
				printBeru( benchmarkMap,isetList,stream,beruout )
			except:
				print colors.FAIL + 'Unexpected error: ', sys.exc_info()[0], colors.END
				sys.exit( 1 )
			print colors.OKGREEN + "Created %s.%s" % ( beruout,colors.END )

	if not run and beruout:
		sys.exit( 0 )

	#                #
	# Run Benchmarks #
	#                #

	if not vars( ).has_key( 'fileMap' ):
		sys.exit( 0 )

	if fileMap['DESCRIPTION']:
		print colors.MESSAGE
		print fileMap['DESCRIPTION']
		print colors.END

	tmp_dir = tempfile.mkdtemp( dir=os.getcwd( ) )

	actbench = 0
	numbench = len( fileMap['COMBINATIONS'] )
	all_results = []
	warning_list  = []
	warning_num   = 0
	keyInterrupt  = False

	# stop the time
	t_start = datetime.datetime.now( ).replace( microsecond=0 )

        maxrepetition_option = maxrepetition
	for set in fileMap['COMBINATIONS']:
		maxrepetition = maxrepetition_option
		actbench += 1
		print "# Running Benchmark %d of %d: %s with %s" % ( actbench,numbench, set['BENCH'], set['CONF'].id )
		if avoidaberration == 'unset':
			numrep = set['REP_CTRL']
		else:
			numrep = avoidaberration
		skipped = False
		this_result = []
		for actrep in range( numrep ):
			print "\n- Avoiding Discordant Values %d of %d" % ( actrep+1,numrep )

			try:
				resultMap = Run.benchmark( set['CONF'].value, 
                                             set['BENCH'] + ' "' + path + '" "' + tmp_dir + '"', cmdOpt['verbose'] )
			except Interrupt:
				print colors.WARNING + "Skipped" + colors.END
				skipped = True
				message = "Skipped Benchmark %d at iteration %d" % ( actbench,actrep+1 )
				warning_list.append( message )
				break
			except KeyboardInterrupt:
				keyInterrupt = True
                                message = "KeyboardInterrupt at %s" % datetime.datetime.fromtimestamp(time.mktime(datetime.datetime.now( ).timetuple( )))
                                warning_list.append( message )
				break

			if resultMap['ERROR']:
				print colors.FAIL + "ERROR: " + resultMap['MESSAGE'] + colors.END
				sys.exit( 1 )

			if resultMap['WARNING']:
				warning_num += 1
				warning_list.append( resultMap['MESSAGE'] )
				print colors.WARNING + "WARNING: " + resultMap['MESSAGE'] + colors.END
				break

			resultMap['CONF'] = set['CONF']

			this_result.append( resultMap )
			for map in this_result:
                                map['ADISC'] = actrep + 1

			maxrepetition = maxrepetition and maxrepetition - 1

		aberration = not keyInterrupt and not skipped and not resultMap['WARNING'] and _aberration( this_result,numrep,threshold )

		while aberration and maxrepetition:
			actrep += 1
			maxrepetition = maxrepetition and maxrepetition - 1
                        print "\n+ Avoiding Discordant Values %d of %d" % ( actrep+1,numrep )

                        try:
	                        resultMap = Run.benchmark( set['CONF'].value, 
                                            set['BENCH'] + ' "' + path + '" "' + tmp_dir + '"', cmdOpt['verbose'] )
                        except Interrupt:
                                print colors.WARNING + "Skipped" + colors.END
                                skipped = True
                                message = "Skipped Benchmark %d at iteration %d" % ( actbench,actrep+1 )
                                warning_list.append( message )
                                break
                        except KeyboardInterrupt:
				keyInterrupt = True
				message = "KeyboardInterrupt at %s" % datetime.datetime.fromtimestamp(time.mktime(datetime.datetime.now( ).timetuple( )))
				warning_list.append( message )
				break

                        if resultMap['ERROR']:
                                print colors.FAIL + "ERROR: " + resultMap['MESSAGE'] + colors.END
                                sys.exit( 1 )

                        if resultMap['WARNING']:
                                warning_num += 1
                                warning_list.append( resultMap['MESSAGE'] )
                                print colors.WARNING + "WARNING: " + resultMap['MESSAGE'] + colors.END
                                break

			resultMap['CONF'] = set['CONF']

                        this_result.append( resultMap )
			for map in this_result:
				map['ADISC'] = actrep + 1
			aberration = _aberration( this_result,numrep,threshold )

		all_results.append( this_result )
		
		print "\n# Finished Benchmark %d of %d\n" % ( actbench,numbench )

		if keyInterrupt:
			break

	t_stop = datetime.datetime.now( ).replace( microsecond=0 )

	dir = os.path.realpath(os.getcwd( ))
	if not os.path.exists(dir+"/.log"):
		executable = os.mkdir(dir+"/.log")
	stream = open( dir + '/.log/result_map.map','w' )
	print >> stream, pickle.dumps( all_results )

	minimum_list = []
	for list in all_results:
		minimum_map = _minimum_result( list )
		if minimum_map:
			minimum_list.append( minimum_map )

	for tmp_file in glob.glob( os.path.join( tmp_dir,'*' ) ):
		os.remove( tmp_file )
	os.rmdir( tmp_dir )

	#                              #
	# Print csv-File, if required, #
	#                              #

	if output:
		print colors.OKPINK + "Creating csv-file %s.%s" % ( output,colors.END )
		output = [output]
		stream = __get_write_stream( output )
		output = output[0]
		if stream:
			if not language:
				benchOut.setLanguage( 'de' )
			benchOut.print_results( minimum_list,stream,True )
			print >> stream,"Total execution time;%s" % ( t_stop-t_start )
			print >> stream,"Warnings;%s" % warning_num
			print >> stream,"Created on;%s" % datetime.datetime.fromtimestamp( time.mktime( datetime.datetime.now( ).timetuple( ) ) )
			print >> stream,"Location;%s" % socket.gethostname( )
			print >> stream,"By;Fraunhofer SCAI (c) 2011 Benchmark-Framework"
			print >> stream,"Author;Robin Rehrmann"
			stream.close( )
			print colors.OKGREEN + "Created %s.%s" % ( output,colors.END )
		if not console and not errorreport:
			sys.exit( 0 )
		elif not language:
			benchOut.setLanguage( 'en' )

	#                                                                                                   #
	# Print out results to console, unless file was printed to csv and only error-report was requested. #
	#                                                                                                   #

	if not output or not errorreport or console:
		print '\n\n'
		benchOut.print_results( minimum_list )
		print '\n\n'
		print '# Warnings: %s\n' % warning_num
		print '# Execution Time: %s\n' % ( t_stop-t_start )
		print '# Location: %s\n' % socket.gethostname( )
		print 'Finished at %s\n' % datetime.datetime.fromtimestamp( time.mktime( datetime.datetime.now( ).timetuple( ) ) )

	if keyInterrupt:
		print colors.FAIL + '\nInterrupted by user.\n' + colors.END

	#                                          #
	# Print out the error report, if required. #
	#                                          #

	if errorreport:
		if errorreport == 'stdout':
			stream = sys.stdout
		else:
			errorreport = [errorreport]
			stream = __get_write_stream( errorreport )

		if stream:
			print >> stream,'# ERROR REPORT %s #\n' % datetime.datetime.fromtimestamp( time.mktime( datetime.datetime.now( ).timetuple( ) ) )

			for error in warning_list:
				print >> stream,error
			print >> stream, '\n'
			stream.close( )

