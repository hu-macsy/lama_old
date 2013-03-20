from Benchmark import Benchmark
from InputSet import InputSet
from LaunchConfiguration import LaunchConfiguration

import yaml
import string

def addSkope( map,gBenchList,gBenchNameList,gISetList,gConfList,timeMin,repMin,repCtrl ):

	# local lists
	locISetList = []
	locConfList = []

	keywordList = map.keys( )

	if 'INPUT_SETS' in map:
		keywordList.remove( 'INPUT_SETS' )
		for iset in map['INPUT_SETS']:
			gISetList.append( iset )
			locISetList.append( iset )

	if 'LAUNCH_CONF' in map:
		keywordList.remove( 'LAUNCH_CONF' )
		iter = map['LAUNCH_CONF'].iteritems( )
		for pair in iter:
			conf = LaunchConfiguration( pair[0],pair[1] )
                        gConfList.append( conf )
			locConfList.append( conf )

	if 'TIME_MIN' in map:
		keywordList.remove( 'TIME_MIN' )
		timeMin = map['TIME_MIN']

	if 'REP_MIN' in map:
		keywordList.remove( 'REP_MIN' )
		repMin = map['REP_MIN']

	if 'REP_CTRL' in map:
                keywordList.remove( 'REP_CTRL' )
		repCtrl = map['REP_CTRL']

	if locISetList or locConfList:
		for glob_bench in gBenchNameList:
			if ( 'BENCHMARKS' in map and glob_bench not in map['BENCHMARKS'] ) or ( 'BENCHMARKS' not in map ):
				if not locISetList:
					bench = Benchmark( glob_bench,gISetList,locConfList,timeMin,repMin,repCtrl )
				elif not locConfList:
					bench = Benchmark( glob_bench,locISetList,gConfList,timeMin,repMin,repCtrl )
				else:
					bench = Benchmark( glob_bench,locISetList,locConfList,timeMin,repMin,repCtrl )
				gBenchList.append( bench )

	if 'BENCHMARKS' in map:
                keywordList.remove( 'BENCHMARKS' )
		for benchName in map['BENCHMARKS']:
			if benchName in gBenchNameList:
				bench = Benchmark( benchName,locISetList,locConfList,timeMin,repMin,repCtrl )
				if locISetList == []:
					bench.setInputSetList( gISetList )
				if locConfList == []:
					bench.setConfigurationList( gConfList )
				gBenchList.append( bench )
			else:
				gBenchNameList.append( benchName )
				bench = Benchmark( benchName,gISetList,gConfList,timeMin,repMin,repCtrl )
				gBenchList.append( bench )

	for keyword in keywordList:
		addSkope( map[keyword],gBenchList,list( gBenchNameList ),list( gISetList ),list( gConfList ),timeMin,repMin,repCtrl )

def config( filename ):
	stream = file( filename,'r' )
	map = yaml.load( stream )
	stream.close( )

	benchmarkList = []

	if 'CMD' in map:
		cmd = map['CMD']
		del map['CMD']
	else:
		cmd = ''

	if 'DESCRIPTION' in map:
		description = map['DESCRIPTION']
		del map['DESCRIPTION']
	else:
		description = ''

	if 'DEFINITIONS' in map:
		del map['DEFINITIONS']

	addSkope( map,benchmarkList,[],[],[],0,1,3 )

	combiList = []
	for benchmark in benchmarkList:
		benchmark.combinations( combiList )

	fileMap                 = {}

	fileMap['CMD']          = cmd
	fileMap['DESCRIPTION']  = description
	fileMap['COMBINATIONS'] = combiList

	return fileMap
