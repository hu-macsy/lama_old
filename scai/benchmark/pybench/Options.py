import getopt
import sys

OPT_HELP_MAP = {
	'f': '''
  -f,--file=FILE

  The beru-file Name for creating the BenchmarkRunner. If no file is set all,
  registered Benchmarks will be run in cobination with all registered InputSets
  and all sensfull numbers of threads.
             ''',
	'v': '''
  -v,--verbose    

  Print verbose information about issued benchmark commands.
             ''',
	'o': '''
  -o,--output=FILE

  The csv-file to write to.
             ''',
	'b': '''
  -b,--beru-out=FILE

  Writes the BenchmarkRunner to a beru-file from which it can be started again.
  Does not run benchmarks.
             ''',
	'r': '''
  -r,--global-repetition=VAL

  The number of repetitions for all Benchmarks. This value defines, how often
  the execution of function execute( ) will be repeated.
             ''',
	'a': '''
  -a, --avoid-aberration=VAL

  The number of repetitions for all Benchmarks. This value defines, how often a
  whole Benchmark will be minimally repeated, to find discordant values. This
  framework does execute a Benchmark at least VAL times, but more often, if dis-
  cordant values could be detected, but not identified.
             ''',
	'm': '''
  -m, --maxrepetition=VAL

  The maximum number of repetitions to avoid aberrations. This value may be set,
  if benchmarks seem to be repeated infinite. The repeating of the benchmark
  will then stop after VAL iterations, where VAL is a number greater than zero.
    However, every benchmark will at least be repeated REP_CTRL times, so if VAL
  was smaller than REP_CTRL, the appropriate benchmark will be stopped after
  REP_CTRL iterations, only.
	     ''',
	'i': '''
  -i, --itemlist=LIST

  Modification of the itemlist. The itemlist holds the values for which the
  aberration will be checked. These are the average execution time, the minimal
  execution time and the maximum execution time. To check for aberrations of other
  values like the setup time too this option changes the items to check.
  Therefore LIST is a list of items being removed from or
  added to the list of items checked. The format of LIST is the following:
    Items to be added are seperated by comma, where the list of items to be
  added starts with a plus. Items to be removed are also seperated by comma, but
  the list begins with a minus. The two lists may be seperated by semicolon.
  Whitespaces are not allowed. These are the abbreviations and their meanings,
  being allowed within LIST:
      E   EXECUTION the execution-time
      MIN MINIMAL EXECUTION the minimal execution-time
      MAX MAXIMAL EXECUTION the maximal execution-time
      S   SETUP the setup-time
      T   TEARDOWN the teardown-time
  E,MIN,MAX are checked by default an can therefor only removed from list, S and T
  can only be added to the list
  The first four options can only be removed from list, the last four can only
  be added to the list.
  Example:
       --itemlist '-MIN,MAX;+S,T'
	     ''',
	't': '''
  -t, --min-time=VAL

  The time running execute( ). execute( ) will be run
             ''',
	'd': '''
  -d, --decimals=VAL

  The number of decimals of the floating point numbers. Defaults to none.
             ''',
	'T': '''
  -T, --threshold=VAL

  The threshold, when testing for aberrations, with VAL in [0,1]. Defaults to
  0.05
             ''',
	'l': '''
  -l, --language=LAN

  The language for printing out. Used for formatting the numbers.
      'en' Uses dot as decimal marker.
      'de' Uses comma as decimal marker.
  Default is 'en'. If parameter -o is set, default is 'de'.
             ''',
	'c': '''
  -c, --columns=COL

  The columns of the output, where COL is a string containing the columns in the
  right order, sepereated by ','. The IDs of the columns are the following:
      n     The name of the benchmark.
      e     The runtime of the benchmark in seconds.
      em    The runtime of the benchmark in milliseconds.
      emax  The maximum runtime of the benchmark in seconds.
      emaxm The maximum runtime of the benchmark in milliseconds.
      emin  The minimum runtime of the benchmark in seconds.
      eminm The minimum runtime of the benchmark in milliseconds.
      s     The setup time of the benchmark in seconds.
      sm    The setup time of the benchmark in milliseconds.
      t     The tearDown time of the benchmark in seconds.
      tm    The tearDown time of the benchmark in milliseconds.
      b     The bandwidth of the benchmark in gigabyte.
      bm    The bandwidth of the benchmark in megabyte.
      bt    The bandwidth of the benchmark in terabyte (on 64-bit systems,
            only).
      bp    The bandwidth of the benchmark in petabyte (on 64-bit systems,
            only).
      be    The bandwidth of the benchmark in exabyte (on 64-bit systems, only).
      f     The number of floating point operations of the benchmark in giga-
            flop.
      fm    The number of floating point operations of the benchmark in mega-
            flop.
      ft    The number of floating point operations of the benchmark in tera-
            flop (on 64-bit systems, only).
      fp    The number of floating point operations of the benchmark in peta-
            flop (on 64-bit systems, only).
      fe    The number of floating point operations of the benchmark in exa-
            flop (on 64-bit systems, only).
      nt    The number of threads used for the benchmark.
      nir   The number of repetitions of the function execute( ) within the
            benchmark.
      nor   The number of repetitions of the benchmarksetting, due to avoidance
            of discordant values.
      nms   The time, the function execute( ) has been run in milliseconds.
      ns    The time, the function execute( ) has been run in seconds.
      lc    The launch configuration in which this benchmark run.
  Default is n,e,s,t,b,f
             ''',
	'p': '''
  -p, --path=PATH

  The path to input files. Defaults to ./input/.
             ''',
	'C': '''
  -C, --console

  Prints result always on console, even if it is printed to file. The default
  language for printing to console is 'en', even though the default for printing
  to file is 'de'.
             ''',
	'R': '''
  -R, --run

  Runs benchmarks, even if they are written out to beru-file.
             ''',
	'e': '''
  -e, --error-report=FILE

  Prints out the errors of the run after finishing. FILE can be either 'stdout'
  to print to console, or a filename, to print to file.
             ''',
	'B': '''
  -B, --list-benchmarks

  List all registered Benchmarks with IDs and Names and exit. If a Benchmark is
  parametered, it will be listed as
      BENCHMARK_ID: <parametered>
             ''',
	'I': '''
  -I, --list-input-sets

  List all registered InputSets and exit.
             '''
}

def help( ):
	print "Execution: BenchmarkRunner [OPTIONS]\n"
	print "-v, --verbose               Prints verbose information about issued commands."
	print "-f, --file=FILE             The beru-file Name for creating the BenchmarkRunner."
	print "-o, --output=FILE           The csv-file to write to."
	print "-b, --beru-out=FILE         Writes the BenchmarkRunner to a beru-file from which"
	print "                            it can be started again."
	print "-r, --global-repetition=VAL The number of repetitions for execute( ) of all"
	print "                            Benchmarks."
	print "-a, --avoid-aberration=VAL  The number of repetitions for all Benchmarks."
	print "-m, --maxrepetition=VAL     The maximum number of repetitions to avoid aber-"
        print "                            rations."
	print "-i, --itemlist=LIST         Modification of the itemlist."
	print "-t, --min-time=VAL          The time of running execute( )."
	print "-d, --decimals=VAL          The number of decimals of the floating point numbers."
	print "                            Defaults to none."
	print "-T, --threshold=VAL         The threshold, when testing for aberrations, with VAL"
	print "                            in [0,1]. Defaults to 0.05"
	print "-l, --language=LAN          The language for printing out."
	print "-c, --columns=COL           The columns of the output, where the Argument is a"
	print "                            string containing the columns in the right order,"
	print "                            seperated by ','."
	print "-p, --path=PATH             The path to input files. Defaults to ./input/."
	print "-C, --console               Prints result always on console."
	print "-R, --run                   Runs benchmarks."
	print "-e, --error-report=FILE     Prints out the errors of the run after finishing."
	print "-B, --list-benchmarks       List all registered Benchmarks and exit."
	print "-I, --list-input-sets       List all registered InputSets and exit."
	print "-h, --help                  Display this usage information and exit."
	print "-H, --opthelp=OPT           Print long help of OPT and exit."
	print "-L, --longhelp              Display long help and exit.\n"
#	print "BENCHMARK_IDS               The Ids of the benchmarks to be executed.\n\n"

def opthelp( key ):
	global OPT_HELP_MAP

	if key in( 'f','file' ):
		print OPT_HELP_MAP['f']
	elif key in( 'o','output' ):
		print OPT_HELP_MAP['o']
	elif key in( 'v','verbose' ):
		print OPT_HELP_MAP['v']
	elif key in( 'b','beru-out' ):
                print OPT_HELP_MAP['b']
        elif key in( 'r','global-repetition' ):
                print OPT_HELP_MAP['r']
        elif key in( 'a','avoid-aberration' ):
                print OPT_HELP_MAP['a']
	elif key in( 'm','maxrepetition' ):
		print OPT_HELP_MAP['m']
	elif key in( 'i','itemlist' ):
		print OPT_HELP_MAP['i']
        elif key in( 't','min-time' ):
                print OPT_HELP_MAP['t']
        elif key in( 'd','decimals' ):
                print OPT_HELP_MAP['d']
	elif key in( 'T','threshold' ):
                print OPT_HELP_MAP['T']
        elif key in( 'l','language' ):
                print OPT_HELP_MAP['l']
        elif key in( 'c','columns' ):
                print OPT_HELP_MAP['c']
        elif key in( 'p','path' ):
                print OPT_HELP_MAP['p']
        elif key in( 'C','console' ):
                print OPT_HELP_MAP['C']
        elif key in( 'R','run' ):
                print OPT_HELP_MAP['R']
        elif key in( 'e','error-report' ):
                print OPT_HELP_MAP['e']
        elif key in( 'B','list-benchmarks' ):
                print OPT_HELP_MAP['B']
        elif key in( 'I','list-input-sets' ):
                print OPT_HELP_MAP['I']
	else:
		print "No help for Option '%s' available." % key

def longHelp( ):
	print "Execution: BenchmarkRunner.py [OPTIONS]\n"
	print "-f, --file=FILE             The beru-file Name for creating the BenchmarkRunner."
	print "                            If no file is set all, registered Benchmarks will be"
	print "                            run in cobination with all registered InputSets and"
	print "                            all sensfull numbers of threads."
	print "-o, --output=FILE           The csv-file to write to."
	print "-b, --beru-out=FILE         Writes the BenchmarkRunner to a beru-file from which"
	print "                            it can be started again. Does not run benchmarks."
	print "-r, --global-repetition=VAL The number of repetitions for all Benchmarks. This"
	print "                            value defines, how often the execution of function"
	print "                            execute( ) will be repeated."
	print "-a, --avoid-aberration=VAL  The number of repetitions for all Benchmarks. This"
	print "                            value defines, how often a whole Benchmark will be"
	print "                            minimally repeated, to find discordant values. This"
	print "                            framework does execute a Benchmark at least VAL"
	print "                            times, but more often, if discordant values could be"
	print "                            detected, but not identified."
	print "-m, --maxrepetition=VAL     The maximum number of repetitions to avoid aber-"
	print "                            rations. This value may be set, if benchmarks seem to"
	print "                            be repeated infinite. The repeating of the benchmark"
	print "                            will then stop after VAL iterations, where VAL is a"
	print "                            number greater than zero."
	print "                            However, every benchmark will at least be repeated"
	print "                            REP_CTRL times, so if VAL was smaller than REP_CTRL,"
	print "                            the appropriate benchmark will be stopped after"
	print "                            REP_CTRL iterations, only."
	print "-i, --itemlist=LIST         Modification of the itemlist. The itemlist holds the"
	print "                            values for which the aberration will be checked. This"
	print "                            is The setup-time, the minimal execution-time, the"
	print "                            median execution-time and the teardown-time. For some"
	print "                            reasons it might be usefull, not to check one of"
	print "                            these for aberrations or check other values like the"
	print "                            maximum execution-time, too. Therefore LIST is a list"
	print "                            of items being removed from or added to the list of"
	print "                            items checked. The format of LIST is the following:"
	print "                                Items to be added are seperated by comma, where"
	print "                            the list of items to be added starts with a plus."
	print "                            Items to be removed are also seperated by comma, but"
	print "                            the list begins with a minus. The two lists may be"
	print "                            seperated by semicolon. Whitespaces are not allowed."
	print "                                These are the abbreviations and their meanings,"
	print "                            being allowed within LIST:"
	print "                                S   SETUP the setup-time"
	print "                                E   EXECUTION the execution-time"
	print "                                MIN MINIMAL EXECUTION the minimal execution-time"
	print "                                T   TEARDOWN the teardown-time"
	print "                                MAX MAXIMAL EXECUTION the maximal execution-time"
	print "                                F   FLOPS the number of flops"
	print "                                B   BANDWIDTH the number of processed bytes"
	print "                                TOT TOTAL EXECUTION the total execution-time"
	print "                            The first four options can only be removed from list,"
	print "                            the last four can only be added to the list."
	print "                            Example:"
	print "                                --itemlist '-S,T;+F,B'"
	print "-t, --min-time=VAL          The time running execute( ). execute( ) will be run"
	print "                            at least VAL seconds."
	print "-d, --decimals=VAL          The number of decimals of the floating point numbers."
	print "	                           Defaults to none."
        print "-T, --threshold=VAL         The threshold, when testing for aberrations, with VAL"
        print "                            in [0,1]. Defaults to 0.05"
	print "-l, --language=LAN          The language for printing out. Used for formatting"
	print "                            the numbers."
	print "                                en Uses dot as decimal marker."
	print "                                de Uses comma as decimal marker."
	print "                            Default is 'en'. If parameter -o is set, default is"
	print "                            'de'."
	print "-c, --columns=COL           The columns of the output, where COL is a string"
	print "                            containing the columns in the right order, sepereated"
	print "                            by ','. The IDs of the columns are the following:"
        print "                                n     The name of the benchmark"
	print "                                e     The runtime of the benchmark in seconds."
	print "                                em    The runtime of the benchmark in milli-"
	print "                                      seconds."
	print "                                emax  The maximum runtime of the benchmark in"
	print "                                      seconds."
	print "	                               emaxm The maximum runtime of the benchmark in"
	print "                                      milliseconds."
	print "	                               emin  The minimum runtime of the benchmark in"
	print "                                      seconds."
	print "                                eminm The minimum runtime of the benchmark in"
	print "                                      milliseconds."
	print "                                s     The setup time of the benchmark in seconds."
	print "                                sm    The setup time of the benchmark in milli-"
	print "                                      seconds."
	print "                                t     The tearDown time of the benchmark in"
	print "                                       seconds."
	print "                                tm    The tearDown time of the benchmark in"
	print "                                      milliseconds."
	print "                                b     The bandwidth of the benchmark in gigabyte"
	print "                                bm    The bandwidth of the benchmark in megabyte"
	print "                                bt    The bandwidth of the benchmark in terabyte"
	print "                                      (on 64-bit systems, only)"
	print "                                bp    The bandwidth of the benchmark in petabyte"
	print "                                      (on 64-bit systems, only)"
	print "                                be    The bandwidth of the benchmark in exabyte"
	print "                                      (on 64-bit systems, only)"
	print "                                f     The number of floating point operations of"
	print "                                      the benchmark in gigaflop."
	print "                                fm    The number of floating point operations of"
	print "                                      the benchmark in megaflop."
	print "                                ft    The number of floating point operations of"
	print "                                      the benchmark in teraflop. (on 64-bit"
	print "                                      systems, only)"
	print "                                fp    The number of floating point operations of"
	print "                                      the benchmark in petaflop. (on 64-bit"
	print "                                      systems, only)"
	print "                                fe    The number of floating point operations of"
	print "                                      the benchmark in exaflop. (on 64-bit"
	print "                                      systems, only)"
	print "                                nt    The number of threads used for the bench-"
	print "                                      mark."
	print "                                nir   The number of repetitions of the function"
	print "                                      execute( ) within the benchmark."
	print "                                nor   The number of repetitions of the benchmark-"
	print "                                      setting, due to avoidance of discordant"
	print "                                      values."
	print "                                nms   The time, the function execute( ) has been"
	print "                                      run in milliseconds."
	print "                                ns    The time, the function execute( ) has been"
        print "                                      run in seconds."
	print "                                lc    The launch configuration in which this"
	print "                                      benchmark run."
	print "	                           Default is: n,e,s,t,b,f"
	print "-p, --path=PATH             The path to input files. Defaults to ./input/."
	print "-C, --console               Prints result always on console, even if it is"
	print "                            printed to file. The default language for printing to"
	print "                            console is 'en', even though the default for printing"
	print "                            to file is 'de'."
	print "-R, --run                   Runs benchmarks, even if they are written out to"
	print "                            beru-file."
	print "-e, --error-report=FILE     Prints out the errors of the run after finishing."
	print "                            FILE can be either 'stdout' to print to console, or a"
	print "                            filename, to print to file."
        print "-B, --list-benchmarks       List all registered Benchmarks with IDs and Names and"
	print "                            exit. If a Benchmark is parametered, it will be"
	print "                            listed as"
	print "                                BENCHMARK_ID: <parametered>"
        print "-I, --list-input-sets       List all registered InputSets and exit."
	print "-h, --help                  Display usage information and exit."
        print "-H, --opthelp=OPT           Print long help of OPT and exit. OPT is supposed to"
	print "                            be either the short or the long option without the"
	print "                            appending en-dash."
	print "                              For example:"
	print "                                 printing out the help of option 'path', means"
	print "                                 OPT=p or OPT=path"
	print "-L, --longhelp              Display this help and exit.\n"
#	print "BENCHMARK_IDS   The Ids of the benchmarks to be executed. If this program "
#        print "    was called with an inputfile, the given benchmark will be added to the"
#        print "    list of benchmarks with all the InputSets and and thread numbers, defined"
#        print "    in the file. Otherwise, only the given benchmarks will be run with all"
#        print "    registered InputSets and all possible numbers of threads. So the input-"
#        print "    file may be used, to define a list of InputSets, to be run with different"
#        print "    benchmarks, given as parameters on the console.\n"
	print "Values set on the commandline do always overwrite values set in the inputfile.\n"

def get_opt( options ):
	if len( options ) == 0:
		help( )
		sys.exit( 1 )
	try:
		opts, args = getopt.getopt( options, "f:o:b:r:a:m:i:d:t:T:l:c:p:CRe:BIhvH:L", 
                                            ["file=", "output=", "beru-out=", "global-repetition=",
                                             "avoid-aberration=", "maxrepetition=", "itemlist=", 
                                             "min-time=", "decimals=", "threshold=", "language=", 
                                             "columns=", "path=", "console", "run", "error-report=",
                                             "list-benchmarks", "list-input-sets", "help", 
                                             "verbose", "opthelp=","longhelp"] )
	except getopt.GetoptError, err:
		print str( err )
		help( )
		sys.exit( 1 )

	if len( args ) > 0:
		print "Unknown Parameters: %s" % ( args )
		sys.exit( 1 )

	map = {}
	map['file']              = ""
	map['output']            = ""
	map['beru-out']          = ""
	map['console']           = bool( 0 )
	map['verbose']           = bool( 0 )
	map['run']               = bool( 0 )
	map['error-report']      = ""
	map['list-benchmarks']   = bool( 0 )
	map['list-input-sets']   = bool( 0 )
       
	for o,a in opts:
		if o in( '-f','--file' ):
			map['file'] = a
		elif o in( '-o','--output' ):
			map['output'] = a
		elif o in( '-b','--beru-out' ):
			map['beru-out'] = a
		elif o in( '-r','--global-repetition' ):
			map['global-repetition'] = a
		elif o in( '-a','--avoid-aborration' ):
			map['avoid-aberration'] = a
		elif o in( '-m','--maxrepetition' ):
			map['maxrepetition'] = a
		elif o in( '-i','--itemlist' ):
			map['itemlist'] = a
		elif o in( '-t','--min-time' ):
			map['min-time'] = a
		elif o in( '-d','--decimals' ):
			map['decimals'] = a
		elif o in( '-T','--threshold' ):
			map['threshold'] = a
		elif o in( '-l','--language' ):
			map['language'] = a
		elif o in( '-c','--columns' ):
			map['columns'] = a
		elif o in( '-p','--path' ):
			map['path'] = a
		elif o in( '-C','--console' ):
			map['console'] = bool( 1 )
		elif o in( '-v','--verbose' ):
			map['verbose'] = bool( 1 )
		elif o in( '-R','--run' ):
			map['run'] = bool( 1 )
		elif o in( '-e','--error-report' ):
			map['error-report'] = a
		elif o in( '-B','--list-benchmarks' ):
			map['list-benchmarks'] = bool( 1 )
		elif o in ( '-I','--list-input-sets'  ):
			map['list-input-sets'] = bool( 1 )
		elif o in( '-h','--help' ):
			help( )
			sys.exit( )
		elif o in( '-H','--opthelp' ):
			opthelp( a )
			sys.exit( )
		elif o in( '-L','--longhelp' ):
			longHelp( )
			sys.exit( )
		else:
			assert False, "unhandled option"

	return map
