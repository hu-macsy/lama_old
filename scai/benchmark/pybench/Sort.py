
def lessNumThreads( x,y ):
	if x['THREADS'] != 'Unthreadded' and y['THREADS'] == 'Unthreadded':
		return 1
	if x['THREADS'] == 'Unthreadded' and y['THREADS'] != 'Unthreadded':
		return -1
	if x['THREADS'] == 'Unthreadded' and y['THREADS'] == 'Unthreadded':
		return 0
	return cmp( int( x['THREADS'] ),int( y['THREADS'] ) )

def compareValueTypeSize( x,y ):
	return cmp( x['VTS'],y['VTS'] )

def compareInputSetIds( x,y ):
	return cmp( x['ISET'],y['ISET'] )

def compareGroupIds( x,y ):
	return cmp( x['GROUP'],y['GROUP'] )

