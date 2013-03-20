from Enum import enum
import locale
import os

if os.name in ( 'nt','ce' ):
    lan_map = {}
    lan_map['de'] = 'German'
    lan_map['en'] = 'English'
    lan_map['C']  = 'C'
else:
    lan_map = {}
    lan_map['de'] = 'de_DE'
    lan_map['en'] = 'en_GB'
    lan_map['C']  = 'C'

glob_lan = ''

def format_num( val,prec,lan ):

	global glob_lan

	if not isinstance( val,float ):
		if isinstance( val,int ) and prec != -1:
			return val + ' '*( prec+1 )
		return val
	if val == 0:
		return '0' + ' '*( prec+1 )

	if lan != glob_lan:
		locale.setlocale( locale.LC_NUMERIC,lan_map[lan] )
		glob_lan = lan
	try:
		if prec == -1:
	                return  locale.format( "%f",( val ),True )
		return locale.format( "%.*f",( prec,val ),True )
	except (ValueError, TypeError):
		return val


class AbstractColumn:
	def __init__( self,other=None ):
		if other == None:
			self.precision = -1
			self.language  = 'C'
			self.name      = ''
		elif isinstance( other,str ):
			self.precision = -1
        	        self.language  = 'C'
	                self.name      = other
		elif isinstance( other,AbstractColumn ):
			self.precision = other.precision
        	        self.language  = other.language
	                self.name      = other.name
		else:
			message = 'AbstractColumn.__init__( ): Do not know what to do with \'other\', since it is neither string nor AbstractColumn.'
			raise TypeError( message )

	def getColumnName( self ):
		return self.name

	def getData( self,map ):
		raise AttributeError( "getData( ) is not implemented." )

	def setLanguage( self,language ):
		if language in ['C','en','de']:
			self.language = language
		else:
			message = "AbstractColumn::setLanguage: Unsupported value for language: '" + language + "'"
			raise ValueError( message )

	def getLanguage( self ):
		return self.language

	def setPrecision( self,precision ):
		if precision < -1:
			raise ValueError( "AbstractColumn.setPrecision( ): precision must be greater than -1." )
		self.precision = precision

	def getPrecision( self ):
		return self.precision

	def copy( self ):
		return AbstractColumn( self )

	def formateToString( self,d ):
		raise Error( "not implemented, yet." )

	def getAdjustment( self ):
		return 'r'

	
def byteUnitToString( byteUnit ):
	if byteUnit == ColumnBandwidth.ByteUnit.EXA_BYTE:
		return '(Eb/s)'
	if byteUnit == ColumnBandwidth.ByteUnit.PETA_BYTE:
		return '(Pb/s)'
	if byteUnit == ColumnBandwidth.ByteUnit.TERA_BYTE:
		return '(Tb/s)'
        if byteUnit == ColumnBandwidth.ByteUnit.GIGA_BYTE:
		return '(Gb/s)'
        if byteUnit == ColumnBandwidth.ByteUnit.MEGA_BYTE:
		return '(Mb/s)'
	raise ValueError( "Column::byteUnitToString: unspecified value for Bytes. ("+`byteUnit`+")" )


class ColumnBandwidth( AbstractColumn ):
	ByteUnit = enum(
		    EXA_BYTE   = 1000000000000000000,
      		PETA_BYTE  = 1000000000000000,
       		TERA_BYTE  = 1000000000000,
		    GIGA_BYTE  = 1000000000,
       		MEGA_BYTE  = 1000000
	)

	def __init__( self,byteUnit ):
		if isinstance( byteUnit, int ) or isinstance( byteUnit, long ):
			AbstractColumn.__init__( self,'Bandwidth '+byteUnitToString( byteUnit ) )
			self.byteUnit = byteUnit
		elif isinstance( byteUnit, ColumnBandwidth ):
			AbstractColumn.__init__( self,byteUnit )
			self.byteUnit = byteUnit.byteUnit
		else:
			message = 'ColumnBandwidth.__init__( ): Do not know what to do with \'other = \', since it is neither int nor ColumnBandwidth.'
			raise TypeError( message )

	def getData( self,map ):
		return format_num( float( map['BANDWIDTH'] )/self.byteUnit,self.precision,self.language )

	def copy( self ):
		return ColumnBandwidth( self )

def timeUnitToString( timeUnit ):
	if timeUnit == ColumnTime.TimeUnit.SECONDS:
		return '(s)'
        if timeUnit == ColumnTime.TimeUnit.MILLISECONDS:
		return '(ms)'
	raise ValueError( "Column::timeUnitToString: unspecified value for Time. ("+`timeUnit`+")" )

class ColumnTime( AbstractColumn ):
	TimeUnit = enum(
		SECONDS      = 1,
		MILLISECONDS = 1000
	)

	def __init__( self,timeUnit,name=None ):
		if isinstance( timeUnit,int ):
			if name == None:
				AbstractColumn.__init__( self )
			else:
				AbstractColumn.__init__( self,name )
			self.timeUnit = timeUnit
		elif isinstance( timeUnit,ColumnTime ):
			AbstractColumn.__init__( self,timeUnit )
	                self.timeUnit = timeUnit.timeUnit
		else:
			message = 'ColumnTime.__init__( ): Do not know what to do with \'other\', since it is neither int nor ColumnTime.'
                        raise TypeError( message )

	def getData( self,map ):
		return format_num( float( self.getTime( map ) )*self.timeUnit,self.precision,self.language )

	def getTime( self ):
		raise AttributeError( "getTime( ) is not implemented." )

class ColumnExecutionTime( ColumnTime ):
	def __init__( self,timeUnit ):
		if isinstance( timeUnit,int ):
			ColumnTime.__init__( self,timeUnit,'Runtime '+timeUnitToString( timeUnit ) )
		else:
			ColumnTime.__init__( self,timeUnit )

	def copy( self ):
		return ColumnExecutionTime( self )

	def getTime( self,map ):
		return map['EXECUTION']

class ColumnMaxExecutionTime( ColumnTime ):
	def __init__( self,timeUnit ):
                if isinstance( timeUnit,int ):
			ColumnTime.__init__( self,timeUnit,'Runtime(Max) '+timeUnitToString( timeUnit ) )
                else:
                        ColumnTime.__init__( self,timeUnit )

	def copy( self ):
		return ColumnMaxExecutionTime( self )

	def getTime( self,map ):
		return map['MAXEX']

class ColumnMinExecutionTime( ColumnTime ):
        def __init__( self,timeUnit ):
                if isinstance( timeUnit,int ):
	                ColumnTime.__init__( self,timeUnit,'Runtime(Min) '+timeUnitToString( timeUnit ) )
                else:
                        ColumnTime.__init__( self,timeUnit )

        def copy( self ):
                return ColumnMinExecutionTime( self )

        def getTime( self,map ):
                return map['MINEX']

class ColumnSetupTime( ColumnTime ):
	def __init__( self,timeUnit ):
                if isinstance( timeUnit,int ):
			ColumnTime.__init__( self,timeUnit,'Setup '+timeUnitToString( timeUnit ) )
                else:
                        ColumnTime.__init__( self,timeUnit )

	def copy( self ):
		return ColumnSetupTime( self )

	def getTime( self,map ):
		return map['SETUP']

class ColumnTearDownTime( ColumnTime ):
        def __init__( self,timeUnit ):
                if isinstance( timeUnit,int ):
	                ColumnTime.__init__( self,timeUnit,'Tear Down '+timeUnitToString( timeUnit ) )
                else:
                        ColumnTime.__init__( self,timeUnit )

        def copy( self ):
                return ColumnTearDownTime( self )

        def getTime( self,map ):
                return map['TEARDOWN']

class ColumnExecutedTime( ColumnTime ):
	def __init__( self,timeUnit ):
                if isinstance( timeUnit,int ):
			ColumnTime.__init__( self,timeUnit,'Running execute( ) '+timeUnitToString( timeUnit ) )
                else:
                        ColumnTime.__init__( self,timeUnit )

	def copy( self ):
		return ColumnExecutedTime( self )

	def getTime( self,map ):
		return map['EXECUTED']

def flopUnitToString( flopUnit ):
	if flopUnit == ColumnFlops.FlopUnit.EXA_FLOP:
		return '(Eflop/s)'
        if flopUnit == ColumnFlops.FlopUnit.PETA_FLOP:
                return '(Pflop/s)'
        if flopUnit == ColumnFlops.FlopUnit.TERA_FLOP:
                return '(Tflop/s)'
        if flopUnit == ColumnFlops.FlopUnit.GIGA_FLOP:
                return '(Gflop/s)'
        if flopUnit == ColumnFlops.FlopUnit.MEGA_FLOP:
                return '(Mflop/s)'
	raise ValueError( "Column::flopUnitToString( ): unspecified value for Flops. ("+`flopUnit`+")" )

class ColumnFlops( AbstractColumn ):
	FlopUnit = enum(
		EXA_FLOP   = 1000000000000000000,
        	PETA_FLOP  = 1000000000000000,
        	TERA_FLOP  = 1000000000000,
		GIGA_FLOP  = 1000000000,
	        MEGA_FLOP  = 1000000
	)

	def __init__( self,flopUnit ):
		if isinstance( flopUnit,int ) or isinstance( flopUnit, long ):
			AbstractColumn.__init__( self,'Flops '+flopUnitToString( flopUnit ) )
			self.flopUnit = flopUnit
		elif isinstance( flopUnit,ColumnFlops ):
			AbstractColumn.__init__( self,flopUnit )
	                self.flopUnit = flopUnit.flopUnit
		else:
			message = 'ColumnFlops.__init__( ): Do not know what to do with \'other\', since it is neither int nor ColumnFlops.'
                        raise TypeError( message )

	def copy( self ):
		return ColumnFlops( self )

	def getData( self,map ):
		return format_num( float( map['FLOPS'] )/self.flopUnit,self.precision,self.language )

class ColumnIdName( AbstractColumn ):
	def __init__( self,other=None ):
		if other==None:
			AbstractColumn.__init__( self,'Name' )
		else:
			AbstractColumn.__init__( self,other )

	def copy( self ):
		return ColumnIdName( self )

	def getData( self,map ):
		return map['NAME']

	def getAdjustment( self ):
		return 'l'

class ColumnNumThreads( AbstractColumn ):
	def __init__( self,other=None ):
                if other==None:
			AbstractColumn.__init__( self,'Num Threads' )
		else:
	                AbstractColumn.__init__( self,other )

	def copy( self ):
		return ColumnNumThreads( self )

	def getData( self,map ):
		if map['THREADS'] == 'Unthreadded':
			return 0
		return int( map['THREADS'] )

class ColumnInnerRepetition( AbstractColumn ):
	def __init__( self,other=None ):
		if other==None:
			AbstractColumn.__init__( self,'Repeated execute( )' )
		else:
			AbstractColumn.__init__( self,other )

	def copy( self ):
		return ColumnInnerRepetition( self )

	def getData( self,map ):
		return int( map['REPETITIONS'] )

class ColumnOuterRepetition( AbstractColumn ):
	def __init__( self,other=None ):
		if other==None:
                        AbstractColumn.__init__( self,'Repeated Benchmark' )
                else:
                        AbstractColumn.__init__( self,other )

	def copy( self ):
		return ColumnOuterRepetition( self )

	def getData( self,map ):
		return int( map['ADISC'] )

class ColumnConfLaunch( AbstractColumn ):
	def __init__( self,other=None ):
		if other==None:
			AbstractColumn.__init__( self,'Launch Configuration' )
		else:
			AbstractColumn.__init__( self,other )

	def copy( self ):
		return ColumnConfLaunch( self )

	def getData( self,map ):
		return map['CONF'].id

