from Column import *
import texttable
import sys
import Sort

def getTerminalSize( ):
        def ioctl_GWINSZ(fd):
                try:
                        import fcntl, termios, struct, os
                        cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,'1234'))
                except:
                        return None
                return cr

        cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
        if not cr:
                try:
                        fd = os.open(os.ctermid(), os.O_RDONLY)
                        cr = ioctl_GWINSZ(fd)
                        os.close(fd)
                except:
                        pass
        if not cr:
                try:
                        cr = (env['LINES'], env['COLUMNS'])
                except:
                        cr = (25, 80)
        return int(cr[1]), int(cr[0])

class BenchmarkOut:
    def __init__( self,columnList=None ):
        self.columnList = []
        if columnList==None:
            self.columnList.append( ColumnIdName( ) )
            self.columnList.append( ColumnExecutionTime( ColumnTime.TimeUnit.SECONDS ) )
            self.columnList.append( ColumnSetupTime( ColumnTime.TimeUnit.SECONDS ) )
            self.columnList.append( ColumnTearDownTime( ColumnTime.TimeUnit.SECONDS ) )
            self.columnList.append( ColumnBandwidth( ColumnBandwidth.ByteUnit.GIGA_BYTE ) )
            self.columnList.append( ColumnFlops( ColumnFlops.FlopUnit.GIGA_FLOP ) )
        else:
            for col in columnList:
                self.columnList.append( col.copy( ) )

	def setColumns( self,columnList ):
		if len( columnList ) == 0:
			raise ValueError( 'BenchmarkOut.setColumns( ): Cannot print zero Columns.' )

		del self.columnList
		self.columnList = []

		for col in columnList:
			self.columnList.append( col.copy( ) )

    def getColumns( self ):
        return list( self.columnList )

    def setLanguage( self,language ):
        for col in self.columnList:
            col.setLanguage( language )

    def setPrecision( self,precision ):
        for col in self.columnList:
            col.setPrecision( precision )

    def print_results( self,benchmarkResults,stream=sys.stdout,csv=bool(0) ):
        """print_results( self,benchmarkResults,stream=sys.stdout,csv=False )\n\tPrints result as table to stream. The given parameters are:\n\t\tbenchmarkResults: A list of maps, holding the result of each benchmark.\n\n\nstream         : The stream to write to. Defaults to stdout.\n\t\tcsv            : Flag, whether to write to csv. Defaults to False.\n"""
        # sort list
        benchmarkResults.sort( Sort.lessNumThreads )
        benchmarkResults.sort( Sort.compareValueTypeSize )
        benchmarkResults.sort( Sort.compareInputSetIds )
        benchmarkResults.sort( Sort.compareGroupIds )

        if csv:
			table = texttable.Texttable( 0 ) # initialize table with infinite width
			table.set_deco( table.VLINES )
			table.set_chars( ['',';','',''] )
        else:
            width = os.popen( 'stty size','r' ).read( ).split( )[1]
            table = texttable.Texttable( width ) # initialize table with width of terminal

        col_width   = []
        adjustments = []
        vertical    = []
        header      = []

		# Compute length of columns. The final length of the columns is the widst
		# length of a column for a value.

        if csv:
            adjustments.append( 'l' )
            vertical.append( 'c' )
            header.append( "Group" )
            col_width.append( 5 )
            adjustments.append( 'l' )
            vertical.append( 'c' )
            header.append( "InputSet" )
            col_width.append( 5 )


        for col in self.columnList:
            adjustments.append( col.getAdjustment( ) )
            vertical.append( 'c' )
            header.append( col.getColumnName( ) )
            col_width.append( len( col.getColumnName( ) ) )


        for map in benchmarkResults:
            for i in range( len( self.columnList ) ):
                l = len( str( self.columnList[i].getData( map ) ) )
                if csv:
                    if l > col_width[i+2]:
                        col_width[i+2] = l
                else:
                    if l > col_width[i]:
                        col_width[i] = l
            if csv:
                lg = len( map['GROUP'].strip() )
                li = len( map['ISET'].strip() )
                if lg > col_width[0]:
                    col_width[0] = lg
                if li > col_width[1]:
                    col_width[1] = li

		# Sum up to length of table.
        table_length = ( len( self.columnList )-1 )*3
        for wid in col_width:
            table_length += wid

        table.set_cols_width( col_width )  # set the width of each column.
        table.set_cols_align( adjustments ) # set adjustment of each column ('l','c','r')
        table.set_cols_valign( vertical )  # set vertical adjustment of column ('t','c','b')

		# var to check, whether we enter a new group.
        groupId  = ''
        inputSet = ''

        # for csv just write the header once
        if csv:
            table.reset( )
            table.set_cols_width( col_width )
            table.set_cols_align( adjustments )
            table.header( header )

        # for each map of the list of maps, do loop over the columns and print the
        # results.
        for map in benchmarkResults:
			# Have we entered this group, yet?
            if map['GROUP'] != groupId:
				# if a table already exists, we print it.
                if table.draw( ) != None:
                    print >> stream,table.draw( )
				# a little table, just for the group
                table.reset( )
                if not csv:
                    table.set_cols_width( [table_length] )
                    table.set_cols_align( ['c'] )
                    table.set_chars( ['=','|','+','='] )
                    table.add_row( [map['GROUP']] )
                    print >> stream,table.draw( )
                    table.set_chars( ['-','|','+','='] )
				# reset table of the group and prepare for data.
				# and reset inputset.
                groupId = map['GROUP']
                inputSet = ''
                table.reset( )
                # no need to create new table, here, because a new inpuset will follow.

            # Have we entered this inputset within the group, yet?
            if map['ISET'] != inputSet:
                # if table already exists, we print it.
                if table.draw( ) != None:
                    print >> stream,table.draw( )

                # little table, just for the input set.
                table.reset( )
                if not csv:
                    table.set_cols_width( [table_length] )
                    table.set_cols_align( ['c'] )
                    table.set_cols_valign( ['c'] )
                    table.add_row( [map['ISET']] )
                    print >> stream,table.draw( )

                inputSet = map['ISET']

                table.reset( )
                table.set_cols_width( col_width )
                table.set_cols_align( adjustments )
                if not csv:
                    table.header( header )
			# var, holding a row
            results = []
            # loop over columns and receive data
            if csv:
                results.append( map['GROUP'] )
                results.append( map['ISET'] )
            for col in self.columnList:
                results.append( col.getData( map ) )
            table.add_row( results ) 

		# finally print the last bit of table, if any.
        if table.draw( ) != None:
            print >> stream,table.draw( )
		
