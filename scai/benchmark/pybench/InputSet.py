import string

class InputSet:
	def __init__( self,name ):
		# Check, because the name may has arguments
                if string.find( name,'(' ) > 0 and string.find( name,')' ) > 0:
                        s = string.split( name,'(' )
                        self.name = s[0].strip( )
                        self.arguments = string.split( s[1],')' )[0].strip( )
                elif ( string.find( name,'(' ) > 0 )^( string.find( name,')' ) > 0 ):
                        message = "InputSet.__init__: missing '(' or ')' in %s" % name
                        raise SyntaxError( message )
                else:
                        self.name      = name.strip( )
                        self.arguments = ""

	def hasArgs( self ):
		return self.arguments != ""

	def __repr__( self ):
		if self.hasArgs( ):
			s = "%s( %s )" % ( self.name,self.arguments )
		else:
			s = self.name
		return s

