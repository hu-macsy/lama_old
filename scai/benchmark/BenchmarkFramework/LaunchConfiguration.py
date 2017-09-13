

class LaunchConfiguration:
	def __init__( self,id,value ):
		self.id = id
		self.value = value

	def __repr__( self ):
		s = "%s: %s" % ( self.id,self.value )
		return s

