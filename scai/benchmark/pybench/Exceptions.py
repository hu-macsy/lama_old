
class Interrupt( KeyboardInterrupt ):
	def __init__( self,value=None ):
		if value is None:
			self.value = ''
		else:
			self.value = value

	def __str__( self ):
		return repr( self.value )

