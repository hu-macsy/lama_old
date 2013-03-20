import string

class Benchmark:
	"Benchmark( name,inputSetList,configurationList,time_min,rep_min,rep_ctrl ): Benchmark holding the list of its inputSet, its configuration options, its inner and outer repitition and minimum runtime."
	def __init__( self,name,inputSetList,configurationList,time_min,rep_min,rep_ctrl ):
                self.name = name
		self.inputSetList      = list( inputSetList )
		self.configurationList = list( configurationList )
		self.time_min = time_min
		self.rep_min  = rep_min
		self.rep_ctrl = rep_ctrl

	def setInputSetList( self,inputSetList ):
		"setInputSetList( inputSetList ): Deletes the old inputsetlist of this benchmark and sets the given one as the new one."
		del self.inputSetList
		self.inputSetList = list( inputSetList )

	def setConfigurationList( self,configurationList ):
		"setConfigurationList( configurationList ): Deletes the old configurationlist of this benchmark and sets the given one as the new one."
		del self.configurationList
		self.configurationList = list( configurationList )

	def setTimeMin( self,time_min ):
		self.time_min = time_min

	def setRepMin( self, rep_min ):
		self.rep_min = rep_min

	def setRepCtrl( self,rep_ctrl ):
		self.rep_ctrl = rep_ctrl

	def __eq__( a,b ):
		return a.name==b

	def __repr__( self ):
		s = "%s - %s - %s, " % ( self.name,self.inputSetList,self.configurationList )
		s = s+"%ds, %d in, %d out\n" % ( self.time_min,self.rep_min,self.rep_ctrl )
		return s

	def combinations( self,combiList ):
		"combinations( combiList ): Appends the given list with maps. These have three values:\n\tREP_CTRL: The outer repitition of this benchmark, ment for the python framework \n\tCONF:     The configuration launch of this benchmark, ment to be executed, before running this benchmark \n\tBENCH:    The arguments of this benchmark for the c++-program, in the following order: \n\t\t\tbenchmarkname inputsetname time repitition [arguments]"
		map = { "REP_CTRL": self.rep_ctrl }
		for conf in self.configurationList:
			map['CONF']  = conf
			for iset in self.inputSetList:
				s = "\"%s\" \"%s\" " % ( self.name,iset )
				if  self.time_min%2 == 0 or self.time_min%2 == 1:
					s += "%d %d" % ( self.time_min,self.rep_min )
				else:
					s += "%.3f %d" % ( self.time_min,self.rep_min )
				map['BENCH'] = s
				combiList.append( dict( map ) )
