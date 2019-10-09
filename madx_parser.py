import os
import sys
import re
import math

#===============================================================

class _possibleElementType:
	"""
		Class. Specifies all possible element types
		"""
	def __init__(self):
		"""
			Constructor. Creates list of element types. 
			"""
		self.__names_type = []
		self.__names_type.append("drift")
		self.__names_type.append("aperture")
		self.__names_type.append("sbend")
		self.__names_type.append("rbend")
		self.__names_type.append("quad")
		self.__names_type.append("quadrupole")
		self.__names_type.append("sextupole")
		self.__names_type.append("octupole")
		self.__names_type.append("multipole")
		self.__names_type.append("solenoid")
		self.__names_type.append("kicker")
		self.__names_type.append("hkicker")
		self.__names_type.append("vkicker")
		self.__names_type.append("hkick")
		self.__names_type.append("vkick")
		self.__names_type.append("rfcavity")
		self.__names_type.append("rcollimator")
		self.__names_type.append("marker")
		self.__names_type.append("monitor")
	
	def __del__(self):
		"""
			Method. Deletes element.
			"""
		del self.__names_type
	
	def checkType(self, name_in):
		"""
			Method. Confirms validity of element type.
			"""
		name = name_in.lower()
		if self.__names_type.count(name) == 0:
			print "Error creating lattice element:"
			print "There is no element with type: ", name
			print "Stop."
			sys.exit (0)
		return name_in

#===============================================================

class MADX_LattElement:
	"""
		Class. Represents an arbitrary element in the lattice
		"""
	_typeChecker = _possibleElementType()
	
	def __init__(self, name, Typename):
		"""
			Constructor. Creates element with name, type,
			and parameter dictionary.
			"""
		self.__name = name
		self.__type = self._typeChecker.checkType(Typename)
		self.__par = {}
	
	def __del__(self):
		"""
			Method. Deletes parameters.
			"""
		del self.__par
	
	def getName(self):
		"""
			Method. Returns name of element
			"""
		return self.__name
	
	def setType(self, tp):
		"""
			Method. Sets the type of element without checking.
			"""
		self.__type = tp
	
	def getType(self):
		"""
			Method. Returns type of element
			"""
		return self.__type
	
	def addParameter(self, nameOfPar, parVal):
		"""
			Method. Adds parameter and value to element.
			"""
		self.__par[nameOfPar] = parVal
	
	def hasParameter(self, nameOfPar):
	
		if self.__par.has_key(nameOfPar) == 0:
			return 0

		else:
			return 1
	
	def getParameter(self, nameOfPar):
		"""
			Method. Returns name of parameter.
			"""
		if self.__par.has_key(nameOfPar) == 0:
			print "class MADX_LattElement, method getParameter"
			print "The name of Element = ", self.__name
			print "The type of Element = ", self.__type
			print "The Element's key-val = ", self.__par
			print "This Element does not have Parameter = ", nameOfPar
			print "Stop."
			sys.exit (0)
		return self.__par[nameOfPar]
	
	def getParameters(self):
		"""
			Method. Returns parameter dictionary.
			"""
		return self.__par
	
	def getElements(self):
		"""
			Method. Returns list of elements (only one here)
			"""
		elements = []
		elements.append(self)
		return elements

#====================================================================
class MADX_sequence:
	"""
	Class. Represents an arbitrary line in the lattice.
	"""
	def __init__(self, name):
		"""
		Constructor. Creates sequence with list of elements.
		"""
		self.__name = name
		self.__items = []
		
	def __del__(self):
		"""
		Method. Deletes instance with list of lines and/or elements.
		"""
		del self.__name
		del self.__items

	def getName(self):
		"""
		Method. Returns name of the line.
		"""
		return self.__name

	def getType(self):
		"""
		Method. Returns type: "sequence".
		"""
		return "SEQUENCE"

	def addItem(self, item):
		"""
		Method. Adds an element to this line.
		"""
		self.__items.append(item)

	def getItems(self):
		"""
		Method. Returns list with elements.
		"""
		return self.__items

#====================================================================

class MADX_Parser:
	""" MAD parser """
	
	def __init__(self):
		""" Create instance of the MAD_Parser class """
		self.__madxLines =  []
		self._sequences = []
		self.__accValues = []
		self._accElemDict = {}
		self._sequencename = ""
		self._sequencelength = ""
		self._sequencelist = []
		#the lines to ignore will start with these words
		self.__ingnoreWords = ["title","beam", "none"]


#---------------------------------
	def initialize(self,mad_file_name):
		fl = open(os.path.join(self.__madFilePath, mad_file_name))
		str_local = ""
		for str in fl.readlines():
			#check if the line is a comment
			if str.find("!") == 0:
				str_local = ""
				continue
			#STOP parsing a MAD file if there is a start of mad commands
			if(str.strip().lower().find("use") == 0):
				break
			str_local = "".join([str_local,str.strip()])
			#this part deal with a comment "!" at the end of line
			str0 = str_local
			if str0.rfind("!") > 0:
				str_local = ""
				for i in xrange(str0.rfind("!")):
					str_local = "".join([str_local,str0[i]])
			str_local.strip()
			#this part deal with a comment ";" at the end of line
			str0 = str_local
			if str0.rfind(";") > 0:
				str_local = ""
				for i in xrange(str0.rfind(";")):
					str_local = "".join([str_local,str0[i]])
			str_local.strip()
			#this part deal with a continue sign at the end of line
			if str_local.endswith("&"):
				str0 = str_local
				if str0.rfind("&") > 0:
					str_local = ""
					for i in xrange(str0.rfind("&")):
						str_local = "".join([str_local,str0[i]])
				str_local.strip()
				continue
			else:
				#check the empty line
				if str_local == "":
					str_local = ""
					continue
				#now we have the line to parse and the city to burn (wake up, samurai)
				self._init_string(str_local)
				str_local = ""
		fl.close()
#---------------------------------

	def parse(self,MADXfileName):
		
		self.__init__()
		
		#1-st stage read MAD file into the lines array
		self.__madFilePath = os.path.dirname(MADXfileName)
		fileName = os.path.basename(MADXfileName)
		#the initialize can be recursive if there are nested MAD files

		self.initialize(fileName)
		#-----------------------------------------------------------
		#The madxLines has all lines
		#Let's create values, elements, and accelerator lines arrays
		#-----------------------------------------------------------
		for line in self.__madxLines:
			if(line.getType() == _sequence.getType()):
				sequence = _sequence()
				sequence.parseLine(line.getLine())
				self._sequences.append(sequence)
			if(line.getType() == _variable.getType()):
				var = _variable()
				var.parseLine(line.getLine())
				self.__accValues.append(var)
			if(line.getType() == _element.getType()):
				elem = _element()
				elem.parseLine(line.getLine())
				self.__accElements.append(elem)
		#print "debug size values=",len(self.__accValues)
		#print "debug size elements=",len(self.__accElements)
		#print "debug size accLines=",len(self.__accLines)
		#-----------------------------------------------
		#replace all elem[key] substrings in elements by
		#variables
		#-----------------------------------------------
		accElemDict = {}
		for accElem in self.__accElements:
			accElemDict[accElem.getName()] = accElem
		accElemDictInit = accElemDict.copy()
		doNotStop = True
		while(doNotStop):
			doNotStop = False
			accElemDictCp = accElemDict.copy()
			#print "debug dict size=",len(accElemDictCp)
			for name,accElem in accElemDictCp.iteritems():
				kvs = accElem.getParameters()
				for key,val in kvs.iteritems():
					if val != None:
						resArr = StringFunctions.getElementKeys(val)
						if(len(resArr) == 0 and accElemDict.has_key(name)):
							del accElemDict[name]
						for [el,k] in resArr:
							doNotStop = True
							accElemInside = accElemDictInit[el]
							replVal = accElemInside.getParameters()[k]
							val = StringFunctions.replaceElementKeys(val,el,k,replVal)
					kvs[key] = val
			if(len(accElemDictCp) == len(accElemDict)):
				print "=========== Unresolved AccElements============"
				for name,accElem in accElemDictCp.iteritems():
					print "name=",name,"  params=",accElem.getParameters()
				print "=========== MAD File Problem ==============="
				print "=================STOP======================="
				sys.exit(1)
		#---------------------------------------------------------
		#Elements are ready!
		#Now let's substitute elements parameters in variables' expression.
		#---------------------------------------------------------
		for var in self.__accValues:
			val = var.getExpression()
			resArr = StringFunctions.getElementKeys(val)
			for [el,k] in resArr:
				accElemInside = accElemDictInit[el]
				replVal = accElemInside.getParameters()[k]
				val = StringFunctions.replaceElementKeys(val,el,k,replVal)
			var.setExpression(val)
		#-----------------------------------------------
		#now let's calculate all variables
		#-----------------------------------------------
		#replace all math cos,sin, etc by math.cos, math.sin, etc
		for var in self.__accValues:
			val = var.getExpression()
			val = StringFunctions.replaceMath(val)
			var.setExpression(val)
		#Then let's calculate numerical values.
		#They can be defined recursivelly, so we need iterations
		accVarDict = {}
		for var in self.__accValues:
			accVarDict[var.getName()] = var
		localValDict = {}
		doNotStop = True
		while(doNotStop):
			doNotStop = False
			accVarDictInner = accVarDict.copy()
			#print "debug variables dict. size=",len(accVarDictInner)
			for name,var in accVarDictInner.iteritems():
				str_in = var.getExpression()
				res,val = StringFunctions.calculateString(str_in.lower(),localValDict)
				if(res):
					localValDict[name.lower()] = val
					var.setValue(val)
					del accVarDict[name]
				else:
					doNotStop = True
			if(len(accVarDictInner) == len(accVarDict) and len(accVarDict) > 0):
				print "=========== Unresolved Variables============"
				for name,var in accVarDictInner.iteritems():
					print "name=",name,"  str=",var.getExpression()
				print "=========== MAD File Problem ==============="
				print "=================STOP======================="
				sys.exit(1)
		#-------------------------------------------
		# Now calculate all parameters in key,string_value
		# for accelerator elements
		#--------------------------------------------
		for accElem in self.__accElements:
			kvs = accElem.getParameters()
			kvNums = accElem.getNumParameters()
			for key,val in kvs.iteritems():
				val_out = None
				if val != None:
					res,val_out = StringFunctions.calculateString(val.lower(),localValDict)
					if(not res):
						print "=============MAD File problem ==============",
						print "Problem with acc. element:",accElem.getName()
						print "Parameter name:",key
						print "Can not calculate string:",val
						print "============ STOP =========================="
						sys.exit(1)
				kvNums[key] = val_out
		#---------------------------------------------
		#Let's create all lattice elements (old style)
		#---------------------------------------------
		for accElem in self.__accElements:
			lattElem = MADX_LattElement(accElem.getName(),accElem.getElementType())
			self.__lattElems.append(lattElem)
			kvs = accElem.getNumParameters()
			for key,val in kvs.iteritems():
				lattElem.addParameter(key,val)
		#----------------------------------------------
		#Let's create lattice lines (old style)
		#----------------------------------------------
		lattElemDict = {}
		for elem in self.__lattElems:
			lattElemDict[elem.getName()] = elem
		accLineDict = {}
		for accLine in self.__accLines:
			accLineDict[accLine.getName()] = accLine
		lattLineDict = {}
		for accLine in self.__accLines:
			lattLine = MAD_LattLine(accLine.getName())
			self.__lattLines.append(lattLine)
			lattLineDict[accLine.getName()] = lattLine
		for lattLine in self.__lattLines:
			name = lattLine.getName()
			accLine = accLineDict[name]
			childs = accLine.getComponents()
			for (child,sign) in childs:
				if(lattElemDict.has_key(child) and accLineDict.has_key(child)):
					print "=========== MAD File Problem ==============="
					print "Accelerator line and element have the same name:",child
					print "=================STOP======================="
					print "Lattice Line name=",name
					sys.exit(1)
				if((not lattElemDict.has_key(child)) and (not accLineDict.has_key(child))):
					print "=========== MAD File Problem ==============="
					print "Can not find Accelerator line and element with name:",child
					print "Lattice Line name=",name
					print "=================STOP======================="
					sys.exit(1)
				if(lattElemDict.has_key(child)):
					lattLine.addItem(lattElemDict[child],sign)
				else:
					lattLine.addItem(lattLineDict[child],sign)

		
						
		
	def makeDrift(self, downstreamelem):
	
		# Now we have to creat a drift between elements because MADX
		# sequence does not include drifts.
		
		seqlength = len(self._sequencelist)
		upstreamelemlength = 0
		downstreamelemlength = 0
		if(seqlength == 0):
			startpos = 0
			downstreamelemlength = float(downstreamelem.getParameter("l"))
		else:
			upstreamelem = self._sequencelist[seqlength-1]
			upstreamelemlength = float(upstreamelem.getParameter("l"))
			downstreamelemlength = float(downstreamelem.getParameter("l"))
			startpos = float(upstreamelem.getParameter("position")) + 0.5*upstreamelemlength
		endpos = float(downstreamelem.getParameter("position")) - 0.5*downstreamelemlength

		driftlength = endpos - startpos

		name = "Drift"# + "_" + str(seqlength)
		type = "drift"
		length = 0.0
		strength = 0.0
		
		if(driftlength < -1e-10):
			print "Warning: Drift between, ', upstreamelem.getName(), ' and ', downstreamelem.getName(), ' has negative length."
			print "Setting length to zero."
			lattElem = MADX_LattElement(name, type)
		else:
			lattElem = MADX_LattElement(name, type)
			lattElem.addParameter("l", driftlength)

		return lattElem
			
			

		
	def makeMultiPols(self,line_init):
		kls = []
		poles = []
		skews = []
		tokens = line_init.split("},")
		nvalues = len(tokens)
		for i in range(0,nvalues):
			subtokens = tokens[i].split(":={")
			name = subtokens[0]
			k = subtokens[1].split(',')
			if name == "knl":
				for i in range(len(k)):
					m = k[i].replace('}','')
					kls.append(float(m))
					poles.append(i)
					skews.append(0) 
			if name == "ksl":
				for i in range(len(k)):
					m = k[i].replace('}','')
					kls.append(float(m))
					poles.append(i)
					skews.append(1)
		return kls, poles, skews
		
	def makeAperture(self, downstreamelem):
	
		# Now we have to creat a aperture before and after the element with the MADX label aperture
		type = "apertype"
		name = "aperture"
		lattElem = MADX_LattElement("Aperture", name)
		lattElem.addParameter("l", 0.0)
		dim = downstreamelem.getParameter(name)
		shape_type = downstreamelem.getParameter(type)
		if shape_type == "circle":
			shape = 1
		elif shape_type == "ellipse":
			shape = 2
		elif shape_type == "rectangle":
			shape = 3
		else:
			print "======== Can not create elementwith type:",shape_type
			print "You have to fix the _teapotFactory, aperture class and madx_parser."
			print "Stop."
			sys.exit(1)
		lattElem.addParameter(name, dim)
		lattElem.addParameter(type, shape)
		return lattElem
			
	def getSequenceName(self):
		"""
			Method. Returns name of the sequence
			"""
		return self._sequencename

	def getSequenceList(self):
		"""
			Method. Returns list of elements in the sequence 
			"""
		return self._sequencelist

#-------------------------------------------------------------------
	def _init_string(self,str0):
		""" The method initializes the one string """
		#Now here 5 types of string
		# 0 - unknown type of the line
		# 1 - variables calculations
		# 2 - element definition
		# 3 - MADX sequence definition
		# 4 - call another nested MAD file
		#Delete spaces
		str0=re.sub(r'[ ]',"",str0)
		tp = self._findLineType(str0)
		if tp == 0:
			#print "StrType =0 :",str0
			return
		if tp == 1:
			#print "StrType =1 :",str0
			madxLine = _madxLine()
			madxLine.setLine(str0)
			madxLine.setType("variable")
			self.__madxLines.append(madxLine)
		if tp == 2:
			#print "StrType =2 :",str0
			madxLine = _madxLine()
			madxLine.setLine(str0)
			madxLine.setType("element")
			self.__madxLines.append(madxLine)
		if tp == 3:
			#print "StrType =3 :",str0
			madxLine = _madxLine()
			madxLine.setLine(str0)
			madxLine.setType("sequence")
			self.__madxLines.append(madxLine)
		if tp == 4:
			#print "StrType =4 :",str0
			fl = self._parseNestedFileLine(str0)
			self.initialize(fl)

	def _findLineType(self,line):
		""" Return type of the string """
		#Now here 5 types of string
		# 0 - unknown type of the line
		# 1 - variables calculations
		# 2 - element definition
		# 3 - MADX sequence definition
		# 4 - call another nested MAD file
		# 5 - location of the element (name, at = location;)
		stype = 0
		for word in self.__ingnoreWords:
			if(line.lower().find(word) == 0):
				return stype
		t_match = re.search(r'.*:=.*',line)
		if t_match:
			stype = 1
		t_match = re.search(r'.*=.*',line)
		if stype != 1:
			t_match = re.search(r'[\w]* *:.*',line)
			if t_match:
				stype = 2
		t_match = re.search(r'.*:.*sequence,.*',line.lower())
		if t_match:
			stype = 3
		t_match = re.search(r' *call *file *=',line.lower())
		if t_match:
			stype = 4
		#deal with constant defenition like AAA = 1.0
		if(stype == 0):
			t_match = re.search(r'.*=.*',line)
			if t_match:
				stype = 1

		t_match = re.search(r' *at=*',line.lower())
		if t_match:
			stype = 5
			print("location found")
		return stype

#===============================================================
# 
#=============================================================== 
class StringFunctions:
	"""
	This class defines the set of static string functions.
	"""

	def replaceElementKeys(self, str_in, elem, key, value):
		"""
		Method. It will replace elem[key] in the string expression of this variable.
		"""
		new_val = r'(' + str(value) + r')'
		s = elem+r'\['+key+r'\]'
		patt = re.compile(s)
		str_out = re.sub(patt,new_val,str_in)
		return str_out

	replaceElementKeys = classmethod(replaceElementKeys)

	def getElementKeys(self, str_in):
		"""
		Method. It returns the set of [element,key] pairs for input string.
		"""
		res = []
		patt=re.compile(r'[\w]*\[[\w]*\]')
		s_name_key = re.findall(patt,str_in)
		if len(s_name_key) == 0:
			return res
		patt_elem = re.compile(r'[\w]*(?=\[)')
		patt_key = re.compile(r'(?<=\[)[\w]*(?=\])')
		for s in s_name_key:
			elem = re.findall(patt_elem,s)[0]
			key = re.findall(patt_key,s)[0]
			res.append([elem,key])
		return res

	getElementKeys = classmethod(getElementKeys)

	def calculateString(self, str_in, localDict):
		"""
		Method. It returns a tuple (True,value) if
		the expression can be evaluated and (False,None) otherwise.
		"""
		try:
			val = eval(str_in,globals(),localDict)
			return (True, val)
		except:
			return (False, None)

	calculateString = classmethod(calculateString)

	def replaceMath(self, str_in):
		"""
		Method. It replaces math symbols
		to make them readable for python eval().
		"""
		#replace .e by .0e
		str_out = re.sub("\.e",".0e",str_in)
		#check the math operatons
		str_out = re.sub("sin\(","math.sin(",str_out)
		str_out = re.sub("SIN\(","math.sin(",str_out)
		str_out = re.sub("cos\(","math.cos(",str_out)
		str_out = re.sub("COS\(","math.cos(",str_out)
		str_out = re.sub("tan\(","math.tan(",str_out)
		str_out = re.sub("TAN\(","math.tan(",str_out)
		str_out = re.sub("exp\(","math.exp(",str_out)
		str_out = re.sub("EXP\(","math.exp(",str_out)
		str_out = re.sub("log\(","math.log(",str_out)
		str_out = re.sub("LOG\(","math.log(",str_out)
		str_out = re.sub("acos\(","math.acos(",str_out)
		str_out = re.sub("ACOS\(","math.acos(",str_out)
		str_out = re.sub("asin\(","math.asin(",str_out)
		str_out = re.sub("ASIN\(","math.asin(",str_out)
		str_out = re.sub("atan\(","math.atan(",str_out)
		str_out = re.sub("ATAN\(","math.atan(",str_out)
		str_out = re.sub("sqrt\(","math.sqrt(",str_out)
		str_out = re.sub("SQRT\(","math.sqrt(",str_out)
		return str_out

	replaceMath = classmethod(replaceMath)

#====================================================================
class _madxLine:
	"""
	Class: A line of a MAD file. It can be one of three types:
	variable, accelerator sequence, accelerator element and positioning of the acc elem.
	"""
	def __init__(self):
		"""
		Constructor. Creates a blank MADXfile line class instance.
		"""
		self.__line = ""
		self.__type = None

	def setLine(self, line):
		"""
		Method. Sets a MAD file line class instance.
		"""
		self.__line = line

	def getLine(self):
		"""
		Method. Returns a MAD file line class instance.
		"""
		return self.__line

	def setType(self, ltype):
		"""
		Method. Sets a MAD file line type.
		"""
		self.__type = ltype

	def getType(self):
		"""
		Method. Returns a MAD file line type.
		"""
		return self.__type
#====================================================================

class _sequence:

	"""The MADX accelerator class. It also keeps initial string (line) from MADX file."""

	def __init__(self):
		self._name = None
		self._expression = ""
		self._components = []

	def getType():
		"""
		Method. It is static method of this class.
		It returns the name of the type.
		"""
		return "sequence"

	getType = staticmethod(getType)

	def setName(self, name):
		self._name = name

	def getName(self):
		return self._name

	def setExpression(self,line):
		self._expression = line

	def getExpression(self):
		return self._expression

	def getComponents(self):
		"""Method. It returns the set with components"""
		return self._components

	def parseLine(self,line_init):
		"""
		Method. It does the first parsing of the initial string.
		Parses the MADX file line with lattice sequence definition.
		"""
		#define name of new line
		patt = re.compile(r'[\w]+(?=:)')
		name = re.findall(patt,line_init)[0]
		self.setName(name)
		patt = re.compile(r'(?<==).*')
		s_def = re.findall(patt,line_init)[0]
		self.setExpression(s_def)
		patt = re.compile(r'(?<=\(|,).+?(?=,|\))')
		item_names = re.findall(patt,s_def)
		#=========================================
		#deal with the N*name expressions
		#=========================================
		item_names_new = []
		for it_in in item_names:
			sign = +1
			it = it_in
			if(it.find("-") == 0):
				sign = -1
				it = it_in[1:]			
			patt = re.compile(r'[\d]+?(?=\*)')
			n_rep = re.findall(patt,it)
			if len(n_rep) > 0:
				n_rep = int(n_rep[0])
				patt = re.compile(r'(?<=\*)[\w]+')
				s=re.findall(patt,it)
				s=s[0]
				for i in range(1,n_rep+1):
					item_names_new.append(s)
			else:
				item_names_new.append(it)
		self._components = item_names_new

#====================================================================
class _variable:
	"""
	Class. Holds MAD variables.
	"""
	def __init__(self):
		"""
		Constructor. Creates empty MAD variable.
		"""
		self._name = None
		self._expression = ""
		self._value = None

	def getType():
		"""
		Method. Static method of this class.
		Returns the name of the type.
		"""
		return "variable"

	getType = staticmethod(getType)


	def getValue(self):
		"""
		Method. It returns the numerical value of this variable.
		"""
		return self._value

	def setValue(self, val):
		"""
		Method. It sets the numerical value of this variable.
		"""
		self._value = val

	def setName(self, name):
		self._name = name

	def getName(self):
		return self._name

	def setExpression(self,line):
		self._expression = line

	def getExpression(self):
		return self._expression

	def parseLine(self,line_init):
		"""
		Method. It does the first parsing of the initial string.
		"""
		#divide string onto two parts: name of value and value
		patt = re.compile(r'(?<=:=).*')
		s_val = re.findall(patt,line_init)
		if(len(s_val) > 0):
			patt = re.compile(r'.*(?=:=)')
			s_name  = re.findall(patt,line_init)
			s_val  = s_val[0]
			s_name = s_name[0]
		else:
			#deal with const defenition like AA = 1.2
			patt = re.compile(r'(?<==).*')
			s_val = re.findall(patt,line_init)
			patt = re.compile(r'.*(?==)')
			s_name  = re.findall(patt,line_init)
			s_val  = s_val[0]
			s_name = s_name[0]
		self.setName(s_name)
		self.setExpression(s_val)

#====================================================================
class _element:
	"""
	The MAD element class. It also keeps initial string (line) from MAD file.
	"""

	def __init__(self):
		self._name = None
		self._elementType = None
		self._parameters = {}
		self._numParameters = {}

	def getType():
		"""
		Method. It is static method of this class.
		It returns the name of the type.
		"""
		return "element"

	getType = staticmethod(getType)

	def setName(self, name):
		self._name = name

	def getName(self):
		return self._name

	def setElementType(self, etype):
		self._elementType = etype

	def getElementType(self):
		return self._elementType

	def getParameters(self):
		"""
		Method. It returns the dictionary with key:string pairs
		"""
		return self._parameters

	def getNumParameters(self):
		"""
		Method. It returns the dictionary with {key:(numeric value)} pairs
		"""
		return self._numParameters

	def parseLine(self,line_init):
		"""
		Method. It does the first parsing of the initial string.
		"""
		#divide string onto three parts: name,type and key-values of element
		patt = re.compile(r'[\w]+(?=:)')
		name = re.findall(patt,line_init)[0]
		patt = re.compile(r'(?<=:)[\w]+?(?=\s|,)')
		type = re.findall(patt,line_init+",")[0]
		self.setName(name)
		self.setElementType(type)
		patt = re.compile(r'(?<=\w),.*')
		s_key_vals = re.findall(patt,line_init)
		#==================================================
		#now we have name and type, let's get key-val pairs
		#==================================================
		if len(s_key_vals) < 1:
			return
		s_key_vals = s_key_vals[0] + ","
		patt = re.compile(r'(?<=,).+?(?=,)')
		s_key_vals = re.findall(patt,s_key_vals)
		for s_key_val in s_key_vals:
			patt = re.compile(r'[\w]+?(?==)')
			key =  re.findall(patt,s_key_val)
			if(len(key) == 0):
				#we have stand alone key, e.g. T1
				key = s_key_val
				val = None
			else:
				key = key[0]
				patt = re.compile(r'(?<==).*')
				val = re.findall(patt,s_key_val)[0]
			self._parameters[key] = val
			self._numParameters[key] = val

#====================================================================

#STOP parsing a MAD file if there is a start of mad commands

	def old_parse(self,MADXfileName):
		
		self.__init__()
		
		#1-st stage read MAD file into the lines array
		self.__madFilePath = os.path.dirname(MADXfileName)
		fileName = os.path.basename(MADXfileName)
		#the initialize can be recursive if there are nested MAD files

		fl = open(os.path.join(self.__madFilePath, fileName))

		str_local = ""
		aper_warning = 0
		for str in fl.readlines():
			#print str
			#take off the ";" at end of each line
			#print str
			
			str0 = str
			if str0.rfind(";") > 0:
				str_local = ""
				for i in xrange(str0.rfind(";")):
					str_local = "".join([str_local,str0[i]])
			str_local.strip()
			#check if the line is a comment
			if str.find("!") == 0:
				str_local = ""
				continue
			#check the empty line
			if str == "":
				str_local = ""
				continue
			#Parse element or sequence definition. 
			if re.search(r'[\w]* *:.*',str_local):
				if(str_local.rfind("sequence") >=0):
					tmp = "".join(str_local.split())
					tokens = tmp.split(":")
					self._sequencename = tokens[0]
					self._sequencelength = eval(tokens[1].split("l=")[1].rstrip(";"))
				else:
					elem = self.parseElem(str_local)
					self._accElemDict[elem.getName()] = elem
			#Add elements to the sequence list.
			if(str_local.rfind("at") >= 0):
				if(self._sequencename==""):
					print "Warning, adding elements to sequence with no name."
				tmp = "".join(str_local.split())
				tokens = tmp.split(",")
				elem_name = tokens[0]
				position = eval((tokens[1].lstrip("at=")).rstrip(";"))
				latt_elem = self._accElemDict[elem_name]
				length = float(latt_elem.getParameter("l"))
				latt_elem.addParameter("position", position)
				if(latt_elem.hasParameter("apertype")!=0):
					latt_aper_entry = self.makeAperture(latt_elem)
					latt_aper_entry.addParameter("position", position-length/2.0)
					latt_aper_exit = self.makeAperture(latt_elem)
					latt_aper_exit.addParameter("position", position+length/2.0)
					latt_drift = self.makeDrift(latt_elem)
					self._sequencelist.append(latt_drift)
					self._sequencelist.append(latt_aper_entry)
					self._sequencelist.append(latt_elem)
					self._sequencelist.append(latt_aper_exit)
					aper_warning = aper_warning + 2
				else:
					latt_drift = self.makeDrift(latt_elem)
					self._sequencelist.append(latt_drift)
					self._sequencelist.append(latt_elem)

		if aper_warning >= 1:
			print "Warning, adding", aper_warning ,"aperture nodes to the teapot lattice. That will slow down the simluation."
			print "If the lost of particles on the aperture is not necessary, please use a madx file without the aperture labels."
		
		#If the last element is not at the end of the lattice, make a drift		
		nelems = len(self._sequencelist)
		#print 'number of elems', nelems
		if(nelems == 0):
			print "Warning: Creating empty lattice."
			sys.exit(1)
		endlength = float(self._sequencelength) - (float(self._sequencelist[nelems-1].getParameter("position")) + 0.5*float(self._sequencelist[nelems-1].getParameter("l")))
		if(endlength < -1e-10):
			print "Warning: Lattice parsing resulted in a lattice with length longer than specified by sequence command."
		if(endlength > 0):
			latt_drift = MADX_LattElement("end drift", "drift")
			latt_drift.addParameter("l", endlength)
			self._sequencelist.append(latt_drift)

	def old_parseElem(self,line_init):
		"""
			Method. It does the first parsing of the initial string.
		"""
		name = ""
		type = ""
		length = 0.0
		strength = 0.0
		
		match = re.search(r'[\w.-]+:={.*',line_init)   # search for knl, ksl, APERTURE, ...
		if match is not None:
			line_init = line_init[:-len(match.group())-1]
			match_multi = re.search(r'k+[\w,-]+:={.*}',match.group())
			match_aper = re.search(r'aperture+.*',match.group())
			if match_aper is not None:
				apertype = re.search(r'apertype+.*',match.group())
				apertype = apertype.group()
				if match_multi is not None:
					aper = match_aper.group()[:-len(match_multi.group()) - len(apertype) - 2]
				else:
					aper = match_aper.group()[: - len(apertype) - 1]
		line_init = ''.join(line_init.split())  # remove all whitespaces character
		tokens = line_init.split(",")
		print(line_init)
		nvalues = len(tokens)
		if nvalues >= 1:
			subtokens = tokens[0].split(":")
			name = subtokens[0]
			type = subtokens[1].strip()
			lattElem = MADX_LattElement(name, type)
			for i in range(1, nvalues):
				subtokens = tokens[i].split(":=")
				lattElem.addParameter(subtokens[0],float(subtokens[1]))
			if match is not None:
				if match_multi is not None:
					kls, poles, skews = self.makeMultiPols(match_multi.group())
					lattElem.addParameter("poles", poles)
					lattElem.addParameter("skews", skews)
					lattElem.addParameter("kls", kls)
				if match_aper is not None:
					apertype = apertype.split("=")
					dim = []
					toke = aper.split("={")
					subtok = toke[1].split(',')
					for i in range(len(subtok)):
						m = subtok[i].replace('}','')
						dim.append(float(m))
					lattElem.addParameter("aperture", dim)
					lattElem.addParameter("apertype", apertype[1])
					
		else:
			lattElem = MADX_LattElement(name, type)
			print "Warning: Returning empty lattice element type."

		#Make sure there is always a length parameter.
		if(lattElem.hasParameter("l")==0):
			lattElem.addParameter("l",0.0)
				
		return lattElem
