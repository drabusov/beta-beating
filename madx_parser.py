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
	"""Class. Represents an arbitrary element in the lattice"""
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
			print "class MAD_LattElement, method getParameter"
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
class _variable:
	"Class. Holds MADX variables."
	def __init__(self):
		"Constructor. Creates empty MAD variable."
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

class MADX_Parser:
	""" MAD parser """
	
	def __init__(self):
		""" Create instance of the MAD_Parser class """
		self._accElemDict = {}
		self._varDict = {} # dict of variables
		self._sequencename = ""
		self._sequencelength = ""
		self._sequencelist = []
		#the lines to ignore will start with these words
		self.__ingnoreWords = ["title","beam", "none"]

	def parse(self,MADXfileName):
		
		self.__init__()
		
		#1-st stage read MAD file into the lines array
		self.__madFilePath = os.path.dirname(MADXfileName)
		fileName = os.path.basename(MADXfileName)
		#the initialize can be recursive if there are nested MAD files
		fl = open(os.path.join(self.__madFilePath, fileName))
		str_local = ""
		aper_warning = 0

		for str in fl.readlines():
			loc_flag_elem = True
			#print str
			#take off the ";" at end of each line
			#print str
			
			str0 = str
			str0 = "".join(str.split()) # remove all spaces!
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
			#Parse element/variable/sequence definition. 
			if re.search(r'[\w]* *:.*',str_local):
				if(str_local.rfind("sequence") >=0):
					tokens = str_local.split(":") # is it possible to have more than two tokens here?
					self._sequencename = tokens[0]
					tmp_str = ":".join(tokens[1:])
					aux = [x.split("l=")[-1] for x in tmp_str.split(",") if "l=" in x]
					self._sequencelength = eval(aux[0])
				else:
					print(str_local, loc_flag_elem)
					if ":=" in str_local:
						if ":" not in str_local.split(":="):
							# we have a MADX variable here! 
							var = _variable()
							var.parseLine(str_local)
							self._varDict[var.getName()] = var
							loc_flag_elem = False # in order to avoid the parsing of element
#----------------------------------------------------------------------------------
					if loc_flag_elem:
						if "at=" in str_local:
							tokens = str_local.split(",") 
							tmp = ",".join([x for x in tokens if "at=" not in x]) # remove location from the string (is it necessary?)
						else:
							tmp = str_local
						elem = self.parseElem(tmp) # elem can be defined directly at the location!
						self._accElemDict[elem.getName()] = elem
			#Add elements to the sequence list.
			if(str_local.rfind("at") >= 0):
				if(self._sequencename==""):
					print "Warning, adding elements to sequence with no name."

				if ":" not in str_local:			
					tokens = str_local.split(",") # it can cause problems with {.,.} def inside the line
					elem_name = tokens[0]	# elem name 100% is the first position, otherwise user did a mistake
					tmp = tokens[1:]
					aux = [x.split("at=")[-1] for x in tmp if "at=" in x]				
					position = eval(aux[0]) 

				else:		
					tokens = str_local.split(":") 
					elem_name = tokens[0]	# elem name 100% is the first position, otherwise user did a mistake
					tmp_str = "".join(tokens[1:])
					tmp = tmp_str.split(",")  # it can cause problems with {.,.} def inside the line
					aux = [x.split("at=")[-1] for x in tmp if "at=" in x]				
					position = eval(aux[0]) 
					
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

		print(self._varDict.keys())
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
			
			
	def parseElem(self,line_init):
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
		nvalues = len(tokens)
		if nvalues >= 1:
			subtokens = tokens[0].split(":")
			name = subtokens[0]
			type = subtokens[1].strip()
			# elem can be defined as a child node
			if type not in self._accElemDict.keys():			
				lattElem = MADX_LattElement(name, type)
			else:	
				type_upd = self._accElemDict[type].getType()
				lattElem = MADX_LattElement(name, type_upd)
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
					kls.append(eval(m))
					poles.append(i)
					skews.append(0) 
			if name == "ksl":
				for i in range(len(k)):
					m = k[i].replace('}','')
					kls.append(eval(m))
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

#STOP parsing a MAD file if there is a start of mad commands
