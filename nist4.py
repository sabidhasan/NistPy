# -*- coding: utf-8 -*-

##TO--DO: error catching - use try except stuff rather than checking for empty lists, etc.
##TO--DO: test the code with test cases, write documentation...

#!/usr/bin/env python
#!/Users/abidhasan/anaconda/envs/rdkkit4
from bs4 import BeautifulSoup
import requests
import re
import pandas as pd
import matplotlib.pyplot as plt
 
class NistCompoundById(object):
	'''This class creates a NIST compound using a provided NIST code.
	NIST identifying codes consist of a letter followed by several numbers'''

	#A search result is NistResult object. It's instantiated is False. This is for user 
	#to see whether the object has been instantaited in a friendly manner. Really you can
	#just do type(object); if it's NistCpdById then it is instantiated; if NistResult, it's not.
	INSTANTIATED = True

	@staticmethod
	def extract_CAS(cas_string):
		'''Extracts a CAS number using RegEx from a given string cas_string'''
		reg = re.search(r'\d{2,7}-\d{2}-\d', cas_string)
		if reg:
			return reg.group(0)
		return None
			
	def __init__(self, id, souped_page=None):
		'''souped_page is an already beautiful souped page for a single compound (search page
		on website returns the compound page onyl one hit is found, so there is no point in
		reinstantiating the compound and redownlaoding the webpage, so if souped_page is passed,
		the init function can use that directly.
		
		self.name
		self.id
		self.synonyms
		self.cas
		self.InChI
		self.mol_weight
		self.formula
		self.available_properties
		self.available_spectra
		self.spectra		
		
		self.thermo_data
		self.phase_data
		self.henrys_data
		self.gas_phase
		'''
		if souped_page is not None:
			cpd_url = "http://webbook.nist.gov/cgi/cbook.cgi?ID=%s" % id
			cpd_soup = download_page(cpd_url)
		
			if cpd_soup is None:
				#download failed
				return None
		else:
			#we are given a souped page already!
			cpd_soup = souped_page	
		
		#check for erroneous number
		if "Registry Number Not Found" in cpd_soup:
			print "<Compound %s not found in registry>" % id
			return None

		#Set compound name
		self.name = cpd_soup.find('h1').text
		self.id = id
		
		#Find other compound properties		
		cpd_raw_props = cpd_soup.find_all('li')
		
		self.synonyms = None
		self.cas = None
		self.InChI = None
		self.mol_weight = None
		self.formula = None
		self.available_properties = {}
		self.available_spectra = {}		
		
		for prop in cpd_raw_props:
			property_text = prop.get_text()
			
			if "Formula" in property_text:
				self.formula = property_text[9:]
			elif "Molecular weight" in property_text:
				self.mol_weight = float(property_text[18:])
			elif "InChI" in property_text:
				self.InChI = property_text
			elif "CAS" in property_text:
				self.cas = self.extract_CAS(property_text)
			elif "Other names:" in property_text:
				#synonyms should start at "1" because the first one is literally "Synonyms:"
				self.synonyms = property_text.split(";")[1:]
			elif "Other data available:" in property_text:
				all_properties = ['Gas phase thermochemistry data', 'Phase change data', "Henry's Law data", \
				'Gas phase ion energetics data', 'IR Spectrum', 'Mass spectrum (electron ionization)', \
				'Constants of diatomic molecules', 'UV/Visible spectrum']
				spectral_properties = ['IR Spectrum', 'Mass spectrum (electron ionization)', 'UV/Visible spectrum']
				
				#Loop through the properties
				for i in prop.find_all('a'):
					if i.text in all_properties:
						self.available_properties[i.text] = i.get('href')
					if i.text in spectral_properties:
						self.available_spectra[i.text] = i.get('href')

				######TESTING CODE REMOVE LATER
				#len_diff = len(prop.find_all('a'))-  len(self.available_properties)
				#print "%s properties were not used for %s..." % (len_diff, self.name)	

		self.spectra = {'uv': None, 'mass': None, 'ir': None}
		
	def __str__(self):
		return 'Compound %s' % self.name
		
	def __repr__(self):
		return "NistCompoundById('%s')" % self.id
			
	def _get_data(self, parameter, table=0):
		'''This (hidden) function recives a parameter (eg. gas phase or henry's law
		and then parses the compounds webpage for that data. Returns pandas dataframe
		int table parameter decides which table to pick on the page (0 is first)'''
		
		#check to see if data is available
		if self.available_properties.get(parameter) is None:
			print "<%s data is not available for compound %s>" % (parameter, self.name)
			return None

		#if available, then get its page
		prop_url = "http://webbook.nist.gov/" + self.available_properties[parameter]
		prop_soup = download_page(prop_url)
		
		#Get the rows in the table
		data_dict = { 'Quantity': [], 'Value': [], 'Unit': [] }
		
		#use [table] to get the n-th table on the page
		prop_table = prop_soup.find_all("table", { "class" : "data" })[table].find_all('tr')
			
		#Loop through the rows, 
		for row in prop_table:
			if len(row.find_all('td')) > 0:
				data_dict['Quantity'].append(row.find_all('td')[0].text)
				data_dict['Unit'].append(row.find_all('td')[2].text)
				data_dict['Value'].append(row.find_all('td')[1].text)
		return data_dict	
	
	def get_thermo_data(self):
		thermo = NistCompoundById._get_data(self, "Gas phase thermochemistry data")
		self.thermo_data = pd.DataFrame(thermo)
		return self.thermo_data
	
	def get_phase_data(self):
		phase = NistCompoundById._get_data(self, "Phase change data")
		self.phase_data = pd.DataFrame(phase)
		return self.phase_data

	def get_henrys_law(self):
		henrys = NistCompoundById._get_data(self, "Henry's Law data")

		#Value is actually methods for Henry's Law, so not needed. Delete from dict
		del henrys['Unit']
		henrys = pd.DataFrame(henrys)		
		henrys.columns = ['k°H (mol/kg*bar', 'd(ln(kH))/d(1/T) (K)']
		self.henrys_data = henrys
		return self.henrys_data

	def get_gas_ion_energetics(self):
		gas_phase = NistCompoundById._get_data(self, "Gas phase ion energetics data", 1)
		del gas_phase['Unit']
		del gas_phase['Value']

		self.gas_phase = pd.DataFrame(gas_phase)
		return self.gas_phase
		
	def get_spectra(self, spectral_type):
		'''This method fetches spectral data for the given spectrum type provided as string (ir, mass, or uv)
		a Pandas DataFrame containing the data if successful
		'''

		spectral_defs = {'mass' : 'Mass spectrum (electron ionization)',
		'uv' : 'UV/Visible spectrum',
		'ir' : 'IR Spectrum'}
		
		#check to see if data is available
		if self.available_properties.get(parameter) is None:
			print "<%s data is not available for compound %s>" % (spectral_type, self.name)
			return None

		#check to see validity of spectral type argument
		if spectral_defs.get(spectral_type) is None:
			print "<Invalid spectral type specified for compound %s. Supported types are 'uv', 'ir' and 'mass'>" % self.name
			return None
		
		#if spectral data already exists, then return in; no need to download twice
		if self.spectra[spectral_type] is not None:
			return self.spectra[spectral_type]
		
		#'uv', 'mass', 'ir' are user-friendly; real DB uses Mass Spectrum (electrion ionization), etc. so 'parameter' stores that
		parameter = spectral_defs[spectral_type]
		

		#if data is available, then get the resulting page
		spectral_url = "http://webbook.nist.gov" + self.available_properties[parameter]
		spectral_page = requests.get(spectral_url)
		spectral_page_content = spectral_page.content
		
		#Check for error
		if spectral_page.status_code != 200:
			print "<Error downloading web page. Error status code %s>" % (spectral_page.status_code)
			return None
		
		#Find and download the spectrum
		if spectral_type == "uv":
			dl_link = "http://webbook.nist.gov/cgi/cbook.cgi?JCAMP=%s&amp;Index=0&amp;Type=UVVis" % self.id
		elif spectral_type == "ir":
			#TO--DO: there are other IR spectra as well!
			dl_link = "http://webbook.nist.gov/cgi/cbook.cgi?JCAMP=%s&Index=0&Type=IR" % self.id			
		elif spectral_type == "mass":	
			dl_link = "http://webbook.nist.gov/cgi/cbook.cgi?JCAMP=%s&Index=0&Type=Mass" % self.id
					   
		#Create and read into dictionary
		xs = []
		ys = []
		raw_data = ""
		spectral_file_data = requests.get(dl_link, stream=False)
		for chunk in spectral_file_data.iter_content(chunk_size=1024*1024*1024):
			raw_data += chunk

		for line in raw_data.split('\n'):
			if len(line) > 0:
				try:
					#This catches any lines that commence with non data characters
					int(line[0])
				except ValueError:
					continue
				
				if spectral_type == 'uv':
					#spectrum is is UV
					if len(line.split(',')) == 2:
						xs.append(float(line.split(',')[0]))
						ys.append(float(line.split(',')[1]))
				elif spectral_type == 'mass':
					#spectra is mass, so data is laid out differently
					for data in line.split(' '):
						xs.append(float(data.split(',')[0]))
						ys.append(float(data.split(',')[1]))
				else:
					#spectral type is IR. TO--DO: I dont know what the many columns mean, so only first columns is obtained for now
					data = line.split(' ')
					xs.append(float(data[0]))
					ys.append(float(data[1]))
		
		if spectral_type == 'ir': 
			#normalize the y data, so that it's relative transmission!
			ys = [item/float(max(ys)) for item in ys]
					
		self.spectra[spectral_type] = pd.DataFrame({'x': xs, 'y': ys})
		#print "<Successfully obtained %s spectral data>" % spectral_type
		return self.spectra[spectral_type]
	
	@property
	def list_available_properties(self):
		return self.available_properties.keys()

	@property
	def list_available_spectra(self):
		return self.available_spectra.keys()
		
class NistPlot(object):
	'''NistPlot class allows plotting of spectral data. Expects compounds_to_plot (Python list, or a NistCompoundById object,
	a spectral type ('uv', 'ir' or 'mass') and an optional style tuple containing one or multiple styles
	self.compounds
	self.spectral_type
	self.style
	'''
	def __init__(self, cpds_to_plot, spectral_type, style, graph_title=''):
		#check for provided spectral type (if not valid then throw error)		
		if not(spectral_type in ['mass', 'uv', 'ir']):
			print "<Invalid spectral type specified %s; mass, ir or uv are valid spectral types>" % spectral_type
			return None
		#set spectral_type (uv, ir etc)
		self.spectral_type = spectral_type
		
		if not(isinstance(style, tuple)):
			print "<Graph style is expected to be tuple; got %s instead>" % type(style)
			return None
		self.style = style

		#check for object identity (list or NistCompoundById are supported)
		if isinstance(cpds_to_plot, list):
			#Check for length
			if len(cpds_to_plot) == 0:
				print "<List has no valid compound to plot>"
			#Check for each member
			for item in cpds_to_plot:
				if not(isinstance(item, NistCompoundById)):
					print "<Object %s is not valid for plotting. Expected type NistCompoundById, got %s>" % (item, type(item))
					return None
			#everything is good so let's set self.compounds
			self.compounds = cpds_to_plot
		elif isinstance(cpds_to_plot, NistCompoundById):
			#make into list so it can be iterated later
			self.compounds = [cpds_to_plot]
		else:
			print "<No list of compounds provided to plot. NistPlot supports Python lists or NistCompoundById objects>"
			return None

		#loop over cpds to see if they are all valid Nist Compounds and if there are missing spectral data
		for counter, compound in enumerate(self.compounds):
			if not(isinstance(compound, NistCompoundById)):
				print "<Compound %s is not a valid NistPy object; got class %s instead>" % (compound, type(compound))
				self.compounds.remove(compound)
				#prevent further checking on this cpd!
				continue
			
			if compound.spectra[spectral_type] is None:
				#No data found, so let's grab the spectra
				if compound.get_spectra(spectral_type) is None:
					print "<Spectrum could not be loaded for compound %s>" % compound.name
					continue
		
		self.title = graph_title
		
	def show(self):
		#TO--DO: Which of clf and close is actually useful here?
		plt.clf()
		plt.close()
		for counter, compound in enumerate(self.compounds):
			curr_style = self.style[int(counter % len(self.style))]
			#if mass specturm then show bar otherwise show plot
			if self.spectral_type == 'mass':
				plt.bar(compound.spectra[self.spectral_type]['x'], compound.spectra[self.spectral_type]['y'], 0.25, color=self.style[0][0])
			else:
				
				plt.plot(compound.spectra[self.spectral_type]['x'], compound.spectra[self.spectral_type]['y'], curr_style)
		#TO--DO::: SHOW THE TITLE OF THE GRAPH! LEGEND!!! AXES!!!
		plt.show()	
		
	def __add__(self, other):
		if not(isinstance(other, (NistPlot, NistCompoundById, list) )):
			raise TypeError("<Unsupported operand for concatenation; only NistCompoundById objects, Python lists \
			containing NistCompoundById objects, and NistPy objects are supported>")
			return None
		
		concat_list = self.compounds
		
		new_style = self.style		
				
		if isinstance(other, NistPlot):
			#it's a nistplot object so just add the lists togehter
			concat_list += other.compounds
			new_style += other.style
			
			#if there is a spectral type mismatch then say so!
			if self.spectral_type != other.spectral_type:
				print "Warning: type mismatch detected; one graph is for %s and the other for %s" % (self.spectral_type, other.spectral_type)
				#TO--DO: MAKE INTERWEAVING TUPLE FOR STYLES. IF THERE IS >2 STLYES, THEY SHOULD INTERWEAVE
		elif isinstance(other, NistCompoundById):	
			concat_list.append(other)
		else:
			#it's a list, so lets parse it
			for cpd in other:
				if isinstance(cpd, NistCompoundById):
					concat_list.append(cpd)
				else:
					print "<Some item(s) in supplied list are not valid NistPy object(s). Expected NistCompoundById class, got %s>" % type(cpd)
		
		self.title = self.spectral_type.title() + ' Spectrum for %s compounds' % len(list(set(concat_list)))
		
		return NistPlot(list(set(concat_list)), self.spectral_type, new_style, self.title + other.title)
	
	def __str__(self):
		return 'Plot for %s' % self.compounds
		
	def __repr__(self):
		return "NistPlot(%s, '%s', %s, '%s')" % (self.compounds, self.spectral_type, self.style, self.title)
		
def allNistCompounds():
	all_compounds = []
	with open('species.txt') as nist_file:
		for line in nist_file.readlines():
			all_compounds.append(tuple(line.split('\t')[:2]))
	return all_compounds

def download_page(url):
	'''This function downloads a page given and URL. If fail, it returns a status-code'''	
	page = requests.get(url)
	page_content = page.content
		
	#test for erroneous download
	if page.status_code != 200:
		print "<Error downloading web page. Error status code %s>" % cpd_page.status_code
		return None
		
	return BeautifulSoup(page_content, 'html.parser')

class NistResult(object):
	'''NistResult class is used to return search results. This class is like a
	NistCompoundById but only has id and name. It has the instantiate method
	to elaborate and return NistCompoundById objects'''

	INSTANTIATED = False
	
	def __init__(self, id_num, name, souped_page=None):
		self.id = id_num
		self.name = name
		self.page = souped_page

	def instantiate(self):
		#This will instantiate the compound
		if self.page is None:
			self.page = download_page("http://webbook.nist.gov/cgi/cbook.cgi?ID=%s" % self.id)
			
		return NistCompoundById(self.id, souped_page=self.page)

	def __str__(self):
		return self.name

	def __repr__(self):
		return 'NistResult(%s, %s)' % (self.id, self.name)

def search_nist(options=None, instantiate=False, **kwargs):
	'''Search NIST database for string search_str using the search keyword [formula, name, inchi key, CAS number, 
	structure]. If instantiate is selected, then it returns NistCompoundById object (which has a bunch of data within it)
	Otherwise it returns a NistResult object, which only has the name and ID code and must be called with the
	.instantiate() method before further use'''

	#Error checking
	if len(kwargs) != 1:
		print "<Expected one keyword, obtained %s>" % len(kwargs)
		return None
	
	search_type, search_str = kwargs.keys()[0], kwargs.values()[0]
		
	#Check for value of kwarg	
	types = {
	'formula': "http://webbook.nist.gov/cgi/cbook.cgi?Formula=%s" % search_str,
	'name' : "http://webbook.nist.gov/cgi/cbook.cgi?Name=%s" % search_str,
	'inchi' : "http://webbook.nist.gov/cgi/cbook.cgi?InChI=%s" % search_str,
	'cas' : "http://webbook.nist.gov/cgi/cbook.cgi?ID=%s" % search_str,
	}
	if not(search_type) in types:
		print "<Got unexpected search type. Only these search types allowed %s>" % types.keys()
		return None
		
	#Download the search page, and check to see whether it's a single compound
	search_page = download_page(types[search_type])
	#will be returned to user
	search_results = []			
	
	if "Search Results" in search_page.text:
		#Loop through all results. The results are in an ordered list!
		
		for raw_result in search_page.ol.find_all('li'):
			id_code = re.search(r'ID=([A-Za-z]\d{1,})', raw_result.a['href']).group(1)
			name = raw_result.a.text
			if instantiate == True:
				search_results.append(NistResult(id_code, name).instantiate())
			else:
				search_results.append(NistResult(id_code, name))
	#there is only one search result!
	else:
		#Find the ID (basically it's the ID-like text that appears the most on the page!)
		links_on_page = []
		#loop through all the links
		for link in search_page.find_all('a'):
			curr_href = link.attrs.get('href')
			#sometimes, there is no hfref (eg. it's an anchor or something), so check for NoneType
			if curr_href is not None and re.search(r'ID=([A-Za-z]\d{1,})', curr_href):
				#add the current link's ID code
				#group(0) is the entire search result, so group(1) is the bracketed thing in the search query
				links_on_page.append(re.search(r'ID=([A-Za-z]\d{1,})', curr_href).group(1))
		id_code = max(links_on_page, key=links_on_page.count)
		
		name = search_page.title.text
		
		if instantiate == True:
			search_results.append(NistResult(id_code, name, souped_page=search_page).instantiate())
		else:
			search_results.append(NistResult(id_code, name, souped_page=search_page))

	return search_results	

#TO DO
#nist plot title

#ERROR CHECKING
#creating nist compound by id (invalid and valid ids)
#getting properties (thermo, gas phase, henry law, gas ion energetics)  (existing and otheewise)
#getting spectra (existing and otherwise)
#listing properties
#Plotting: NIstplot - with valid and invalid objects; mixed objects (valid/invalid, NistCompoundById and other Nist objects, lists)
#changing options (title, style)
#Adding NistPlot objects to each other, plotting multiple things on same graph
#ALL NIST COMPOUNDS, looping through, searching using it
#Searching - for all things (cas, inchi, etc.) using or not using instantiate flag


#############################
#     SAMPLE USAGE CODE     #
#############################	
#x  = search_nist(formula='C6H6')
#print x
#print type(x)
#y = x.instantiate()
#print type(y), '\n'
#z = NistPlot(y, 'uv', ('r-',))
#z.show()
#x = NistCompoundById('C108883')
#y = NistCompoundById('C2809690')
#y = NistCompoundById('C1073672')
#z = NistCompoundById('C71432')
#p = NistPlot([x,y], 'ir', ('g-','b-'))
#q = NistPlot(y, 'ir', ('g-','b-'))
#p.show()
#q.show()
#p.show()

#p.show()
#x.get_spectra('mass')
#y = NistPlot(x, 'ir')
#y.show()
#x1 = NistPlot(x, 'uv', 'b-')
#y1 = NistPlot(y, 'uv')
#z = x1 + [y]
#print z.compounds
#z.show()

#x1.show()
#y1.show()

#x.list_available_spectra
#x.list_available_properties

#y = NistCompoundById('C71432')
#print 'done'
#x.get_spectra('mass')
#print 'done'
#a = NistPlot([x], 'mass', ('r-', 'b-', 'g-'))
#b = NistPlot([y], 'uv', ('b-', 'g-'))
#c = a + [y]
#c.show()
#a.show()

		#mass spec, UV and IR *use matplotlub to pklot this!
