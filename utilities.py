import os
import shutil

class Gene:
	'''Class for holding information about one gene'''
	
	def __init__(self, name, chromosome, start, stop, cStart, cStop):
		self.name = name				#the ENSEMBL code for the gene
		self.chromosome = chromosome	#the chromosome the gene is found on
		self.start = start				#the location of the beginning of the gene
		self.stop = stop				#the location of the end of the gene
		self.cStart = cStart			#the start of the gene in the cDNA file (0 for non-coding)
		self.cStop = cStop				#the end of the gene in the cDNA file


class Species:
	'''Class for holding the genes for a specific species'''
	
	def __init__(self, name, folder, filter):
		print "reading data for " + name
		self.name = name		#name of the species
		self.folder = folder	#folder with genetic data
		self.genes = []			#list to hold genes
		self.matches = []		#list to hold matches for each gene
		self.matchvals = []		#list to hold values of matches
		self.filename = ""
		self.readCDNA(filter)
			
	
	def readCDNA(self, filter):
		'''fills in array of genes with data from the cDNA file'''
		
		cDNA = self.folder + "/" + self.name + "/cdna"
		files = os.listdir(cDNA)
		for file in files:
			if file[-11:] == "cdna.all.fa":
				self.filename = file[:-11]
		if self.filename == "":
			print "missing cDNA data file for " + self.name
		cDNA = cDNA + "/" + self.filename + "cdna.all.fa"
		
		cStart = 0
		cStop = 0
		file = open(cDNA, 'rb')
		loop = True
		save = False
		
		#a = 0
		
		while loop:
		
			#a = a+1
			#if a%100000 == 0:
			#	print a
		
			line = file.readline()
			if len(line) == 0 or line[len(line)-1] != '\n':		#end of file
				cStop = file.tell()
				if save:
					self.genes = self.genes + [Gene(name, chromosome, start, stop, cStart, cStop)]
				loop = False
				
			elif line[0] == ">":	#new gene -> save old gene and get new data
				if save:
					self.genes = self.genes + [Gene(name, chromosome, start, stop, cStart, cStop)]
				cStart = file.tell()
				
				#print line[0:len(line)-1]
				if "gene_biotype:protein_coding" in line and "transcript_biotype:protein_coding" in line and "chromosome:" in line:
					save = True
					pos = line.find("chromosome:") + 11
					if pos == 10:
						pos = line.find("supercontig:") + 12
						if pos == 11:
							pos = line.index("scaffold:") + 9
					pos = pos + line[pos:].index(":")
					[chromosome, pos] = self.getCDNAHeadVal(line, pos)	#find the chromosome number
					[start, pos] = self.getCDNAHeadVal(line, pos)		#find the start position
					start = int(start)
					[stop, pos] = self.getCDNAHeadVal(line, pos)			#find the end position
					stop = int(stop)
					pos = line.index("gene:") + 4
					name = self.getCDNAHeadVal(line, pos)[0]				#find the ENSEMBL gene code
					
					if filter != []:
						save = False
						for entry in filter:
							if name == entry:
								save = True
					
				else:
					save = False
				
			else:
				cStop = file.tell()
				
	def getCDNAHeadVal(self, line, pos):
		'''reads one value from the cDNA header line'''
		string = ""
		pos = pos + 1
		while line[pos] != ":" and line[pos] != " ":
			string = string + line[pos]
			pos = pos + 1
		return string, pos
	
	def sequence(self, gene):
		cDNA = self.folder + "/" + self.name + "/cdna/" + self.filename + "cdna.all.fa"
		file = open(cDNA, 'rb')
		sequence = ""
		start = self.genes[gene].cStart
		stop = self.genes[gene].cStop
		file.seek(start)
		while file.tell() < stop:
			line = file.readline()
			if line[len(line)-1] == '\n':
				line = line[:(len(line)-1)]
			sequence = sequence + line
		return sequence






def extractData(source, dest):
	'''extracts the cDNA data files from the rest of the data, and decompresses them'''
	
	names = os.listdir(source)
	for name in names:
		sourceFolder = source + "/" + name + "/cdna"
		destFolder = dest + "/" + name + "/cdna"
		if os.path.exists(sourceFolder):
			print "copying " + name
			files = os.listdir(sourceFolder)
			filename = ""
			for file in files:
				if file[-11:] == "cdna.all.fa":
					filename = file
				if file[-14:] == "cdna.all.fa.gz":
					filename = file
					os.system("sudo gunzip " + sourceFolder + "/" + filename)
					filename = filename[:-3]
			if filename != "":
				if not os.path.exists(destFolder):
					os.makedirs(destFolder)
				shutil.copyfile(sourceFolder + "/" + filename, destFolder + "/" + filename)
	print "done"
	return
	

def generateJobs(speciesList, data, machines, maxGenes):
	'''takes a list of species to be compared and splits it into jobs for several machines'''
	
	names = open(speciesList).read().split("\n")
	nSpecies = len(names)
	genes = []
	
	for i in range(0, nSpecies):
		name = names[i]
		while name[-1] == '\r':	#windows messes everything up
			name = name[:-1]
		names[i] = name
		
		animal = Species(name, data, [])
		if len(animal.genes) == 0:
			print "Error: no genes matched criteria for " + name
			names[i] = ""
		elif len(animal.genes) < maxGenes:
			genes = genes + [len(animal.genes)]
		else:
			genes = genes + [maxGenes]
	
	loop = True		#remove species that don't have any good genes
	while loop:
		loop = False
		for i in range(0, nSpecies):
			if names[i] == "":
				names = names[:i] + names[i+1:]
				nSpecies = nSpecies - 1
				loop = True
				break
	
	comp = [0]*machines
	jobs = []
	for i in range(0, machines):
		jobs = jobs + [[]]
	
	for i in range(0, nSpecies-1):
		for j in range(i+1, nSpecies):
			next = 0
			min = comp[0]
			for k in range(1, machines):
				if comp[k] < min:
					next = k
					min = comp[k]
			comp[next] = comp[next] + genes[i]*genes[j]
			jobs[next] = jobs[next] + [names[i] + "," + names[j]]
	if not os.path.exists("jobs"):
		os.makedirs("jobs")
	for i in range(0, machines):
		if jobs[i] == []:
			return
		file = open("jobs/job" + str(i) + ".txt", 'w')
		for line in jobs[i]:
			file.write(line + '\n')
		file.close()
	return
		
		
		















