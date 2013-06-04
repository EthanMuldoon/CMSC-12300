import random
import pylab
import os
import shutil
import heapq

BASES = "ATCGWSMKRYBDHVN"
#                     A T C G W S M K R Y B D H V N
BASE_MATCH_MATRIX = [[1,0,0,0,1,0,1,0,1,0,0,1,1,1,1], #A	matrix used for baseMatch function
					 [0,1,0,0,1,0,0,1,0,1,1,1,1,0,1], #T
					 [0,0,1,0,0,1,1,0,0,1,1,0,1,1,1], #C
					 [0,0,0,1,0,1,0,1,1,0,1,1,0,1,1]] #G

class Node:
	'''Class for holding the genetic tree data'''
	
	def __init__(self, mutations, children, name):
		self.mutations = mutations
		self.children = children
		self.name = name


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
				if "cdna:known" in line and "gene_biotype:protein_coding" in line and "transcript_biotype:protein_coding" in line:
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


def run(folder):
	'''runs the program'''
	
	animals = os.listdir(folder)
	
	print "comparing animals"
	
	diffMatrix = makeDiffMatrix(animals, folder)
	tree = generateTree(diffMatrix, animals)
	drawTree(tree, len(tree)-1)
	
	return diffMatrix
	



def loadData(names, folder):
	animals = []
	for name in names:
		animals = animals + [Species(name, folder)]

def getFilter(filterFile):
	'''load a filter file to reduce the number of genes being compared'''
	file = open(filterFile)
	filter = []
	for line in file:
		if line[len(line)-1] == '\n':
			line = line[:len(line)-1]
		filter = filter + [line]
	return filter

def makeDiffMatrix(animals, folder):
	totalSpecies = len(animals)
	diffMatrix = []
	for n in range(0, totalSpecies):
		diffMatrix = diffMatrix + [[-1]*totalSpecies]
	for n in range(0, totalSpecies-1):
		for m in range(n+1, totalSpecies):
			diffMatrix[n][m] = min(findAvgDiff(folder + "/" + animals[n] + "/" + animals[m] + ".txt"), findAvgDiff(folder + "/" + animals[m] + "/" + animals[n] + ".txt"))
			
			#print animals[n] + ", " + animals[m] + ":  " + str(findAvgDiff(folder + "/" + animals[n] + "/" + animals[m] + ".txt"))
			#print animals[m] + ", " + animals[n] + ":  " + str(findAvgDiff(folder + "/" + animals[m] + "/" + animals[n] + ".txt"))
			#print "  avg = " + str(diffMatrix[n][m])
			
	return diffMatrix


def speciesCompare(A, B):
	'''finds the number of mutations between two species'''
	
	checked = [False] * len(A.genes)
	differences = []
	matches = []		#eventually I would like to report the k most different genes
						# and what genes they were matched to
	for n in range(0, len(A.genes)):
		if not checked[n]:
			checked[n] = True
			minDiff = -1
			minMatch = 0
			for m in range(0, len(B.genes)):
				diff = geneCompare(A.sequence(n), B.sequence(m))
				if minDiff == -1 or diff < minDiff:
					minDiff = diff
					minMatch = m
			for l in range(n+1, len(A.genes)):
				if A.genes[n].name == A.genes[l].name:
					checked[l] = True
					for m in range(0, len(B.genes)):
						diff = geneCompare(A.sequence(l), B.sequence(m))
						if minDiff == -1 or diff < minDiff:
							minDiff = diff
							minMatch = m
			differences = differences + [minDiff]
			matches = matches + [minMatch]
		else:
			differences = differences + [0]
			matches = matches + [-1]
	totalDiff = sum(differences)
	return totalDiff

def baseMatch(base1, base2):
	'''checks if bases are equal, considering ambiguous notation'''
	#see:  http://en.wikipedia.org/wiki/Nucleic_acid_notation
	
	if base1 == base2:
		return True
	if base1 == "N" or base2 == "N":
		return True
	base1 = BASES.index(base1)
	base2 = BASES.index(base2)
	if base1 < 4:
		return bool(BASE_MATCH_MATRIX[base1][base2])
	if base2 < 4:
		return bool(BASE_MATCH_MATRIX[base2][base1])
	return False

def geneCompare(geneA, geneB):
	'''Compares the differences in a single gene'''
	#see:  text comparison.pdf
	
	N = len(geneA)
	M = len(geneB)
	MAX = max(0, N, M)				#if the strings are completely different, D = MAX
	
	positions = [0]*(2*MAX + 1)		#keeps track of the current position along each diagonal
	new = 0							#temporarily stores new position along one diagonal
	
	for D in range(0, MAX):			#for each step, calculate the new positions along each diagonal
		for k in range(-D, D+1):	#after D steps, we can only reach the D'th diagonal
			#for each diagonal k, chose whether to attach to top, left, or diagonal chain
			if k == -D or (k != -D and positions[k-1] < positions[k+1] and positions[k] < positions[k+1]):
				x = positions[k+1]	
			elif k == D or (k != D and positions[k] <= positions[k-1]):
				x = positions[k-1] + 1
			else:
				x = positions[k] + 1
			y = x-k
			
			#extend chain by as many null edges as possible
			while x < N and y < M and baseMatch(geneA[x], geneB[y]):	#previously geneA[x] == geneB[y]
				x = x+1
				y = y+1
			
			#each iteration of the loop requires the original positions[k-1] value, so we have to 
			#wait until the next iteration to modify the value
			if k != -D:
				positions[k-1] = new
			new = x
			
			#if we have reached the bottom right corner then we are done, return the differences
			if x >= N and y >= M:
				#print str(D) + " differences"
				return D
		positions[D] = new
	return MAX


def generateTree(diffMatrix, names):
	'''generates an evolutionary tree based on the differences between species'''
	
	animals = len(names)
	tree = []
	indices = range(0, animals)
	
	#print indices
	
	for name in names:
		tree = tree + [Node(0, [], name)]
	while animals > 1:
		minRow = 0
		minCol = 1
		minVal = diffMatrix[0][1]
		for row in range(0, animals-1):
			for col in range(row+1, animals):
				if diffMatrix[row][col] < minVal:
					minVal = diffMatrix[row][col]
					minRow = row
					minCol = col
		
		print str(minRow) + "," + str(minCol) + " = " + str(minVal)
		
		#mutations = minVal + (tree[indices[minRow]].mutations + tree[indices[minCol]].mutations)/2
		mutations = minVal
		
		tree = tree + [Node(mutations, [indices[minRow], indices[minCol]], "")]
		indices[minRow] = len(tree)-1
		indices[minCol:] = indices[(minCol+1):]
		
		#print indices
		
		animals = animals - 1
		
		for row in range(0, animals-1):
			for col in range(row+1, animals):	#combine two rows and columns into one and shrink the matrix
				if col == minRow:
					diffMatrix[row][col] = min(diffMatrix[row][col], diffMatrix[row][minCol])
				elif col < minCol and row == minRow:
					diffMatrix[row][col] = min(diffMatrix[row][col], diffMatrix[col][minCol])
				elif row == minRow:
					diffMatrix[row][col] = min(diffMatrix[row][col+1], diffMatrix[minCol][col+1])
				elif col >= minCol and row < minCol:
					diffMatrix[row][col] = diffMatrix[row][col+1]
				elif col > minCol and row >= minCol:
					diffMatrix[row][col] = diffMatrix[row+1][col+1]
	return tree

def drawTree(tree, start):
	'''draws an evolutionary tree'''
	
	maxMutations = tree[start].mutations
	current = [start]
	leafs = []
	while current != []:
		children = tree[current[0]].children
		if children == []:
			leafs = leafs + [current[0]]
			tree[current[0]].visited = True
		else:
			current = current + children
			tree[current[0]].visited = False
		current = current[1:]
	
	lines = range(1, len(leafs)+1)
	for n in range(0, len(lines)):
		lines[n] = lines[n]*1.0
		tree[leafs[n]].line = n
	
	position = 0
	
	pylab.figure()
	while position < maxMutations:
		minVal = maxMutations
		minPos = start
		for n in range(0, len(tree)):
			if tree[n].mutations < minVal and (not tree[n].visited):
				minPos = n
				minVal = tree[n].mutations
		
		for n in range(0, len(lines)):
			if lines[n] != 0:
				pylab.plot([maxMutations-minVal, maxMutations-position], [lines[n]]*2, 'k-')
		position = minVal
		
		#print position
		#print minPos
		#print tree[minPos].children
		
		a = tree[tree[minPos].children[0]].line
		b = tree[tree[minPos].children[1]].line
		
		#print a
		#print b
		#print lines[a]
		#print lines[b]
		
		tree[minPos].line = a
		tree[minPos].visited = True
		pylab.plot([maxMutations-position]*2, [lines[a], lines[b]], 'k-')
		lines[a] = (lines[a] + lines[b])/2
		lines[b] = 0
	
	axis = [-0.1, 1.3*maxMutations, 0.9, len(leafs)+0.1]
	pylab.axis(axis)
	pylab.axis('off')
	x = 1.01*maxMutations
	for y in range(0, len(leafs)):
		pylab.text(x, y+1, tree[leafs[y]].name)
	pylab.show()



def findAvgDiff(file):
	'''calculate the average differences from a species comparison file'''
	
	file = open(file).read().split("\n")
	if file[-1] == "":
		file = file[:-1]
	for n in range(0, len(file)):
		file[n] = file[n].split(",")
		if n > 0:
			file[n][2] = float(file[n][2])
	
	file[1][0] = float(file[1][0])
	file[1][1] = float(file[1][1])
	if file[1][1] > file[1][2]:
		ratio = (file[1][2]/file[1][1])/4	#the 4 is there to ensure we have eliminated most of the false matches
	else:
		ratio = 1
	nCorrect = int(ratio * file[1][2])
	
	sorted = []
	for i in range(2, nCorrect+2):
		heapq.heappush(sorted, -1*file[i][2])
	for i in range(nCorrect+2, len(file)):
		heapq.heappushpop(sorted, -1*file[i][2])
	
	sum = 0
	for i in range(0, nCorrect):
		sum = sum + heapq.heappop(sorted)
	avg = -1*sum/nCorrect
	
	return avg









def test(length, maxMutations, steps):
	'''debug - tests gene comparison'''
	
	parents = [0]
	children = []
	sequence = ""
	for n in range(0, length):
		sequence = sequence + BASES[random.randint(0,3)]
	#sequence = BASES[0:4]*(length/4)+BASES[0:(length%4)]
	
	tree = [Node(0, [], sequence)]
	for n in range(0, steps):
		for m in range(0, len(parents)):
			newChildren = []
			for l in range(0,2):
				sequence = tree[parents[m]].name
				mutations = random.randint(1,maxMutations)
				for k in range(0,mutations):
					type = random.randint(0,2)
					if type == 0: 		#substitution
						position = random.randint(0, len(sequence)-1)
						old = sequence[position]
						new = BASES[0:4]
						old = new.index(old)
						new = new[random.randint(0,2)]
						sequence = sequence[0:position] + new + sequence[(position+1):len(sequence)]
					elif type == 1 or (len(sequence)<(length/2)): 	#addition
						position = random.randint(0, len(sequence))
						new = BASES[random.randint(0,3)]
						sequence = sequence[0:position] + new + sequence[position:len(sequence)]
					else: 				#deletion
						position = random.randint(0, len(sequence)-1)
						sequence = sequence[0:position] + sequence[(position+1):len(sequence)]
				newChildren = newChildren + [len(tree)]
				tree = tree + [Node(mutations, [], sequence)]
			tree[parents[m]].children = [newChildren]
			
			#print newChildren
			
			children = children + newChildren
		parents = children
		children = []
	totalSpecies = len(parents)
	
	#for node in tree:
	#	print node.name
	
	sequences = []
	diffMatrix = []
	names = []
	for n in range(0, totalSpecies):
		sequences = sequences + [tree[parents[n]].name]
		diffMatrix = diffMatrix + [[-1]*totalSpecies]
		names = names + ["test " + str(n)]
		tree[parents[n]].name = names[n]
	
	for n in range(0, totalSpecies-1):
		for m in range(n+1, totalSpecies):
			diffMatrix[n][m] = geneCompare(sequences[n], sequences[m])
			print str(n) + "," + str(m) + " = " + str(diffMatrix[n][m])
			
	#print diffMatrix
	print ""
	
	tree2 = generateTree(diffMatrix, names)
	
	drawTree(tree2, len(tree2)-1)
	














