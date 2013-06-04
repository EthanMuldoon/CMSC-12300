#include<string>
#include<iostream>
#include<fstream>
#include<dirent.h>
#include<vector>
#include<stdlib.h>
#include<thread>
#include<sys/stat.h>
using namespace std;

class Gene
{
private:
	string name;
	string chromosome;
	int start;
	int stop;
	int cStart;
	int cStop;
	
	friend class Species;
public:
	Gene(string name, string chromosome, int start, int stop, int cStart, int cStop)
	{
		this->name = name;
		this->chromosome = chromosome;
		this->start = start;
		this->stop = stop;
		this->cStart = cStart;
		this->cStop = cStop;
	}
};

class Species
{
private:
	string name;
	string folder;
	vector<Gene> genes;
	string filename;
	string cDNApath;
	int numGenes;
public:
	Species(string name, string folder, string filter[])
	{
		this->name = name;
		this->folder = folder;
		readCDNA(filter);
	}
	
	Species(string name, string folder)
	{
		this->name = name;
		this->folder = folder;
		string filter[] = {""};
		readCDNA(filter);
	}
	
	string getPath()
	{
		string cDNA = this->folder + "/" + this->name + "/cdna";
		DIR *dir;
		struct dirent *ent;
		
		if ((dir = opendir (cDNA.c_str())) != NULL) {
			while ((ent = readdir (dir)) != NULL) {
				string entry = ent->d_name;
				if(entry.length() > 11 && entry.substr(entry.length() - 11, entry.length()) == "cdna.all.fa")
				{
					this->filename = entry.substr(0, entry.length() - 11);
				}
				else if(entry.length() > 14 && entry.substr(entry.length() - 14, entry.length()) == "cdna.all.fa.gz")
				{
					this->filename = entry.substr(0, entry.length() - 14);
					string unzip = "sudo gunzip " + cDNA + "/" + entry;
					cout << "unzipping " << entry << endl;
					system(unzip.c_str());
				}
			}
			closedir (dir);
		}
		if(this->filename == "")
		{
			cerr << "missing cDNA data file for " << this->name << endl;
		}
		cDNA = cDNA + "/" + this->filename + "cdna.all.fa";
		this->cDNApath = cDNA;
		return cDNA;
	}
	
	void readCDNA(string filter[])
	{
		cout << "reading data for " << this->name << endl;
		string path = getPath();
		ifstream file(path);
		string line;
		string name;
		string chrom;
		int start = 0;
		int stop = 0;
		int cStart = 0;
		int cStop = 0;
		bool save = false;
		this->numGenes = 0;
		
		while(!file.eof())
		{
			getline(file, line);
			if(line[0] == '>')
			{
				if(save)
				{
					this->genes.push_back(Gene(name, chrom, start, stop, cStart, cStop));
					this->numGenes++;
				}
				cStart = file.tellg();
				// old condition: line.find("cdna:known") != string::npos
				if(line.find("gene_biotype:protein_coding") != string::npos 
					&& line.find("transcript_biotype:protein_coding") != string::npos
					&& line.find("chromosome:") != string::npos)
				{
					save = true;
					int n = line.find("chromosome:") + 11;
					n = line.find(":", n);
					chrom = getCDNAHeaderVal(line, n);
					string temp = getCDNAHeaderVal(line, n);
					start = atoi(temp.c_str());
					temp = getCDNAHeaderVal(line, n);
					stop = atoi(temp.c_str());
					n = line.find("gene:") + 4;
					name = getCDNAHeaderVal(line, n);
					
					if(filter[0] != "")
					{
						save = false;
						int i=0;
						while(filter[i] != "")
						{
							if(name == filter[i])
							{
								save = true;
							}
						}
					}
				}
				else
				{
					save = false;
				}
			}
			else
			{
				cStop = file.tellg();
			}
		}
		file.close();
		
		return;
	}
	
	string getCDNAHeaderVal(string line, int &pos)
	{
		pos++;
		int start = pos;
		while(line[pos] != ':' && line[pos] != ' ')
		{
			pos++;
		}
		string entry = line.substr(start, pos-start);
		return entry;
	}
	
	int getNumGenes()
	{
		return numGenes;
	}
	
	string getName()
	{
		return name;
	}
	
	string sequence(int n)
	{
		ifstream file(cDNApath);
		file.seekg(genes[n].cStart);
		string sequence = "";
		string line;
		while(file.tellg() < genes[n].cStop)
		{
			getline(file, line);
			sequence = sequence + line;
		}
		return sequence;
	}
	
	string getGeneName(int n)
	{
		return genes[n].name;
	}
};

struct GenePair
{
	string geneA;
	string geneB;
	int indexA;
	int indexB;
};

struct Match
{
	int match;
	double diff;
	Match()
	{
		this->match = -1;
		this->diff = -1;
	}
};

const string BASES = "ATCGWSMKRYBDHVN";
static bool BASE_MATCHES[60] = {1,0,0,0,1,0,1,0,1,0,0,1,1,1,1,0,1,0,0,1,0,0,1,0,1,1,1,1,0,1,0,0,1,0,0,1,1,0,0,1,1,0,1,1,1,0,0,0,1,0,1,0,1,1,0,1,1,0,1,1};

bool baseMatch(char baseA, char baseB)
{
	if(baseA == baseB)					//Check if there is an exact match
	{
		return true;
	}
	if(baseA == 'N' || baseB == 'N')	//'N' stands for any base
	{
		return true;
	}
	/*	The data doesn't seem to include anything but ATCG and N
		I am leaving this part out so the program goes faster
		
	int indA = BASES.find(baseA);		//Look up the bases in the matching matrix
	int indB = BASES.find(baseB);
	if(indA < 4)
	{
		return BASE_MATCHES[indA*15+indB];
	}
	if(indB < 4)
	{
		return BASE_MATCHES[indB*15+indA];
	}
	*/
	return false;						//Otherwise, assume the bases don't match
}

int geneCompare(string geneA, string geneB)
{
	if(geneA == geneB)
	{
		return 0;
	}
	int N = geneA.length();
	int M = geneB.length();
	int MAX = max(N, M);
	int DELTA = N-M;
	int fore[4*MAX + 1];
	int back[4*MAX + 1];
	for(int i = 0; i < 4*MAX + 1; i++)
	{
		fore[i] = 0;
		back[i] = N;
	}
	int* forward = &fore[2*MAX];
	int* backward = &back[2*MAX];
	int next = 0;
	int x, y;
	for(int D = 0; D < MAX; D++)
	{
		for(int k = -D; k < D+1; k++)
		{
			if(k == -D || (k != -D && forward[k-1] < forward[k+1] && forward[k] < forward[k+1]))
			{
				x = forward[k+1];		//move down
			}
			else if(k == D || (k != D && forward[k] <= forward[k-1]))
			{
				x = forward[k-1] + 1;	//move right
			}
			else
			{
				x = forward[k] + 1;		//move diagonally
			}
			y = x-k;
			while(x < N && y < M && baseMatch(geneA[x], geneB[y]))
			{
				x++;
				y++;
			}
			if(k != -D)
			{
				forward[k-1] = next;
			}
			next = x;
			if(x >= backward[k] && k > DELTA-D && k < DELTA+D)
			{
				return 2*D - 1;
			}
		}
		forward[D] = next;
		
		for(int k = DELTA+D; k >= DELTA-D; k--)
		{
			if(k == DELTA+D || (k != DELTA+D && backward[k-1] < backward[k+1] && backward[k-1] < backward[k]))
			{
				x = backward[k-1];	//move up
			}
			else if(k == DELTA-D || (k != DELTA-D && backward[k+1] <= backward[k]))
			{
				x = backward[k+1] - 1;		//move left
			}
			else
			{
				x = backward[k] - 1;		//move diagonally
			}
			y = x-k;
			while(x > 0 && y > 0 && baseMatch(geneA[x-1], geneB[y-1]))
			{
				x--;
				y--;
			}
			if(k != DELTA+D)
			{
				backward[k+1] = next;
			}
			next = x;
			if(x <= forward[k] && k >= -D && k <= D)
			{
				return 2*D;
			}
		}
	}
	return MAX;
}

void saveOutput(Species A, Species B, int geneLimit, Match *matchA)
{
	int genes = A.getNumGenes();
	
	if(geneLimit != 0 && genes > geneLimit)
	{
		genes = geneLimit;
	}
	
	umask(0);
	string path = "output";
	mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	path = "output/"+A.getName();
	mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	ofstream output(path + "/" + B.getName() + ".txt");
	
	output << A.getName() << "," << B.getName() << endl;
	output << A.getNumGenes() << "," << B.getNumGenes() << "," << geneLimit << endl;
	for(int n=0; n < genes; n++)
	{
		output << A.getGeneName(n) << "," << B.getGeneName(matchA[n].match) << "," << matchA[n].diff << endl;
	}
	output.close();
	return;
}

/* old single-threaded version
int speciesCompare(Species A, Species B)
{
	if(B.getNumGenes() > A.getNumGenes())
	{
		Species C = A;
		A = B;
		B = C;
	}
	int genesA = A.getNumGenes();
	int genesB = B.getNumGenes();
	bool checked[genesA];
	int differences[genesA];
	int matches[genesA];
	for(int n = 0; n < genesA; n++)
	{
		checked[n] = false;
		differences[n] = 0;
		matches[n] = -1;
	}
	
	for(int n = 0; n < genesA; n++)
	{
		if(!checked[n])
		{
			int minDiff = -1;
			int minMatch = 0;
			for(int l = n; l < genesA; l++)
			{
				if(A.getGeneName(n) == A.getGeneName(l))
				{
					checked[l] = true;
					for(int m = 0; m < genesB; m++)
					{
						int diff = geneCompare(A.sequence(n), B.sequence(m));
						if(minDiff == -1 || diff < minDiff)
						{
							minDiff = diff;
							minMatch = m;
						}
					}
				}
			}
			differences[n] = minDiff;
			matches[n] = minMatch;
		}
	}
	int totalDiff = 0;
	for(int n = 0; n < genesA; n++)
	{
		totalDiff = totalDiff + differences[n];
	}
	return totalDiff;
}
*/



static int QUEUE_LENGTH = 64;
static int N_THREADS = 32;
bool done = false;
int reading = 0;
int writing = 0;
int entries = 0;
mutex mQueue;
mutex mOutput;
condition_variable empty;
condition_variable full;

void scConsumer(GenePair *scQueue, Match *matchA, Match *matchB)
{
	GenePair genes;
	int diff;
	while(true)
	{
		{
			unique_lock<mutex> lock(mQueue);
			while(entries == 0)
			{
				if(done)
				{
					return;
				}
				empty.wait(lock);
			}
			genes = scQueue[reading];
			
			//cout << "r - " << reading;
			//if(reading < 10){ cout << " "; }
			//cout << " - " << entries << endl;	//debug
			
			reading = (reading + 1) % QUEUE_LENGTH;
			entries--;
			if(entries == QUEUE_LENGTH - 1)
			{
				full.notify_one();
			}
		}
		
		diff = geneCompare(genes.geneA, genes.geneB);
		
		{
			unique_lock<mutex> lock(mOutput);
			if(matchA[genes.indexA].diff == -1 || diff < matchA[genes.indexA].diff)
			{
				matchA[genes.indexA].match = genes.indexB;
				matchA[genes.indexA].diff = diff;
			}
			if(matchB[genes.indexB].diff == -1 || diff < matchB[genes.indexB].diff)
			{
				matchB[genes.indexB].match = genes.indexA;
				matchB[genes.indexB].diff = diff;
			}
		}
	}
}

void scProducer(Species A, Species B, GenePair *scQueue, int geneLimit)
{
	GenePair next;
	int genesA = A.getNumGenes();
	int genesB = B.getNumGenes();
	
	if(geneLimit != 0)
	{
		if(genesA > geneLimit)
		{
			genesA = geneLimit;
		}
		if(genesB > geneLimit)
		{
			genesB = geneLimit;
		}
	}
	
	int l = 0;
	
	for(int n = 0; n < genesA; n++)
	{
		if(n/(double) genesA >= l/40.0)	//print a progress bar
		{
			if(l%4 == 0)
			{
				cout << l*2.5 << flush;
			}
			else
			{
				cout << '.' << flush;
			}
			l++;
		}
		
		next.indexA = n;
		next.geneA = A.sequence(n);
		for(int m = 0; m < genesB; m++)
		{
			next.indexB = m;
			next.geneB = B.sequence(m);
			
			{
				unique_lock<mutex> lock(mQueue);
				while(entries == QUEUE_LENGTH)
				{
					full.wait(lock);
				}
				scQueue[writing] = next;
				
				//cout << "w - " << writing;
				//if(writing < 10){ cout << " "; }
				//cout << " - " << entries << endl;	//debug
				
				writing = (writing + 1) % QUEUE_LENGTH;
				entries++;
				empty.notify_one();
			}
		}
	}
	mQueue.lock();
	empty.notify_all();
	mQueue.unlock();
	cout << "100" << endl;
	
	return;
}

void speciesCompare(Species A, Species B, int geneLimit)
{
	cout << "comparing " << A.getName() << " to " << B.getName() << endl;
	int genesA = A.getNumGenes();
	int genesB = B.getNumGenes();
	
	if(geneLimit != 0)
	{
		if(genesA > geneLimit)
		{
			genesA = geneLimit;
		}
		if(genesB > geneLimit)
		{
			genesB = geneLimit;
		}
	}
	
	GenePair *scQueue = new GenePair[QUEUE_LENGTH];
	Match *matchA = new Match[genesA];
	Match *matchB = new Match[genesB];
	
	done = false;
	
	thread reader(scProducer, A, B, scQueue, geneLimit);
	thread *workers = new thread[N_THREADS];
	for(int i=0; i < N_THREADS; i++)
	{
		workers[i] = thread(scConsumer, scQueue, matchA, matchB);
	}

	reader.join();
	done = true;
	for(int i=0; i < N_THREADS; i++)
	{
		workers[i].join();
	}
	
	for(int n = 0; n < genesA; n++)
	{
		matchA[n].diff = matchA[n].diff*100/A.sequence(n).length();
	}
	for(int n = 0; n < genesB; n++)
	{
		matchB[n].diff = matchB[n].diff*100/B.sequence(n).length();
	}
	saveOutput(A, B, geneLimit, matchA);
	saveOutput(B, A, geneLimit, matchB);
	
	return;
}





int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cerr << "usage:  'a.out data_folder job_file gene_limit'" << endl;
		return 0;
	}
	
	string folder (argv[1]);
	string job (argv[2]);
	int geneLimit = atoi(argv[3]);
	
	ifstream file(job);
	string line;
	vector<string> namesA;
	vector<string> namesB;
	
	while(!file.eof())
	{
		getline(file, line);
		int n = line.find(",");
		if(n == string::npos)
		{
			break;
		}
		
		//cout << line.substr(0, n) << " - " << line.substr(n+1, line.length()-n-1) << endl;
		
		string B = line.substr(n+1, string::npos);
		
		while(B[B.length()-1] == '\r')		//deal with possible windows line endings
		{
			B = B.substr(0, B.length()-1);
		}
		namesA.push_back(line.substr(0, n));
		namesB.push_back(B);
	}
	
	for(int n=0; n<namesA.size(); n++)
	{
		Species A = Species(namesA[n], folder);
		Species B = Species(namesB[n], folder);
		
		speciesCompare(A, B, geneLimit);
	}
	
	return 0;
}