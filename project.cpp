#include<string>
#include<iostream>
#include<dirent.h>
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
	
	string getName()
	{
		return name;
	}
};

class Species
{
private:
	string name;
	string folder;
	Gene genes[];
	string filename;
	int numGenes;
public:
	Species(string name, string folder, string filter[])
	{
		this->name = name;
		this->folder = folder;
		readCDNA(filter);
	}
	
	string getFileName()
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
			}
			closedir (dir);
		}
		if(this->filename == "")
		{
			cout << "missing cDNA data file for " << this->name << endl;
		}
		cDNA = cDNA + "/" + this->filename + "cdna.all.fa";
		return cDNA;
	}
	
	void readCDNA(string filter[])
	{
		cout << "reading data for " << this->name << endl;
		string cDNA = getFileName();
		
	}
	
	int getNumGenes()
	{
		return numGenes;
	}
	
	string sequence(int n)
	{
		
	}
	
	string getGeneName(int n)
	{
		return genes[n].getName();
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
		//cout << "exact match" << endl;
		return 0;
	}
	
	int N = geneA.length();
	int M = geneB.length();
	int MAX = max(N, M);
	int DELTA = N-M;
	int fore[2*MAX + 1];
	int back[2*MAX + 1];
	for(int i = 0; i < 2*MAX + 1; i++)
	{
		fore[i] = 0;
		back[i] = N;
	}
	int* forward = &fore[MAX];
	int* backward = &back[MAX];
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
		backward[DELTA-D] = next;
	}
	return MAX;
}

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
			/*	probably redundant
			checked[n] = true;
			for(int m = 0; m < genesB; m++)
			{
				int diff = geneCompare(A.sequence(n), B.sequence(m));
				if(minDiff == -1 || diff < minDiff)
				{
					minDiff = diff;
					minMatch = m;
				}
			}
			*/
			for(int l = n; l < genesB; l++)
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

int main(int argc, char** argv)
{
	string a, b;
	while(true)
	{
		cin >> a;
		cin >> b;
		cout << geneCompare(a, b) << " differences" << endl << endl;
	}
	
	
	return 0;
}