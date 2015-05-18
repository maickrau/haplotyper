//./make_snps.exe Afile.fasta Bfile.fasta mutationProbability
//g++ make_snps.cpp fasta_utils.cpp -o make_snps.exe -std=c++11

#include <fstream>
#include <string>
#include <random>
#include <chrono>

#include "fasta_utils.h"

Genome mutate(Genome genome, double mutationProbability)
{
	char bases[4] = {'A', 'T', 'C', 'G'};
	Genome ret { genome };
    std::mt19937 mt(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> distMutation(0, 1);
    std::uniform_int_distribution<int> distBase(0, 3);
    for (auto& a : ret.bases)
	{
		if (distMutation(mt) < mutationProbability)
		{
			char b;
			do 
			{
				b = bases[distBase(mt)];
			} while (b == a);
			a = b;
		}
	}
	return ret;
}

int main(int argc, char** argv)
{
	Genome A = loadFasta(argv[1]);
	Genome B = mutate(A, std::stod(argv[3]));
	writeFasta(B, argv[2]);
}