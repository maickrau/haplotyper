//g++ simulate_SNPfinding.cpp fasta_utils.cpp variant_utils.cpp -o simulate_SNPfinding.exe -std=c++11
//./simulate_SNPfinding.exe readsFile genomesFile outputFile readLength

#include <unordered_map>
#include <sstream>
#include <iostream>

#include "fasta_utils.h"
#include "variant_utils.h"

int main(int argc, char** argv)
{
	std::vector<Genome> reads = loadFastas(argv[1]);
	std::vector<Genome> genomes = loadFastas(argv[2]);
	size_t readLength = std::stoi(argv[4]);
	std::unordered_map<size_t, std::vector<size_t>> readsAtLocation;
	for (size_t i = 0; i < reads.size(); i++)
	{
		std::istringstream name { reads[i].name };
		std::string unused;
		size_t readLoc;
		name >> unused >> readLoc;
		readsAtLocation[readLoc].push_back(i);
	}
	std::vector<SNPSupport> supports;
	size_t SNPnum = 0;
	for (size_t i = 0; i < genomes[0].bases.size(); i++)
	{
		std::cout << i << "\n";
		bool hasSNP = false;
		for (size_t a = 1; a < genomes.size(); a++)
		{
			if (genomes[a].bases[i] != genomes[a-1].bases[i])
			{
				hasSNP = true;
				break;
			}
		}
		if (hasSNP)
		{
			size_t loc = i-(readLength-1);
			if (i < (readLength-1))
			{
				loc = 0;
			}
			for (; loc <= i; loc++)
			{
				for (auto x : readsAtLocation[loc])
				{
					supports.emplace_back(x, SNPnum, reads[x].bases[i-loc], 1);
				}
			}
			SNPnum++;
		}
	}
	writeSupports(supports, argv[3]);
}