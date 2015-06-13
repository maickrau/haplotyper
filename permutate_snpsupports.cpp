//g++ permutate_snpsupports.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o permutate_snpsupports.exe
//./permutate_snpsupports.exe inputSupportsFile outputSupportsFile renumberingFile

#include <random>
#include <chrono>
#include <algorithm>

#include "variant_utils.h"

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);
	size_t maxRead = supports[0].readNum;
	size_t maxSNP = supports[0].SNPnum;
	for (auto x : supports)
	{
		maxRead = std::max(maxRead, x.readNum);
		maxSNP = std::max(maxSNP, x.SNPnum);
	}
	std::vector<size_t> readPermutation;
	for (size_t i = 0; i < maxRead+1; i++)
	{
		readPermutation.push_back(i);
	}
	std::vector<size_t> SNPPermutation;
	for (size_t i = 0; i < maxSNP+1; i++)
	{
		SNPPermutation.push_back(i);
	}

    std::mt19937 mt(std::chrono::system_clock::now().time_since_epoch().count());
    std::shuffle(readPermutation.begin(), readPermutation.end(), mt);
    std::shuffle(SNPPermutation.begin(), SNPPermutation.end(), mt);

    SupportRenumbering renumbering;
    for (size_t i = 0; i < readPermutation.size(); i++)
    {
    	renumbering.addReadRenumbering(i, readPermutation[i]);
    }
    for (size_t i = 0; i < SNPPermutation.size(); i++)
    {
    	renumbering.addSNPRenumbering(i, SNPPermutation[i]);
    }

    assert(renumbering.checkValidity());
    supports = renumberSupports(supports, renumbering);
    writeSupports(supports, argv[2]);
    writeRenumbering(renumbering, argv[3]);
}