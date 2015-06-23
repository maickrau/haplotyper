//g++ mutate_snpsupports.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o mutate_snpsupports.exe
//./mutate_snpsupports.exe inputSupportsFile outputSupportsFile mutationProbability

#include <random>
#include <chrono>

#include "variant_utils.h"

std::vector<SNPSupport> mutate(std::vector<SNPSupport> supports, double mutationProbability)
{
    std::mt19937 mt(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> distMutation(0, 1);
    std::uniform_int_distribution<int> distBase(0, 3);
	for (auto& x : supports)
	{
		if (distMutation(mt) <= mutationProbability)
		{
			char newVariant;
			do
			{
				newVariant = "ACTG"[distBase(mt)];
			} while (newVariant == x.variant);
			x.variant = newVariant;
		}
	}
	return supports;
}

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);
	supports = mutate(supports, std::stod(argv[3]));
	writeSupports(supports, argv[2]);
}
