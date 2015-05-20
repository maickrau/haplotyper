#include <algorithm>

#include "variant_utils.h"

SNPPosition::SNPPosition(Genome gen1, Genome gen2) : leftFlank(), rightFlank(), variants()
{
	size_t loc = 0;
	while (gen1.bases[loc] == gen2.bases[loc])
	{
		loc++;
	}
	leftFlank = gen1.bases.substr(0, loc);
	rightFlank = gen1.bases.substr(loc+1);
	variants.insert(gen1.bases[loc]);
	variants.insert(gen2.bases[loc]);
	std::reverse(leftFlank.begin(), leftFlank.end());
}

SNPSupport::SNPSupport(size_t readNum, size_t SNPnum, char variant) : readNum(readNum), SNPnum(SNPnum), support(1), variant(variant) 
{
}

SNPSupport::SNPSupport(size_t readNum, size_t SNPnum, char variant, double support) : readNum(readNum), SNPnum(SNPnum), support(support), variant(variant) 
{
}

std::vector<SNPSupport> loadSupports(std::string fileName)
{
	std::vector<SNPSupport> result;
	std::ifstream file { fileName };
	do
	{
		size_t readNum, SNPnum;
		char variant;
		double support;
		file >> readNum >> SNPnum >> variant >> support;
		if (file.good())
		{
			result.emplace_back(readNum, SNPnum, variant, support);
		}
	} while (file.good());
	return result;
}

void writeSupports(std::vector<SNPSupport> supports, std::string fileName)
{
	std::ofstream file { fileName };
	for (auto x : supports)
	{
		file << x.readNum << " " << x.SNPnum << " " << x.variant << " " << x.support << "\n";
	}
}
