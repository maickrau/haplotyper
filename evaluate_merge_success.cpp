//g++ evaluate_merge_success.cpp fasta_utils.cpp variant_utils.cpp -std=c++11 -o evaluate_merge_success.exe
//./evaluate_merge_success.exe allReadsFile mergeRenumbering renumbering1 renumbering2...

#include <iostream>
#include <map>
#include <string>
#include <sstream>

#include "fasta_utils.h"
#include "variant_utils.h"

std::vector<int> getGenomeIds(std::string fileName)
{
	std::vector<Genome> genomes = loadFastas(fileName);
	std::map<std::string, int> names;
	int nextId = 0;
	std::vector<int> ret;
	for (auto x : genomes)
	{
		std::istringstream s { x.name };
		std::string name;
		s >> name;
		if (names.find(name) == names.end())
		{
			names.insert({name, nextId});
			nextId++;
		}
		ret.push_back(names[name]);
	}
	return ret;
}

int main(int argc, char** argv)
{
	std::vector<int> genomeIds = getGenomeIds(argv[1]);
	SupportRenumbering mergeRenumbering = loadRenumbering(argv[2]);
	SupportRenumbering previousRenumberings = SupportRenumbering::identity(genomeIds.size(), 1);
	for (int i = 3; i < argc; i++)
	{
		SupportRenumbering add = loadRenumbering(argv[i]);
		previousRenumberings = previousRenumberings.merge(add);
	}
	std::vector<std::set<int>> renumberingsBeforeMerge;
	renumberingsBeforeMerge.resize(previousRenumberings.readSize());
	size_t totalAfterMerges = 0;
	for (size_t i = 0; i < genomeIds.size(); i++)
	{
		if (previousRenumberings.hasReadRenumbering(i))
		{
			totalAfterMerges++;
			size_t targetRead = previousRenumberings.getReadRenumbering(i);
			renumberingsBeforeMerge[targetRead].insert(genomeIds[i]);
		}
	}
	std::vector<std::set<int>> renumberings;
	std::vector<size_t> countAfterMerging;
	renumberings.resize(mergeRenumbering.readSize());
	countAfterMerging.resize(mergeRenumbering.readSize(), 0);
	for (size_t i = 0; i < mergeRenumbering.readSize(); i++)
	{
		countAfterMerging[mergeRenumbering.getReadRenumbering(i)]++;
	}
	size_t wrongBefore = 0;
	for (size_t i = 0; i < genomeIds.size(); i++)
	{
		if (previousRenumberings.hasReadRenumbering(i))
		{
			size_t previousTargetRead = previousRenumberings.getReadRenumbering(i);
			if (renumberingsBeforeMerge[previousTargetRead].size() == 1)
			{
				size_t targetRead = mergeRenumbering.getReadRenumbering(previousTargetRead);
				renumberings[targetRead].insert(genomeIds[i]);
			}
			else
			{
				wrongBefore++;
			}
		}
	}
	size_t total = 0;
	size_t wrong = 0;
	for (size_t i = 0; i < renumberings.size(); i++)
	{
		if (countAfterMerging[i] > 1)
		{
			total++;
		}
		if (renumberings[i].size() > 1)
		{
			wrong++;
		}
	}
	std::cout << mergeRenumbering.readSize() << "\n";
	std::cout << total << "\n";
	std::cout << wrong << "\n";
	std::cout << 1.0-(double)wrong/(double)total << "\n";
	std::cout << wrongBefore << "\n";
	std::cout << totalAfterMerges << "\n";
	std::cout << 1.0-(double)wrongBefore/(double)totalAfterMerges << "\n";
}