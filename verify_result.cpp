//g++ verify_result.cpp fasta_utils.cpp -std=c++11 -o verify_result.exe
//./verify_result.exe resultsFile allReadsFile

#include <iostream>
#include <algorithm>
#include <vector>
#include <map>

#include "fasta_utils.h"

std::vector<size_t> loadAssignments(std::string fileName)
{
	std::ifstream resultsFile { fileName };
	std::vector<size_t> assignments;
	while (resultsFile.good())
	{
		size_t assignment;
		resultsFile >> assignment;
		if (resultsFile.good())
		{
			assignments.push_back(assignment);
		}
	}
	assignments.pop_back();
	return assignments;
}

std::vector<size_t> getGenomeIds(std::vector<Genome> genomes)
{
	std::vector<size_t> result;
	std::map<char, size_t> assigned;
	for (auto x : genomes)
	{
		auto found = assigned.find(x.name[0]);
		if (found == assigned.end())
		{
			assigned.emplace(x.name[0], assigned.size());
		}
		result.push_back(assigned[x.name[0]]);
	}
	return result;
}

std::vector<std::vector<size_t>> getAllPermutations(std::vector<size_t> start, std::vector<size_t> itemsLeft)
{
	if (itemsLeft.empty())
	{
		return std::vector<std::vector<size_t>> { start };
	}
	std::vector<std::vector<size_t>> ret;
	for (size_t i = 0; i < itemsLeft.size(); i++)
	{
		std::vector<size_t> nextStart = start;
		nextStart.push_back(itemsLeft[i]);
		std::vector<size_t> nextItemsLeft = itemsLeft;
		nextItemsLeft.erase(nextItemsLeft.begin()+i);

		std::vector<std::vector<size_t>> parts = getAllPermutations(nextStart, nextItemsLeft);
		ret.insert(ret.end(), parts.begin(), parts.end());
	}
	return ret;
}

std::vector<std::vector<size_t>> getAllPermutations(size_t size)
{
	std::vector<size_t> startItems;
	for (size_t i = 0; i < size; i++)
	{
		startItems.push_back(i);
	}
	return getAllPermutations({}, startItems);
}

size_t calculateResult(std::vector<size_t> genomes, std::vector<size_t> assignments)
{
	std::vector<std::vector<size_t>> permutations = getAllPermutations(*max_element(genomes.begin(), genomes.end())+1);
	size_t minValue = -1;
	for (auto permutation : permutations)
	{
		size_t thisValue = 0;
		for (size_t i = 0; i < genomes.size(); i++)
		{
			if (genomes[i] != permutation[assignments[i]])
			{
				thisValue++;
			}
		}
		minValue = std::min(minValue, thisValue);
	}
	return minValue;
}

int main(int argc, char** argv)
{
	std::vector<size_t> genomes = getGenomeIds(loadFastas(argv[2]));
	std::vector<size_t> assignments = loadAssignments(argv[1]);
	size_t result = calculateResult(genomes, assignments);
	std::cout << result;
}