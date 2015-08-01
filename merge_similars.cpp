//g++ merge_similars.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o merge_similars.exe
//./merge_similars.exe inputSupportsFile coverageLimit outputSupportsFile

#include <unordered_map>
#include <algorithm>
#include <iostream>

#include "variant_utils.h"

size_t findSNPWithHighestNecessaryCoverage(std::vector<SNPSupport> supports, size_t minCoverage)
{
	size_t maxSNP = 0;
	for (auto x : supports)
	{
		maxSNP = std::max(maxSNP, x.SNPnum);
	}
	maxSNP++;
	std::vector<size_t> coverages;
	coverages.resize(maxSNP, 0);
	for (auto x : supports)
	{
		coverages[x.SNPnum]++;
	}
	size_t highestIndex = 0;
	for (size_t i = 1; i < coverages.size(); i++)
	{
		if (coverages[i] > coverages[highestIndex])
		{
			highestIndex = i;
		}
	}
	if (coverages[highestIndex] >= minCoverage)
	{
		return highestIndex;
	}
	return -1;
}

std::vector<SNPLine> filterLines(std::vector<SNPLine> lines, size_t SNPposition)
{
	std::vector<SNPLine> ret;
	for (auto x : lines)
	{
		if (std::any_of(x.variantsAtLocations.begin(), x.variantsAtLocations.end(), [SNPposition](std::pair<size_t, char> v) { return v.first == SNPposition; }))
		{
			ret.push_back(x);
		}
	}
	return ret;
}

size_t lineDifference(SNPLine left, SNPLine right)
{
	std::set<size_t> leftSNPs;
	for (auto x : left.variantsAtLocations)
	{
		leftSNPs.insert(x.first);
	}
	size_t result = leftSNPs.size();
	for (auto x : right.variantsAtLocations)
	{
		if (leftSNPs.count(x.first) > 0)
		{
			if (left.variantAt(x.first) != x.second)
			{
				result += leftSNPs.size();
			}
			else
			{
				result--;
			}
		}
		else
		{
			result++;
		}
	}
	return result;

}

std::pair<size_t, size_t> findMostSimilarRows(std::vector<SNPSupport> supports, size_t SNPposition)
{
	std::vector<SNPLine> lines = makeLines(supports);
	lines = filterLines(lines, SNPposition);
	size_t bestLeft = 0;
	size_t bestRight = 0;
	size_t bestDifference = -1;
	for (size_t i = 0; i < lines.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			size_t difference = lineDifference(lines[i], lines[j]);
			if (difference < bestDifference)
			{
				bestDifference = difference;
				bestLeft = lines[i].readNum;
				bestRight = lines[j].readNum;
			}
		}
	}
	return std::pair<size_t, size_t> { bestLeft, bestRight };
}

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);
	size_t coverageLimit = std::stoi(argv[2]);
	size_t mergedRows = 0;
	while (true)
	{
		size_t SNPposition = findSNPWithHighestNecessaryCoverage(supports, coverageLimit);
		if (SNPposition == -1)
		{
			break;
		}
		std::pair<size_t, size_t> similars = findMostSimilarRows(supports, SNPposition);
		supports = mergeRows(supports, similars.first, similars.second);
		mergedRows++;
	}
	std::cerr << "merged " << mergedRows << " rows\n";
	writeSupports(supports, argv[3]);
}