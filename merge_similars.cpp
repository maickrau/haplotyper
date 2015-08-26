//g++ merge_similars.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o merge_similars.exe
//./merge_similars.exe inputSupportsFile necessaryCoverageLimit totalCoverageLimit outputSupportsFile renumberingFile

#include <unordered_map>
#include <algorithm>
#include <iostream>

#include "variant_utils.h"

size_t findSNPWithHighestNecessaryCoverage(const std::vector<SNPSupport>& supports, size_t minCoverage)
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

size_t findSNPWithHighestTotalCoverage(const std::vector<SNPSupport>& supports, size_t minCoverage)
{
	size_t maxSNP = 0;
	size_t maxRead = 0;
	for (auto x : supports)
	{
		maxSNP = std::max(maxSNP, x.SNPnum);
		maxRead = std::max(maxRead, x.readNum);
	}
	maxSNP++;
	maxRead++;
	std::vector<std::pair<size_t, size_t>> rowExtents;
	rowExtents.resize(maxRead, {-1, 0});
	for (auto x : supports)
	{
		rowExtents[x.readNum].first = std::min(rowExtents[x.readNum].first, x.SNPnum);
		rowExtents[x.readNum].second = std::max(rowExtents[x.readNum].second, x.SNPnum);
	}
	size_t maxIndex = 0;
	size_t maxCoverage = 0;
	for (size_t i = 0; i < maxSNP; i++)
	{
		size_t currentCoverage = 0;
		for (size_t j = 0; j < maxRead; j++)
		{
			if (rowExtents[j].first <= i && rowExtents[j].second >= i)
			{
				currentCoverage++;
			}
		}
		if (currentCoverage > maxCoverage)
		{
			maxCoverage = currentCoverage;
			maxIndex = i;
		}
	}
	if (maxCoverage >= minCoverage)
	{
		return maxIndex;
	}
	return -1;
}

std::vector<SNPLine> filterLines(const std::vector<SNPLine>& lines, size_t SNPposition)
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

std::vector<SNPLine> filterLinesAny(const std::vector<SNPLine>& lines, size_t SNPposition)
{
	std::vector<SNPLine> ret;
	for (auto x : lines)
	{
		if (x.variantsAtLocations[0].first <= SNPposition && x.variantsAtLocations.back().first >= SNPposition)
		{
			ret.push_back(x);
		}
	}
	return ret;
}

double lineDifference(SNPLine left, SNPLine right)
{
	std::set<size_t> leftSNPs;
	for (auto x : left.variantsAtLocations)
	{
		leftSNPs.insert(x.first);
	}
	double result = 0;
	for (size_t i = 0; i < left.supportsAtLocations.size(); i++)
	{
		result += left.supportsAtLocations[i];
	}
	for (auto x : right.variantsAtLocations)
	{
		if (leftSNPs.count(x.first) > 0)
		{
			if (left.variantAt(x.first) != x.second)
			{
				result += right.supportAt(x.first)*2+left.supportAt(x.first);
			}
			else
			{
				result -= left.supportAt(x.first);
			}
		}
		else
		{
			result += right.supportAt(x.first);
		}
	}
	return result;

}

template <typename RowFilter>
std::pair<size_t, size_t> findMostSimilarRows(const std::vector<SNPSupport>& supports, size_t SNPposition, RowFilter filter)
{
	std::vector<SNPLine> lines = makeLines(supports);
	lines = filter(lines, SNPposition);
	size_t bestLeft = lines[1].readNum;
	size_t bestRight = lines[0].readNum;
	double bestDifference = lineDifference(lines[1], lines[0]);
	for (size_t i = 0; i < lines.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			double difference = lineDifference(lines[i], lines[j]);
			assert(lineDifference(lines[j], lines[i]) == difference);
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
	size_t necessaryCoverageLimit = std::stoi(argv[2]);
	size_t totalCoverageLimit = std::stoi(argv[3]);
	size_t mergedRows = 0;
	SupportRenumbering renumbering;
	size_t maxSNP = 0;
	size_t maxRead = 0;
	for (auto x : supports)
	{
		maxSNP = std::max(maxSNP, x.SNPnum);
		maxRead = std::max(maxRead, x.readNum);
	}
	maxSNP++;
	maxRead++;
	for (size_t i = 0; i < maxSNP; i++)
	{
		renumbering.addSNPRenumbering(i, i);
	}
	for (size_t i = 0; i < maxRead; i++)
	{
		renumbering.addReadRenumbering(i, i);
	}
	std::cerr << maxRead << " lines\n";
	while (true)
	{
		size_t SNPposition = findSNPWithHighestNecessaryCoverage(supports, necessaryCoverageLimit);
		if (SNPposition == -1)
		{
			break;
		}
		std::pair<size_t, size_t> similars = findMostSimilarRows(supports, SNPposition, filterLines);
		supports = mergeRowsForceMerge(supports, similars.first, similars.second);
		for (size_t i = 0; i < maxRead; i++)
		{
			if (renumbering.getReadRenumbering(i) == similars.second)
			{
				renumbering.overwriteReadRenumbering(i, similars.first);
			}
		}
		mergedRows++;
	}
	while (true)
	{
		size_t SNPposition = findSNPWithHighestTotalCoverage(supports, totalCoverageLimit);
		if (SNPposition == -1)
		{
			break;
		}
		std::pair<size_t, size_t> similars = findMostSimilarRows(supports, SNPposition, filterLinesAny);
		supports = mergeRowsForceMerge(supports, similars.first, similars.second);
		for (size_t i = 0; i < maxRead; i++)
		{
			if (renumbering.getReadRenumbering(i) == similars.second)
			{
				renumbering.overwriteReadRenumbering(i, similars.first);
			}
		}
		mergedRows++;
	}
	std::cerr << "merged " << mergedRows << " rows\n";
	writeSupports(supports, argv[4]);
	writeRenumbering(renumbering, argv[5]);
}