//g++ remove_outliers.cpp variant_utils.cpp fasta_utils.cpp -o remove_outliers.exe -std=c++11
//./remove_outliers.exe inputFile outputFile maxAccidentalCoverage

#include <iostream>

#include "variant_utils.h"

size_t findBiggestAccidentalCoverage(std::vector<SNPSupport> supports, size_t minCoverage)
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
	std::vector<std::vector<bool>> used;
	used.resize(maxSNP);
	for (size_t i = 0; i < used.size(); i++)
	{
		used[i].resize(maxRead, false);
	}
	for (auto x : supports)
	{
		used[x.SNPnum][x.readNum] = true;
		rowExtents[x.readNum].first = std::min(rowExtents[x.readNum].first, x.SNPnum);
		rowExtents[x.readNum].second = std::max(rowExtents[x.readNum].second, x.SNPnum);
	}
	size_t biggestIndex = 0;
	size_t biggestCoverage = 0;
	for (size_t i = 0; i < used.size(); i++)
	{
		size_t coverage = 0;
		for (size_t j = 0; j < used[i].size(); j++)
		{
			if (!used[i][j] && rowExtents[j].first < i && rowExtents[j].second > i)
			{
				coverage++;
			}
		}
		if (coverage > biggestCoverage)
		{
			biggestCoverage = coverage;
			biggestIndex = i;
		}
	}
	if (biggestCoverage >= minCoverage)
	{
		return biggestIndex;
	}
	return -1;
}

std::vector<SNPLine> filterLines(std::vector<SNPLine> lines, size_t SNPnum)
{
	std::vector<SNPLine> ret;
	for (auto x : lines)
	{
		bool hasLeft = false;
		bool hasRight = false;
		for (auto y : x.variantsAtLocations)
		{
			if (y.first == SNPnum)
			{
				hasRight = false;
				hasLeft = false;
				break;
			}
			if (y.first < SNPnum)
			{
				hasLeft = true;
			}
			if (y.first > SNPnum)
			{
				hasRight = true;
			}
		}
		if (hasLeft && hasRight)
		{
			ret.push_back(x);
		}
	}
	return ret;
}

std::pair<size_t, bool> findEasiestRemovableLine(std::vector<SNPLine> lines, size_t SNPnum)
{
	std::pair<size_t, bool> easiestRemovable { -1, false };
	size_t easiestRemovableSize = -1;
	for (auto x : lines)
	{
		size_t leftSize = 0;
		size_t rightSize = 0;
		for (auto y : x.variantsAtLocations)
		{
			if (y.first < SNPnum)
			{
				leftSize++;
			}
			if (y.first > SNPnum)
			{
				rightSize++;
			}
		}
		if (leftSize < easiestRemovableSize)
		{
			easiestRemovable = std::pair<size_t, bool> { x.readNum, false };
			easiestRemovableSize = leftSize;
		}
		if (rightSize < easiestRemovableSize)
		{
			easiestRemovable = std::pair<size_t, bool> { x.readNum, true };
			easiestRemovableSize = rightSize;
		}
	}
	return easiestRemovable;
}

std::vector<SNPSupport> removeOutliers(std::vector<SNPSupport> supports, size_t readNum, size_t SNPnum, bool right)
{
	std::vector<SNPSupport> ret;
	for (auto x : supports)
	{
		if (x.readNum != readNum)
		{
			ret.push_back(x);
		}
		else
		{
			if (x.SNPnum < SNPnum && right)
			{
				ret.push_back(x);
			}
			else if (x.SNPnum > SNPnum && !right)
			{
				ret.push_back(x);
			}
		}
	}
	return ret;
}

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);
	size_t sizeStart = supports.size();
	size_t limit = std::stol(argv[3]);
	size_t foundSNP = findBiggestAccidentalCoverage(supports, limit);
	size_t oldSize = sizeStart;
	while (foundSNP != -1)
	{
		std::vector<SNPLine> lines = makeLines(supports);
		lines = filterLines(lines, foundSNP);
		std::pair<size_t, bool> easiestRemovable = findEasiestRemovableLine(lines, foundSNP);
		supports = removeOutliers(supports, easiestRemovable.first, foundSNP, easiestRemovable.second);
		size_t currentSize = supports.size();
		oldSize = currentSize;
		foundSNP = findBiggestAccidentalCoverage(supports, limit);
	}
	writeSupports(supports, argv[2]);
	size_t sizeEnd = supports.size();
	std::cerr << "removed " << sizeStart-sizeEnd << " outliers, from " << sizeStart << " to " << sizeEnd << "\n";
}