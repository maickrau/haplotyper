//./collapse_multiples.exe inputFile outputFile renumberingFile
//g++ collapse_multiples.cpp variant_utils.cpp -std=c++11 -o collapse_multiples.exe

#include <cassert>
#include <algorithm>
#include <iostream>
#include <unordered_map>

#include "variant_utils.h"

std::pair<SupportRenumbering, std::vector<SNPSupport>> mergeSupports(std::vector<SNPSupport> supports)
{
	SupportRenumbering renumbering;
	size_t maxSNP = 0;
	for (auto x : supports)
	{
		maxSNP = std::max(maxSNP, x.SNPnum);
	}
	maxSNP++;
	for (size_t i = 0; i < maxSNP; i++)
	{
		renumbering.addSNPRenumbering(i, i);
	}

	std::cout << "start\n";
	std::sort(supports.begin(), supports.end(), [](SNPSupport left, SNPSupport right) { return left.SNPnum < right.SNPnum; });
	std::stable_sort(supports.begin(), supports.end(), [](SNPSupport left, SNPSupport right) { return left.readNum < right.readNum; });
	std::vector<SNPLine> rows;
	size_t maxReadNum = 0;
	for (auto x : supports)
	{
		maxReadNum = std::max(maxReadNum, x.readNum);
	}
	rows.resize(maxReadNum+1);
	std::cout << "make lines\n";
	size_t lastReadStart = 0;
	for (size_t i = 1; i < supports.size(); i++)
	{
		if (supports[i].readNum != supports[i-1].readNum)
		{
			rows[supports[lastReadStart].readNum] = SNPLine {supports.begin()+lastReadStart, supports.begin()+i, supports[lastReadStart].readNum};
			lastReadStart = i;
		}
	}
	rows[supports[lastReadStart].readNum] = SNPLine {supports.begin()+lastReadStart, supports.end(), supports[lastReadStart].readNum};
	std::cout << "made lines\n";
	size_t lastRead = 0;
	std::vector<SNPLine> merged;
	merged.push_back(rows[0]);
	double currentSupport = 1;
	renumbering.addReadRenumbering(0, 0);
	std::cout << "merge\n";
	for (size_t i = 1; i < rows.size(); i++)
	{
		bool exists = false;
		for (size_t a = 0; a < merged.size(); a++)
		{
			if (rows[i] == merged[a])
			{
				exists = true;
				merged[a].merge(rows[i]);
				renumbering.addReadRenumbering(i, a);
				break;
			}
		}
		if (!exists)
		{
			merged.push_back(rows[i]);
			renumbering.addReadRenumbering(i, merged.size()-1);
		}
	}
	std::cout << "merged from " << rows.size() << " rows to " << merged.size() << " rows\n";
	std::cout << "get snpsupports\n";
	std::vector<SNPSupport> ret;
	for (auto x : merged)
	{
		std::vector<SNPSupport> newSupports = x.toSupports();
		ret.insert(ret.end(), newSupports.begin(), newSupports.end());
	}
	std::stable_sort(ret.begin(), ret.end(), [](SNPSupport left, SNPSupport right) { return left.SNPnum < right.SNPnum; });
	std::cout << "return\n";
	return std::pair<SupportRenumbering, std::vector<SNPSupport>> { renumbering, ret };
}

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);
	auto merged = mergeSupports(supports);
	merged.second = renumberSupports(merged.second, merged.first);
	writeSupports(merged.second, argv[2]);
	writeRenumbering(merged.first, argv[3]);
}