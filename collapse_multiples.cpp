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
	size_t maxRead = 0;
	for (auto x : supports)
	{
		maxSNP = std::max(maxSNP, x.SNPnum);
		maxRead = std::max(maxRead, x.readNum);
	}
	maxSNP++;
	maxRead++;
	renumbering = SupportRenumbering::identity(maxRead, maxSNP);

	std::cout << "start\n";
	std::sort(supports.begin(), supports.end(), [](SNPSupport left, SNPSupport right) { return left.SNPnum < right.SNPnum; });
	std::stable_sort(supports.begin(), supports.end(), [](SNPSupport left, SNPSupport right) { return left.readNum < right.readNum; });
	std::vector<SNPLine> rows = makeLines(supports);
	std::cout << "made lines\n";
	size_t lastRead = 0;
	std::vector<SNPLine> merged;
	merged.push_back(rows[0]);
	double currentSupport = 1;
	for (int i = 0; i < rows.size(); i++)
	{
		for (int j = 0; j < i; j++)
		{
			assert(rows[i].readNum != rows[j].readNum);
		}
	}
	std::cout << "merge\n";
	for (size_t i = 0; i < rows.size(); i++)
	{
		bool exists = false;
		for (size_t a = 0; a < merged.size(); a++)
		{
			if (rows[i] == merged[a])
			{
				exists = true;
				merged[a].merge(rows[i]);
				renumbering.overwriteReadRenumbering(rows[i].readNum, a);
				break;
			}
		}
		if (!exists)
		{
			merged.push_back(rows[i]);
			merged.back().readNum = merged.size()-1;
			renumbering.overwriteReadRenumbering(rows[i].readNum, merged.size()-1);
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
//	merged.second = renumberSupports(merged.second, merged.first);
	writeSupports(merged.second, argv[2]);
	writeRenumbering(merged.first, argv[3]);
}