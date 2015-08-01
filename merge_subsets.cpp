//g++ merge_subsets.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o merge_subsets.exe
//./merge_subsets.exe inputSupportsFile outputSupportsFile renumberingFile

#include <algorithm>
#include <iostream>

#include "variant_utils.h"

std::pair<std::vector<SNPSupport>, SupportRenumbering> mergeSubsets(std::vector<SNPSupport> supports)
{
	size_t maxSNP = 0;
	for (auto x : supports)
	{
		maxSNP = std::max(maxSNP, x.SNPnum);
	}
	maxSNP++;
	std::cout << "lines ";
	std::vector<SNPLine> lines = makeLines(supports);
	for (auto x : lines)
	{
		assert(x.variantsAtLocations.size() > 0);
	}
	std::cout << lines.size() << "\n";
	std::sort(lines.begin(), lines.end(), [](SNPLine left, SNPLine right) { return left.variantsAtLocations.size() > right.variantsAtLocations.size(); });
	SupportRenumbering renumbering;
	std::vector<SNPLine> merged;
	std::cout << "merge ";
	for (size_t i = 0; i < lines.size(); i++)
	{
		bool wasMerged = false;
		for (size_t a = 0; a < merged.size(); a++)
		{
			if (merged[a].contains(lines[i]))
			{
				merged[a].mergeSubset(lines[i]);
				wasMerged = true;
				renumbering.addReadRenumbering(lines[i].readNum, a);
				break;
			}
		}
		if (!wasMerged)
		{
			renumbering.addReadRenumbering(lines[i].readNum, merged.size());
			merged.push_back(lines[i]);
		}
	}
	std::cout << "to " << merged.size() << "\n";
	std::vector<SNPSupport> result;
	for (auto x : merged)
	{
		std::vector<SNPSupport> part = x.toSupports();
		result.insert(result.end(), part.begin(), part.end());
	}
	for (size_t i = 0; i < maxSNP; i++)
	{
		renumbering.addSNPRenumbering(i, i);
	}
	return std::pair<std::vector<SNPSupport>, SupportRenumbering> { result, renumbering };
}

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);
	std::pair<std::vector<SNPSupport>, SupportRenumbering> result = mergeSubsets(supports);
	writeSupports(result.first, argv[2]);
	writeRenumbering(result.second, argv[3]);
}