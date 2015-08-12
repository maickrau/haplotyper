//g++ prettify_sort.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o prettify_sort.exe
//./prettify_sort.exe inputSupports outputSupports renumbering

#include <algorithm>

#include "variant_utils.h"

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);
	size_t maxSNP = 0;
	size_t maxRead = 0;
	for (auto x : supports)
	{
		maxSNP = std::max(maxSNP, x.SNPnum);
		maxRead = std::max(maxRead, x.readNum);
	}
	maxSNP++;
	maxRead++;
	std::vector<std::pair<size_t, size_t>> rowStart;
	rowStart.resize(maxRead, {-1, 0});
	for (size_t i = 0; i < maxRead; i++)
	{
		rowStart[i].second = i;
	}
	for (auto x : supports)
	{
		rowStart[x.readNum].first = std::min(rowStart[x.readNum].first, x.SNPnum);
	}
	std::sort(rowStart.begin(), rowStart.end(), [](std::pair<size_t, size_t> left, std::pair<size_t, size_t> right) { return left.first < right.first; });
	SupportRenumbering renumbering;
	for (size_t i = 0; i < maxRead; i++)
	{
		renumbering.addReadRenumbering(rowStart[i].second, i);
	}
	for (size_t i = 0; i < maxSNP; i++)
	{
		renumbering.addSNPRenumbering(i, i);
	}
	std::vector<SNPSupport> renumbered = renumberSupports(supports, renumbering);
	writeSupports(renumbered, argv[2]);
	writeRenumbering(renumbering, argv[3]);
}