//g++ active_rows.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o active_rows.exe

#include <iostream>

#include "variant_utils.h"

std::vector<std::pair<size_t, size_t>> necessaryAndIncidentalActives(std::vector<SNPSupport> supports)
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
	std::vector<std::vector<bool>> used;
	used.resize(maxSNP);
	for (size_t i = 0; i < maxSNP; i++)
	{
		used[i].resize(maxRead, false);
	}
	for (auto x : supports)
	{
		used[x.SNPnum][x.readNum] = true;
	}
	std::vector<std::pair<size_t, size_t>> ret;
	for (size_t i = 0; i < maxSNP; i++)
	{
		ret.emplace_back(0, 0);
		for (size_t j = 0; j < maxRead; j++)
		{
			if (used[i][j])
			{
				ret.back().first++;
				continue;
			}
			bool hasLeft = false;
			for (size_t k = 0; k < i; k++)
			{
				if (used[k][j])
				{
					hasLeft = true;
				}
			}
			for (size_t k = i+1; k < maxSNP && hasLeft; k++)
			{
				if (used[k][j])
				{
					ret.back().second++;
					hasLeft = false;
				}
			}
		}
	}
	return ret;
}

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);
	std::vector<std::pair<size_t, size_t>> actives = necessaryAndIncidentalActives(supports);
	for (size_t i = 0; i < actives.size(); i++)
	{
		std::cout << actives[i].first << "\t" << actives[i].second << "\t" << actives[i].first+actives[i].second << "\n";
	}
}