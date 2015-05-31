//g++ visualize_snpsupport.cpp variant_utils.cpp fasta_utils.cpp -o visualize_snpsupport.exe -std=c++11
//./visualize_snpsupport.exe supportsFile outputFile

#include "variant_utils.h"

int main(int argc, char** argv)
{
	std::vector<SNPSupport> support = loadSupports(argv[1]);
	//[rownum][snpnum]
	std::vector<std::vector<char>> supports;
	size_t maxRow = 0;
	size_t maxSNPnum = 0;
	for (auto x : support)
	{
		maxRow = std::max(maxRow, x.readNum);
		maxSNPnum = std::max(maxSNPnum, x.SNPnum);
	}
	supports.resize(maxRow+1);
	for (size_t i = 0; i < supports.size(); i++)
	{
		supports[i].resize(maxSNPnum+1, 0);
	}
	for (auto x : support)
	{
		supports[x.readNum][x.SNPnum] = x.variant;
	}

	std::ofstream file {argv[2]};
	for (size_t i = 0; i < supports.size(); i++)
	{
		for (size_t j = 0; j < supports[i].size(); j++)
		{
			if (supports[i][j] != 0)
			{
				file << supports[i][j];
			}
			file << "\t";
		}
		file << "\n";
	}
}