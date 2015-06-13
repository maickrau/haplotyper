//g++ remove_zero_columns.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o remove_zero_columns.exe
//./remove_zero_columns.exe inputSupportsFile outputSupportsFile outputRenumberingFile

#include "variant_utils.h"

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);
	std::vector<bool> readUsed;
	std::vector<bool> SNPused;
	for (auto x : supports)
	{
		if (readUsed.size() <= x.readNum)
		{
			readUsed.resize(x.readNum+1, false);
		}
		if (SNPused.size() <= x.SNPnum)
		{
			SNPused.resize(x.SNPnum+1, false);
		}
		readUsed[x.readNum] = true;
		SNPused[x.SNPnum] = true;
	}

	SupportRenumbering renumbering;
	size_t usedReads = 0;
	for (size_t i = 0; i < readUsed.size(); i++)
	{
		if (readUsed[i])
		{
			renumbering.addReadRenumbering(i, usedReads);
			usedReads++;
		}
	}
	size_t usedSNPs = 0;
	for (size_t i = 0; i < SNPused.size(); i++)
	{
		if (SNPused[i])
		{
			renumbering.addSNPRenumbering(i, usedSNPs);
			usedSNPs++;
		}
	}

	std::vector<SNPSupport> result = renumberSupports(supports, renumbering);
	writeSupports(result, argv[2]);
	writeRenumbering(renumbering, argv[3]);
}