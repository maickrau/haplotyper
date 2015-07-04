//g++ identify_outliers.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o identify_outliers.exe
//./identify_outliers.exe inputFile outputFile numOfOutliers

#include <tuple>
#include <cassert>

#include "variant_utils.h"

std::vector<SNPSupport> markOutliers(std::vector<SNPSupport> in, int numOfOutliers)
{
	size_t numRows = 0;
	size_t numSNPs = 0;
	for (auto x : in)
	{
		numRows = std::max(numRows, x.readNum);
		numSNPs = std::max(numSNPs, x.SNPnum);
	}
	numRows++;
	numSNPs++;

	std::vector<SNPSupport> result = in;
	for (int i = 0; i < numOfOutliers; i++)
	{
		std::vector<std::tuple<size_t, size_t, size_t, size_t>> extents;
		extents.resize(numSNPs, std::make_tuple(-1, -1, 0, 0));
		for (auto x : result)
		{
			if (x.variant == 'X')
			{
				continue;
			}
			if (x.readNum <= std::get<0>(extents[x.SNPnum]))
			{
				std::get<1>(extents[x.SNPnum]) = std::get<0>(extents[x.SNPnum]);
				std::get<0>(extents[x.SNPnum]) = x.readNum;
			}
			if (x.readNum >= std::get<3>(extents[x.SNPnum]))
			{
				std::get<2>(extents[x.SNPnum]) = std::get<3>(extents[x.SNPnum]);
				std::get<3>(extents[x.SNPnum]) = x.readNum;
			}
		}
		for (auto& x : extents)
		{
			if (std::get<1>(x) == -1)
			{
				std::get<1>(x) = std::get<0>(x);
			}
			if (std::get<2>(x) == 0)
			{
				std::get<2>(x) = std::get<3>(x);
			}
		}
		size_t maxSNP = 0;
		size_t maxLength = 0;
		bool top = true;
		for (size_t i = 0; i < extents.size(); i++)
		{
			if (std::get<1>(extents[i])-std::get<0>(extents[i]) > maxLength)
			{
				maxLength = std::get<1>(extents[i])-std::get<0>(extents[i]);
				maxSNP = i;
				top = false;
			}
			if (std::get<3>(extents[i])-std::get<2>(extents[i]) > maxLength)
			{
				maxLength = std::get<3>(extents[i])-std::get<2>(extents[i]);
				maxSNP = i;
				top = true;
			}
		}
		size_t maxRead = 0;
		if (top)
		{
			maxRead = std::get<3>(extents[maxSNP]);
		}
		else
		{
			maxRead = std::get<0>(extents[maxSNP]);
		}
		bool found = false;
		for (auto& x : result)
		{
			if (x.readNum == maxRead && x.SNPnum == maxSNP)
			{
				x.variant = 'X';
				found = true;
				break;
			}
		}
		assert(found);
	}
	return result;
}

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);
	supports = markOutliers(supports, std::stoi(argv[3]));
	writeSupports(supports, argv[2]);
}