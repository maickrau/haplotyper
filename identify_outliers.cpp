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
		std::vector<std::tuple<size_t, size_t, size_t, size_t>> rowExtents;
		rowExtents.resize(numRows, std::make_tuple(-1, -1, 0, 0));
		for (auto x : result)
		{
			if (x.variant == 'X')
			{
				continue;
			}
			if (x.SNPnum <= std::get<0>(rowExtents[x.readNum]))
			{
				std::get<1>(rowExtents[x.readNum]) = std::get<0>(rowExtents[x.readNum]);
				std::get<0>(rowExtents[x.readNum]) = x.SNPnum;
			}
			else if (x.SNPnum <= std::get<1>(rowExtents[x.readNum]))
			{
				std::get<1>(rowExtents[x.readNum]) = x.SNPnum;
			}
			if (x.SNPnum >= std::get<3>(rowExtents[x.readNum]))
			{
				std::get<2>(rowExtents[x.readNum]) = std::get<3>(rowExtents[x.readNum]);
				std::get<3>(rowExtents[x.readNum]) = x.SNPnum;
			}
			else if (x.SNPnum >= std::get<2>(rowExtents[x.readNum]))
			{
				std::get<2>(rowExtents[x.readNum]) = x.SNPnum;
			}
		}
		for (auto& x : rowExtents)
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
		size_t maxRead = 0;
		size_t maxLength = 0;
		bool top = true;
		for (size_t i = 0; i < rowExtents.size(); i++)
		{
			if (std::get<1>(rowExtents[i])-std::get<0>(rowExtents[i]) > maxLength)
			{
				maxLength = std::get<1>(rowExtents[i])-std::get<0>(rowExtents[i]);
				maxRead = i;
				top = false;
			}
			if (std::get<3>(rowExtents[i])-std::get<2>(rowExtents[i]) > maxLength)
			{
				maxLength = std::get<3>(rowExtents[i])-std::get<2>(rowExtents[i]);
				maxRead = i;
				top = true;
			}
		}
		size_t maxSNP = 0;
		if (top)
		{
			maxSNP = std::get<3>(rowExtents[maxRead]);
		}
		else
		{
			maxSNP = std::get<0>(rowExtents[maxRead]);
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