//g++ variant_utils.cpp fasta_utils.cpp proper_active_rows_per_snp.cpp -std=c++11 -o proper_active_rows_per_snp.exe

#include <iostream>

#include "variant_utils.h"

std::vector<size_t> calculateProperActivesPerSNP(const std::vector<SNPSupport>& supports)
{
	size_t numRows = 0;
	size_t numSNPs = 0;
	for (auto x : supports)
	{
		numRows = std::max(numRows, x.readNum);
		numSNPs = std::max(numSNPs, x.SNPnum);
	}
	numRows++;
	numSNPs++;
	std::vector<std::pair<size_t, size_t>> rowExtent;
	rowExtent.resize(numRows, std::pair<size_t, size_t>{-1, 0});
	for (auto x : supports)
	{
		rowExtent[x.readNum].first = std::min(rowExtent[x.readNum].first, x.SNPnum);
		rowExtent[x.readNum].second = std::max(rowExtent[x.readNum].second, x.SNPnum);
	}
	std::vector<size_t> ret;
	for (size_t i = 0; i < numSNPs; i++)
	{
		size_t num = 0;
		for (auto x : rowExtent)
		{
			if (x.first <= i && x.second >= i)
			{
				num++;
			}
		}
		ret.push_back(num);
	}
	return ret;
}

std::vector<size_t> calculateUsedActivesPerSNP(const std::vector<SNPSupport>& supports)
{
	size_t numRows = 0;
	size_t numSNPs = 0;
	for (auto x : supports)
	{
		numRows = std::max(numRows, x.readNum);
		numSNPs = std::max(numSNPs, x.SNPnum);
	}
	numRows++;
	numSNPs++;
	std::vector<std::pair<size_t, size_t>> columnExtent;
	columnExtent.resize(numSNPs, std::pair<size_t, size_t> {-1, 0});
	for (auto x : supports)
	{
		columnExtent[x.SNPnum].second = std::max(x.readNum, columnExtent[x.SNPnum].second);
		columnExtent[x.SNPnum].first = std::min(x.readNum, columnExtent[x.SNPnum].first);
	}
	for (size_t i = 1; i < columnExtent.size(); i++)
	{
		columnExtent[i].second = std::max(columnExtent[i-1].second, columnExtent[i].second);
	}
	for (size_t i = columnExtent.size()-2; i < columnExtent.size(); i--)
	{
		columnExtent[i].first = std::min(columnExtent[i+1].first, columnExtent[i].first);
	}
	std::vector<size_t> ret;
	for (size_t i = 0; i < columnExtent.size(); i++)
	{
		ret.push_back(columnExtent[i].second-columnExtent[i].first+1);
	}
	return ret;
}

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);
	std::vector<size_t> properPerSnp = calculateProperActivesPerSNP(supports);
	std::vector<size_t> usedPerSnp = calculateUsedActivesPerSNP(supports);
	for (size_t i = 0; i < properPerSnp.size(); i++)
	{
		std::cout << properPerSnp[i] << "\t" << usedPerSnp[i] << "\n";
	}
}
