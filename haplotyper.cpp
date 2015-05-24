#include <algorithm>
#include <iostream>
#include <cassert>
#include <chrono>

#include "variant_utils.h"
#include "haplotyper.h"

//must not return permutations, otherwise will produce about k! times more partitions than necessary
//[start, end)
std::vector<Partition> Partition::getAllPartitions(size_t start, size_t end, size_t k)
{
	std::vector<Partition> ret;
	std::vector<size_t> newPartition;
	newPartition.resize(end-start, 0);
	std::vector<size_t> rowsOccupied;
	rowsOccupied.resize(k+1, 0);
	rowsOccupied[0] = end-start;
	if (end == start+1)
	{
		ret.emplace_back(newPartition.begin(), newPartition.end());
		ret.back().minRow = start;
		ret.back().maxRow = end;
		ret.back().k = k;
		return ret;
	}
	size_t loc = end-start-1;
	bool repeat = true;
	while (repeat)
	{
		ret.emplace_back(newPartition.begin(), newPartition.end());
		ret.back().minRow = start;
		ret.back().maxRow = end;
		ret.back().k = k;
		assert(rowsOccupied[newPartition[loc]] > 0);
		assert(rowsOccupied[newPartition[loc]] <= end-start);
		rowsOccupied[newPartition[loc]]--;
		newPartition[loc]++;
		rowsOccupied[newPartition[loc]]++;
		assert(rowsOccupied[newPartition[loc]] > 0);
		assert(rowsOccupied[newPartition[loc]] <= end-start);
		if (rowsOccupied[newPartition[loc]-1] == 0 || newPartition[loc] == k)
		{
			while (rowsOccupied[newPartition[loc]-1] == 0 || newPartition[loc] == k)
			{
				assert(rowsOccupied[newPartition[loc]] > 0);
				assert(rowsOccupied[newPartition[loc]] <= end-start);
				rowsOccupied[newPartition[loc]]--;
				newPartition[loc] = 0;
				rowsOccupied[newPartition[loc]]++;
				assert(rowsOccupied[newPartition[loc]] > 0);
				assert(rowsOccupied[newPartition[loc]] <= end-start);
				loc--;
				if (loc == 0)
				{
					repeat = false;
					break;
				}
				assert(rowsOccupied[newPartition[loc]] > 0);
				assert(rowsOccupied[newPartition[loc]] <= end-start);
				rowsOccupied[newPartition[loc]]--;
				newPartition[loc]++;
				rowsOccupied[newPartition[loc]]++;
				assert(rowsOccupied[newPartition[loc]] > 0);
				assert(rowsOccupied[newPartition[loc]] <= end-start);
			}
			loc = end-start-1;
		}
	}
	return ret;
}

Partition Partition::merge(Partition second)
{
	Partition ret { *this };
	assert(maxRow >= second.minRow);
	assert(minRow <= second.maxRow);
	if (second.minRow < minRow)
	{
		ret.assignments.insert(ret.assignments.begin(), second.assignments.begin(), second.assignments.begin()+(minRow-second.minRow));
		ret.minRow = second.minRow;
	}
	if (second.maxRow > maxRow)
	{
		ret.assignments.insert(ret.assignments.end(), second.assignments.begin()+(maxRow-second.minRow), second.assignments.end());
		ret.maxRow = second.maxRow;
	}
	return ret;
}

double Partition::wCost(Column col, char variant, size_t haplotype)
{
	double sum = 0;
	assert(col.minRow == minRow);
	assert(col.maxRow >= maxRow);
	assert(col.maxRow <= maxRow);
	for (size_t i = col.minRow; i < col.maxRow; i++)
	{
		if (col.variants[i] != variant && assignments[i-minRow] == haplotype)
		{
			sum += col.costs[i];
		}
	}
	return sum;
}

double Partition::deltaCost(Column col)
{
	double sum = 0;
	for (size_t i = 0; i < k; i++)
	{
		double partSum = std::min(std::min(wCost(col, 'A', i), wCost(col, 'T', i)), std::min(wCost(col, 'C', i), wCost(col, 'G', i)));
		sum += partSum;
	}
	return sum;
}

bool Partition::extends(Partition second)
{
	size_t intersectionStart = std::max(minRow, second.minRow);
	size_t intersectionEnd = std::min(maxRow, second.maxRow);
	Partition thisSubset = filter(intersectionStart, intersectionEnd);
	Partition thatSubset = second.filter(intersectionStart, intersectionEnd);
	thisSubset.unpermutate();
	thatSubset.unpermutate();
	return std::equal(thisSubset.assignments.begin(), thisSubset.assignments.end(), thatSubset.assignments.begin());
}

Partition Partition::filter(size_t newMinRow, size_t newMaxRow)
{
	newMinRow = std::max(newMinRow, minRow);
	newMaxRow = std::min(newMaxRow, maxRow);
	Partition ret {assignments.begin()+(newMinRow-minRow), assignments.end()-(maxRow-newMaxRow)};
	ret.minRow = newMinRow;
	ret.maxRow = newMaxRow;
	ret.k = k;
	return ret;
}

void Partition::unpermutate()
{
	std::vector<size_t> mapping;
	mapping.resize(k, -1);
	size_t nextNum = 0;
	for (auto& x : assignments)
	{
		assert(x < k);
		if (mapping[x] == -1)
		{
			mapping[x] = nextNum;
			nextNum++;
		}
		x = mapping[x];
	}
}

template <typename Iterator>
Partition::Partition(Iterator start, Iterator end) :
	assignments(start, end),
	minRow(),
	maxRow(),
	k()
{
}

//returns optimal partition and its score
std::pair<Partition, double> haplotype(std::vector<SNPSupport> supports, size_t k)
{
	std::sort(supports.begin(), supports.end(), [](SNPSupport left, SNPSupport right) { return left.SNPnum < right.SNPnum; });

	//handle first column separately
	size_t firstSNPEnd = 0;
	size_t minCurrentReadNum = supports[0].readNum;
	size_t maxCurrentReadNum = supports[0].readNum;
	for (size_t i = 1; i < supports.size(); i++)
	{
		if (supports[i].SNPnum != supports[i-1].SNPnum)
		{
			firstSNPEnd = i;
			break;
		}
		else
		{
			minCurrentReadNum = std::min(minCurrentReadNum, supports[i].readNum);
			maxCurrentReadNum = std::max(maxCurrentReadNum, supports[i].readNum);
		}
	}
	//Column oldColumn { supports.begin(), supports.end(), supports[0].SNPnum };
	Column oldColumn { supports.begin(), supports.begin()+firstSNPEnd, supports[0].SNPnum };
	std::vector<Partition> oldRowPartitions = Partition::getAllPartitions(minCurrentReadNum, maxCurrentReadNum+1, k);
	std::vector<double> oldRowCosts;
	std::vector<Partition> oldOptimalPartitions = oldRowPartitions;
	for (auto x : oldRowPartitions)
	{
		oldRowCosts.push_back(x.deltaCost(oldColumn));
	}

	//each column, including last one
	minCurrentReadNum = supports[firstSNPEnd].readNum;
	maxCurrentReadNum = supports[firstSNPEnd].readNum;
	size_t thisSNPStart = firstSNPEnd;
	auto lastColumnTime = std::chrono::system_clock::now();
	for (size_t i = firstSNPEnd+1; i <= supports.size(); i++)
	{
		if (i == supports.size() || supports[i].SNPnum != supports[i-1].SNPnum)
		{
			auto newColumnTime = std::chrono::system_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::duration<int,std::milli>>(newColumnTime-lastColumnTime);
			lastColumnTime = newColumnTime;
			std::cerr << " " << diff.count() << "ms\n";
			std::cerr << "column " << supports[i-1].SNPnum;
			//Column col { supports.begin(), supports.end(), supports[thisSNPStart].SNPnum };
			Column col { supports.begin()+thisSNPStart, supports.begin()+i, supports[thisSNPStart].SNPnum };
			std::vector<Partition> newRowPartitions = Partition::getAllPartitions(minCurrentReadNum, maxCurrentReadNum+1, k);
			std::vector<double> newRowCosts;
			std::vector<Partition> newOptimalPartitions;
			for (auto x : newRowPartitions)
			{
				double cost = 0;
				bool hasCost = false;
				size_t optimalExtensionIndex = 0;
				for (size_t a = 0; a < oldRowPartitions.size(); a++)
				{
					if (x.extends(oldRowPartitions[a]))
					{
						if (!hasCost)
						{
							cost = oldRowCosts[a];
							hasCost = true;
							optimalExtensionIndex = a;
						}
						if (oldRowCosts[a] < cost)
						{
							optimalExtensionIndex = a;
						}
						cost = std::min(cost, oldRowCosts[a]);
					}
				}
				assert(hasCost);
				newRowCosts.push_back(cost+x.deltaCost(col));
				newOptimalPartitions.push_back(oldOptimalPartitions[optimalExtensionIndex].merge(x));
			}
			oldRowPartitions = newRowPartitions;
			oldRowCosts = newRowCosts;
			oldOptimalPartitions = newOptimalPartitions;
			if (i == supports.size())
			{
				break;
			}
			thisSNPStart = i;
			minCurrentReadNum = supports[i].readNum;
			maxCurrentReadNum = supports[i].readNum;
		}
		minCurrentReadNum = std::min(minCurrentReadNum, supports[i].readNum);
		maxCurrentReadNum = std::max(maxCurrentReadNum, supports[i].readNum);
	}

	//return only value instead of solution
	size_t optimalResultIndex = 0;
	for (size_t i = 1; i < oldRowCosts.size(); i++)
	{
		if (oldRowCosts[i] < oldRowCosts[optimalResultIndex])
		{
			optimalResultIndex = i;
		}
	}
	return std::pair<Partition, double>{oldOptimalPartitions[optimalResultIndex], oldRowCosts[optimalResultIndex]};
}
