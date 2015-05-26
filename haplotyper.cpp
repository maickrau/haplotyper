#include <algorithm>
#include <iostream>
#include <cassert>
#include <chrono>
#include <valarray>

#include "variant_utils.h"
#include "haplotyper.h"

class PartitionContainer
{
public:
	size_t insertPartition(Partition partition);
	size_t extendPartition(size_t partitionNum, Partition extension);
	Partition getPartition(size_t partitionNum);
	void clearUnused(std::vector<size_t> used);
private:
	std::vector<size_t> unusedIndices;
	std::vector<Partition> partitions;
};

size_t PartitionContainer::insertPartition(Partition partition)
{
	partitions.push_back(partition);
	return partitions.size()-1;
}

size_t PartitionContainer::extendPartition(size_t partitionNum, Partition extension)
{
	Partition left = getPartition(partitionNum);
	Partition extend = left.merge(extension);
	if (unusedIndices.size() > 0)
	{
		partitions[unusedIndices.back()] = extend;
		size_t ret = unusedIndices.back();
		unusedIndices.pop_back();
		return ret;
	}
	partitions.push_back(extend);
	return partitions.size()-1;
}

Partition PartitionContainer::getPartition(size_t partitionNum)
{
	return partitions[partitionNum];
}

void PartitionContainer::clearUnused(std::vector<size_t> used)
{
	std::vector<bool> isUsed;
	isUsed.resize(partitions.size(), false);
	for (auto x : used)
	{
		isUsed[x] = true;
	}
	//mark unused as used so they don't get inserted twice into unusedIndices
	for (auto x : unusedIndices)
	{
		isUsed[x] = true;
	}
	for (size_t i = 0; i < partitions.size(); i++)
	{
		if (!isUsed[i])
		{
			unusedIndices.push_back(i);
		}
	}
}

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
/*	if (maxRow < second.minRow)
	{
		for (size_t i = maxRow; i < second.minRow; i++)
		{
			ret.assignments.push_back(0);
		}
		ret.maxRow = second.minRow;
	}*/
	assert(maxRow >= second.minRow);
	assert(minRow <= second.minRow);
	
	//preserve numbering of overlapping haplotypes
	//eg {0, 1, 2}
	//      {0, 1, 0}
	//=>
	//   {0, 1, 2, 1}
	std::vector<size_t> numbering;
	numbering.resize(k+1, -1);
	std::set<size_t> unusedNumberings;
	for (size_t i = 0; i < k; i++)
	{
		unusedNumberings.insert(i);
	}
	for (size_t i = std::max(minRow, second.minRow); i < std::min(maxRow, second.maxRow); i++)
	{
		numbering[second.assignments[i-second.minRow]] = ret.assignments[i-ret.minRow];
		auto found = unusedNumberings.find(second.assignments[i-second.minRow]);
		if (found != unusedNumberings.end())
		{
			unusedNumberings.erase(found);
		}
	}
	for (size_t i = 0; i < k; i++)
	{
		if (numbering[i] == -1)
		{
			auto found = unusedNumberings.lower_bound(0);
			assert(found != unusedNumberings.end());
			numbering[i] = *found;
			unusedNumberings.erase(found);
		}
	}
	
	for (size_t i = ret.maxRow; i < second.maxRow; i++)
	{
		ret.assignments.push_back(numbering[second.assignments[i-second.minRow]]);
	}
	ret.maxRow = second.maxRow;

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
		if (col.variants[i] != 0 && col.variants[i] != variant && assignments[i-minRow] == haplotype)
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
	if (intersectionEnd <= intersectionStart)
	{
		return true;
	}
	std::valarray<size_t> mapping(-1, k);
	std::valarray<size_t> reverseMapping(-1, k);
	for (size_t i = intersectionStart; i < intersectionEnd; i++)
	{
		if (mapping[assignments[i-minRow]] == -1)
		{
			mapping[assignments[i-minRow]] = second.assignments[i-second.minRow];
		}
		if (mapping[assignments[i-minRow]] != second.assignments[i-second.minRow])
		{
			return false;
		}
		if (reverseMapping[second.assignments[i-second.minRow]] == -1)
		{
			reverseMapping[second.assignments[i-second.minRow]] = assignments[i-minRow];
		}
		if (reverseMapping[second.assignments[i-second.minRow]] != assignments[i-minRow])
		{
			return false;
		}
	}
	return true;
}

//basically <
bool partitionCompare(Partition left, Partition right)
{
	for (size_t i = 0; i < left.assignments.size(); i++)
	{
		if (left.assignments[i] < right.assignments[i])
		{
			return true;
		}
		if (left.assignments[i] > right.assignments[i])
		{
			return false;
		}
	}
	return false;
}

//first index is new partition index, second index (inner vector) is old partition indices
std::vector<std::vector<size_t>> findExtensions(std::vector<Partition> lastRow, std::vector<Partition> newRow)
{
	std::vector<std::vector<size_t>> ret;
	ret.resize(newRow.size());
	size_t intersectionStart = std::max(lastRow[0].minRow, newRow[0].minRow);
	size_t intersectionEnd = std::min(lastRow[0].maxRow, newRow[0].maxRow);
	if (intersectionEnd <= intersectionStart)
	{
		for (size_t i = 0; i < newRow.size(); i++)
		{
			for (size_t j = 0; j < lastRow.size(); j++)
			{
				ret[i].push_back(j);
			}
		}
		return ret;
	}
	std::vector<std::pair<size_t, Partition>> newLastRow;
	std::vector<std::pair<size_t, Partition>> newNewRow;
	if (intersectionStart != lastRow[0].minRow || intersectionEnd != lastRow[0].maxRow)
	{
		for (size_t i = 0; i < lastRow.size(); i++)
		{
			Partition insertion = lastRow[i].filter(intersectionStart, intersectionEnd);
			insertion.unpermutate();
			newLastRow.emplace_back(i, insertion);
		}
	}
	else
	{
		for (size_t i = 0; i < lastRow.size(); i++)
		{
			newLastRow.emplace_back(i, lastRow[i]);
		}
	}
	if (intersectionStart != newRow[0].minRow || intersectionEnd != newRow[0].maxRow)
	{
		for (size_t i = 0; i < newRow.size(); i++)
		{
			Partition insertion = newRow[i].filter(intersectionStart, intersectionEnd);
			insertion.unpermutate();
			newNewRow.emplace_back(i, insertion);
		}
	}
	else
	{
		for (size_t i = 0; i < newRow.size(); i++)
		{
			newNewRow.emplace_back(i, newRow[i]);
		}
	}
	std::sort(newLastRow.begin(), newLastRow.end(), [](std::pair<size_t, Partition> left, std::pair<size_t, Partition> right) { return partitionCompare(left.second, right.second); });
	std::sort(newNewRow.begin(), newNewRow.end(), [](std::pair<size_t, Partition> left, std::pair<size_t, Partition> right) { return partitionCompare(left.second, right.second); });

	size_t lastIndex = 0;
	size_t newIndex = 0;
	while (lastIndex < newLastRow.size() && newIndex < newNewRow.size())
	{
		if (partitionCompare(newLastRow[lastIndex].second, newNewRow[newIndex].second))
		{
			lastIndex++;
		}
		else if (partitionCompare(newNewRow[newIndex].second, newLastRow[lastIndex].second))
		{
			newIndex++;
		}
		else
		{
			size_t newIndexEnd = newIndex+1;
			size_t lastIndexEnd = lastIndex+1;
			while (newIndexEnd < newNewRow.size() && !partitionCompare(newNewRow[newIndexEnd-1].second, newNewRow[newIndexEnd].second))
			{
				newIndexEnd++;
			}
			while (lastIndexEnd < newLastRow.size() && !partitionCompare(newLastRow[lastIndexEnd-1].second, newLastRow[lastIndexEnd].second))
			{
				lastIndexEnd++;
			}
			for (size_t i = newIndex; i < newIndexEnd; i++)
			{
				for (size_t j = lastIndex; j < lastIndexEnd; j++)
				{
					ret[newNewRow[i].first].push_back(newLastRow[j].first);
				}
			}
			newIndex = newIndexEnd;
			lastIndex = lastIndexEnd;
		}
	}
	return ret;
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

std::vector<std::pair<size_t, size_t>> getActiveRows(std::vector<SNPSupport> supports)
{
	std::vector<std::pair<size_t, size_t>> ret;
	std::vector<bool> hasValue;
	for (auto x : supports)
	{
		if (ret.size() <= x.SNPnum)
		{
			ret.resize(x.SNPnum+1);
			hasValue.resize(x.SNPnum+1, false);
		}
		if (!hasValue[x.SNPnum])
		{
			ret[x.SNPnum].first = x.readNum;
			ret[x.SNPnum].second = x.readNum;
			hasValue[x.SNPnum] = true;
		}
		ret[x.SNPnum].first = std::min(ret[x.SNPnum].first, x.readNum);
		ret[x.SNPnum].second = std::max(ret[x.SNPnum].second, x.readNum);
	}
	for (size_t i = 0; i < ret.size(); i++)
	{
		if (!hasValue[i])
		{
			ret[i].first = -1;
			ret[i].second = 0;
		}
	}
	for (size_t i = 1; i < ret.size(); i++)
	{
		ret[i].second = std::max(ret[i].second, ret[i-1].second);
	}
	for (size_t i = ret.size()-2; i > 0; i--)
	{
		ret[i].first = std::min(ret[i].first, ret[i+1].first);
	}
	return ret;
}

//returns optimal partition and its score
std::pair<Partition, double> haplotype(std::vector<SNPSupport> supports, size_t k)
{
	std::cerr << "sort\n";
	std::sort(supports.begin(), supports.end(), [](SNPSupport left, SNPSupport right) { return left.SNPnum < right.SNPnum; });

	std::cerr << "active columns\n";
	std::vector<std::pair<size_t, size_t>> activeRowsPerColumn = getActiveRows(supports);

	std::cerr << "column " << supports[0].SNPnum << " (" << activeRowsPerColumn[supports[0].SNPnum].first << "-" << activeRowsPerColumn[supports[0].SNPnum].second << ")\n";
	//handle first column separately
	size_t firstSNPEnd = 0;
	for (size_t i = 1; i < supports.size(); i++)
	{
		if (supports[i].SNPnum != supports[i-1].SNPnum)
		{
			firstSNPEnd = i;
			break;
		}
	}
	PartitionContainer optimalPartitions;
	//Column oldColumn { supports.begin(), supports.end(), supports[0].SNPnum };
	Column oldColumn { supports.begin(), supports.begin()+firstSNPEnd, supports[0].SNPnum, activeRowsPerColumn[supports[0].SNPnum].first, activeRowsPerColumn[supports[0].SNPnum].second+1 };
	std::vector<Partition> oldRowPartitions = Partition::getAllPartitions(activeRowsPerColumn[supports[0].SNPnum].first, activeRowsPerColumn[supports[0].SNPnum].second+1, k);
	std::vector<double> oldRowCosts;
	std::vector<size_t> oldOptimalPartitions;
	for (size_t i = 0; i < oldRowPartitions.size(); i++)
	{
		oldOptimalPartitions.push_back(optimalPartitions.insertPartition(oldRowPartitions[i]));
	}
	for (auto x : oldRowPartitions)
	{
		oldRowCosts.push_back(x.deltaCost(oldColumn));
	}

	std::cerr << "start real columns\n";
	//each column, including last one
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
			std::cerr << "column " << supports[i-1].SNPnum << " (" << activeRowsPerColumn[supports[thisSNPStart].SNPnum].first << "-" << activeRowsPerColumn[supports[thisSNPStart].SNPnum].second << ")";
			//Column col { supports.begin(), supports.end(), supports[thisSNPStart].SNPnum };
			Column col { supports.begin()+thisSNPStart, supports.begin()+i, supports[thisSNPStart].SNPnum, activeRowsPerColumn[supports[thisSNPStart].SNPnum].first, activeRowsPerColumn[supports[thisSNPStart].SNPnum].second+1 };
			std::vector<Partition> newRowPartitions = Partition::getAllPartitions(activeRowsPerColumn[supports[thisSNPStart].SNPnum].first, activeRowsPerColumn[supports[thisSNPStart].SNPnum].second+1, k);
			std::vector<double> newRowCosts;
			std::vector<size_t> newOptimalPartitions;
			std::vector<std::vector<size_t>> extensions = findExtensions(oldRowPartitions, newRowPartitions);
			for (size_t i = 0; i < extensions.size(); i++)
			{
				double cost = 0;
				bool hasCost = false;
				size_t optimalExtensionIndex = 0;
				for (size_t a = 0; a < extensions[i].size(); a++)
				{
						if (!hasCost)
						{
							cost = oldRowCosts[extensions[i][a]];
							hasCost = true;
							optimalExtensionIndex = extensions[i][a];
						}
						if (oldRowCosts[extensions[i][a]] < cost)
						{
							optimalExtensionIndex = extensions[i][a];
						}
						cost = std::min(cost, oldRowCosts[extensions[i][a]]);
				}
				assert(hasCost);
				newRowCosts.push_back(cost+newRowPartitions[i].deltaCost(col));
				newOptimalPartitions.push_back(optimalPartitions.extendPartition(oldOptimalPartitions[optimalExtensionIndex], newRowPartitions[i]));
			}
			oldRowPartitions = newRowPartitions;
			oldRowCosts = newRowCosts;
			oldOptimalPartitions = newOptimalPartitions;
			if (i == supports.size())
			{
				break;
			}
			thisSNPStart = i;
		}
	}

	size_t optimalResultIndex = 0;
	for (size_t i = 1; i < oldRowCosts.size(); i++)
	{
		if (oldRowCosts[i] < oldRowCosts[optimalResultIndex])
		{
			optimalResultIndex = i;
		}
	}
	return std::pair<Partition, double>{optimalPartitions.getPartition(oldOptimalPartitions[optimalResultIndex]), oldRowCosts[optimalResultIndex]};
}
