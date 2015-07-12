#include <algorithm>
#include <iostream>
#include <cassert>
#include <chrono>
#include <valarray>
#include <tuple>
#include <cmath>
#include <cstring>

#include "variant_utils.h"
#include "haplotyper.h"

Partition::Partition() {}

PartitionAssignments::PartitionAssignmentElement::PartitionAssignmentElement(PartitionAssignments& container, size_t pos)	:
	container(container),
	pos(pos)
{
}

namespace std
{
	template <>
	void iter_swap(PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElement> left, PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElement> right)
	{
		size_t value = *left;
		*left = (size_t)*right;
		*right = value;
	}

}

void swap(PartitionAssignments::PartitionAssignmentElement& left, PartitionAssignments::PartitionAssignmentElement& right)
{
	size_t value = (size_t)left;
	left = (size_t)right;
	right = value;
}

PartitionAssignments::PartitionAssignmentElement::operator size_t() const
{
	size_t bytePos = pos*container.log2k/8;
	size_t bitPos = (pos*container.log2k)%8;
	size_t ret = 0;
	for (size_t i = 0; i < container.log2k; i++)
	{
		ret <<= 1;
		ret |= (container.data[bytePos]&(1<<bitPos))>>bitPos;
		bitPos++;
		if (bitPos == 8)
		{
			bitPos = 0;
			bytePos++;
		}
	}
	return ret;
}

PartitionAssignments::PartitionAssignmentElementConst::operator size_t() const
{
	return (size_t)el;
}

PartitionAssignments::PartitionAssignmentElementConst::PartitionAssignmentElementConst(const PartitionAssignments& container, size_t pos) :
	el(const_cast<PartitionAssignments&>(container), pos)
{
}

PartitionAssignments::PartitionAssignmentElementConst::PartitionAssignmentElementConst(PartitionAssignmentElement el) :
	el(el)
{
}

PartitionAssignments::PartitionAssignmentElement& PartitionAssignments::PartitionAssignmentElement::operator=(size_t value)
{
	assert(container.k > 0);
	assert(container.log2k > 0);
	size_t bytePos = ((pos+1)*container.log2k-1)/8;
	size_t bitPos = ((pos+1)*container.log2k-1)%8;
	size_t ret = 0;
	for (size_t i = 0; i < container.log2k; i++)
	{
		if (value & 1)
		{
			container.data[bytePos] |= 1<<bitPos;
		}
		else
		{
			container.data[bytePos] &= (-1) ^ (1<<bitPos);
		}
		value >>= 1;
		bitPos--;
		if (bitPos > 8)
		{
			bitPos = 7;
			bytePos--;
		}
	}
	return *this;
}

PartitionAssignments::PartitionAssignments() : k(0), log2k(0), actualSize(0), data()
{
}

PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElement> PartitionAssignments::begin()
{
	assert(log2k > 0);
	assert(k > 0);
	return PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElement>(*this, 0);
}

PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElement> PartitionAssignments::end()
{
	assert(log2k > 0);
	assert(k > 0);
	return PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElement>(*this, actualSize);
}

PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElementConst> PartitionAssignments::begin() const
{
	assert(log2k > 0);
	assert(k > 0);
	return PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElementConst>(const_cast<PartitionAssignments&>(*this), 0);
}

PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElementConst> PartitionAssignments::end() const
{
	assert(log2k > 0);
	assert(k > 0);
	return PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElementConst>(const_cast<PartitionAssignments&>(*this), actualSize);
}

void PartitionAssignments::push_back(size_t value)
{
	assert(log2k > 0);
	assert(k > 0);
	if (capacity() <= size())
	{
		extendCapacity(size()*2, 0);
	}
	(*this)[size()] = value;
	actualSize++;
}

void PartitionAssignments::extendCapacity(size_t newCapacity, size_t defaultValue)
{
	assert(log2k > 0);
	assert(k > 0);
	assert(defaultValue == 0);
	size_t newDataSize = (newCapacity*log2k/8+7)+1;
	data.resize(newDataSize, defaultValue);
}

void PartitionAssignments::setk(size_t newk)
{
	k = newk;
	log2k = ceil(log2(k));
	assert(data.size() == 0);
}

void PartitionAssignments::resize(size_t newSize, size_t defaultValue)
{
	assert(log2k > 0);
	assert(k > 0);
	assert(defaultValue == 0);
	size_t newDataSize = (newSize*log2k/8+7)+1;
	data.resize(newDataSize, defaultValue);
	actualSize = newSize;
}

size_t PartitionAssignments::size() const
{
	assert(log2k > 0);
	assert(k > 0);
	return actualSize;
}

size_t PartitionAssignments::capacity() const
{
	return data.size()*8/log2k;
}

PartitionAssignments::PartitionAssignmentElement PartitionAssignments::operator[](size_t pos)
{
	return PartitionAssignmentElement(*this, pos);
}

PartitionAssignments::PartitionAssignmentElementConst PartitionAssignments::operator[](size_t pos) const
{
	return PartitionAssignmentElementConst(*this, pos);
}



//preserve numbering of overlapping haplotypes
//eg {0, 1, 2}
//      {0, 1, 0}
//=>
//   {0, 1, 2, 1}
std::vector<size_t> getNumbering(std::vector<size_t> left, std::vector<size_t> right, size_t k)
{
	assert(left.size() == right.size());
	std::vector<size_t> numbering;
	numbering.resize(k, -1);
	std::set<size_t> unusedNumberings;
	for (size_t i = 0; i < k; i++)
	{
		unusedNumberings.insert(i);
	}
	for (size_t i = 0; i < left.size(); i++)
	{
		assert(left[i] < k);
		assert(right[i] < k);
		std::set<size_t>::iterator found = unusedNumberings.find(left[i]);
		if (found != unusedNumberings.end())
		{
			unusedNumberings.erase(found);
		}
		else
		{
			assert(numbering[right[i]] == left[i]);
		}
		numbering[right[i]] = left[i];
	}
	for (size_t i = 0; i < k; i++)
	{
		if (numbering[i] == -1)
		{
			std::set<size_t>::iterator found = unusedNumberings.lower_bound(0);
			assert(found != unusedNumberings.end());
			numbering[i] = *found;
			unusedNumberings.erase(found);
		}
	}
	return numbering;
}

void verifyNumbering(std::vector<size_t> numbering, size_t k)
{
	std::vector<bool> used;
	used.resize(k, false);
	for (auto x : numbering)
	{
		assert(x < k);
		assert(!used[x]);
		used[x] = true;
	}
}

class PartitionContainer
{
public:
	PartitionContainer(size_t k) : k(k) {};
	size_t insertPartition(Partition partition);
	size_t extendPartition(size_t partitionNum, Partition extension);
	Partition getPartition(size_t partitionNum);
	void clearUnused(std::vector<size_t> used);
private:
	size_t k;
	//haplotype, row number, previous index (rownumber-1's index)
	std::vector<std::tuple<size_t, size_t, size_t>> partitionsLinkedList;
	std::vector<size_t> unusedIndices;
};

size_t PartitionContainer::insertPartition(Partition partition)
{
	size_t lastIndex = partitionsLinkedList.size();
	partitionsLinkedList.emplace_back(partition.assignments[0], partition.minRow, -1);
	for (size_t i = 1; i < partition.assignments.size(); i++)
	{
		partitionsLinkedList.emplace_back(partition.assignments[i], partition.minRow+i, lastIndex+i-1);
	}
	return partitionsLinkedList.size()-1;
}

size_t PartitionContainer::extendPartition(size_t partitionNum, Partition extension)
{
	assert(std::get<1>(partitionsLinkedList[partitionNum]) >= extension.minRow);
	assert(std::get<1>(partitionsLinkedList[partitionNum]) <= extension.maxRow);

	//find numbering
	std::vector<size_t> left;
	size_t pos = partitionNum;
	while (std::get<1>(partitionsLinkedList[pos]) >= extension.minRow)
	{
		left.push_back(std::get<0>(partitionsLinkedList[pos]));
		if (std::get<2>(partitionsLinkedList[pos]) == -1)
		{
			break;
		}
		size_t oldRowNum = std::get<1>(partitionsLinkedList[pos]);
		pos = std::get<2>(partitionsLinkedList[pos]);
		size_t newRowNum = std::get<1>(partitionsLinkedList[pos]);
		assert(newRowNum == oldRowNum-1);
		assert(oldRowNum > 0);
	}
	std::reverse(left.begin(), left.end());
	assert(std::get<1>(partitionsLinkedList[partitionNum])+1 >= left.size());
	assert(std::get<1>(partitionsLinkedList[partitionNum]) >= extension.minRow);
	assert(std::get<1>(partitionsLinkedList[partitionNum])-left.size()+1 >= extension.minRow);
	std::vector<size_t> right { extension.assignments.begin()+(std::get<1>(partitionsLinkedList[partitionNum])+1-left.size()-extension.minRow), extension.assignments.begin()+(std::get<1>(partitionsLinkedList[partitionNum])-extension.minRow+1) };
	assert(left.size() == right.size());
	std::vector<size_t> numbering = getNumbering(left, right, k);

	Partition leftPartition{left.begin(), left.end(), k};
	Partition rightPartition{right.begin(), right.end(), k};
	assert(leftPartition.extends(rightPartition));
	verifyNumbering(numbering, k);

	//insert part after overlap
	pos = partitionNum;
	size_t extensionIndex = std::get<1>(partitionsLinkedList[partitionNum])-extension.minRow+1;
	while (extensionIndex < extension.assignments.size())
	{
		if (unusedIndices.size() > 0)
		{
			size_t insertPos = unusedIndices.back();
			unusedIndices.pop_back();
			std::get<0>(partitionsLinkedList[insertPos]) = numbering[extension.assignments[extensionIndex]];
			std::get<1>(partitionsLinkedList[insertPos]) = extension.minRow+extensionIndex;
			std::get<2>(partitionsLinkedList[insertPos]) = pos;
			pos = insertPos;
		}
		else
		{
			partitionsLinkedList.emplace_back(numbering[extension.assignments[extensionIndex]], extension.minRow+extensionIndex, pos);
			pos = partitionsLinkedList.size()-1;
		}
		extensionIndex++;
	}
	return pos;
}

Partition PartitionContainer::getPartition(size_t partitionNum)
{
	Partition ret;
	ret.maxRow = std::get<1>(partitionsLinkedList[partitionNum]);
	ret.minRow = ret.maxRow;
	ret.setk(k);
	while (std::get<2>(partitionsLinkedList[partitionNum]) != -1)
	{
		ret.assignments.push_back(std::get<0>(partitionsLinkedList[partitionNum]));
		size_t oldMinRow = ret.minRow;
		partitionNum = std::get<2>(partitionsLinkedList[partitionNum]);
		ret.minRow = std::get<1>(partitionsLinkedList[partitionNum]);
		size_t newMinRow = ret.minRow;
		assert(newMinRow == oldMinRow-1);
	}
	ret.assignments.push_back(std::get<0>(partitionsLinkedList[partitionNum]));
	ret.minRow = std::get<1>(partitionsLinkedList[partitionNum]);
	std::reverse(ret.assignments.begin(), ret.assignments.end());
	return ret;
}

void PartitionContainer::clearUnused(std::vector<size_t> used)
{
	std::vector<bool> isUsed;
	isUsed.resize(partitionsLinkedList.size(), false);
	for (auto x : used)
	{
		while (std::get<2>(partitionsLinkedList[x]) != -1)
		{
			isUsed[x] = true;
			x = std::get<2>(partitionsLinkedList[x]);
		}
		isUsed[x] = true;
	}
	for (auto x : unusedIndices)
	{
		isUsed[x] = true;
	}
	for (size_t i = 0; i < partitionsLinkedList.size(); i++)
	{
		if (!isUsed[i])
		{
			unusedIndices.push_back(i);
		}
	}
}

size_t Partition::getk() const
{
	return k;
}

void Partition::setk(size_t newK)
{
	k = newK;
	assignments.setk(k);
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
		ret.emplace_back(newPartition.begin(), newPartition.end(), k);
		ret.back().minRow = start;
		ret.back().maxRow = end;
		return ret;
	}
	size_t loc = end-start-1;
	bool repeat = true;
	while (repeat)
	{
		ret.emplace_back(newPartition.begin(), newPartition.end(), k);
		ret.back().minRow = start;
		ret.back().maxRow = end;
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

Partition Partition::merge(Partition second) const
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

	size_t intersectionStart = std::max(minRow, second.minRow);
	size_t intersectionEnd = std::min(maxRow, second.maxRow);
	std::vector<size_t> leftPart { assignments.begin()+(intersectionStart-minRow), assignments.begin()+(intersectionEnd-minRow) };
	std::vector<size_t> rightPart { second.assignments.begin()+(intersectionStart-second.minRow), second.assignments.begin()+(intersectionEnd-second.minRow) };
	std::vector<size_t> numbering = getNumbering(leftPart, rightPart, k);
	
	verifyNumbering(numbering, k);

	for (size_t i = ret.maxRow; i < second.maxRow; i++)
	{
		ret.assignments.push_back(numbering[second.assignments[i-second.minRow]]);
	}
	ret.maxRow = second.maxRow;

	return ret;
}

double Partition::deltaCost(const Column& col) const
{
	std::vector<std::array<size_t, 4>> costs;
	std::vector<size_t> costSum;
	costs.resize(k, {0, 0, 0, 0});
	costSum.resize(k, 0);
	for (size_t i = col.minRow; i < col.maxRow; i++)
	{
		switch(col.variants[i])
		{
		case 'A':
			costs[assignments[i-minRow]][0] += col.costs[i];
			costSum[assignments[i-minRow]] += col.costs[i];
			break;
		case 'T':
			costs[assignments[i-minRow]][1] += col.costs[i];
			costSum[assignments[i-minRow]] += col.costs[i];
			break;
		case 'C':
			costs[assignments[i-minRow]][2] += col.costs[i];
			costSum[assignments[i-minRow]] += col.costs[i];
			break;
		case 'G':
			costs[assignments[i-minRow]][3] += col.costs[i];
			costSum[assignments[i-minRow]] += col.costs[i];
			break;
		default:
			break;
		}
	}
	size_t totalCost = 0;
	for (size_t i = 0; i < k; i++)
	{
		totalCost += costSum[i]-std::max(std::max(costs[i][0], costs[i][1]), std::max(costs[i][2], costs[i][3]));
	}
	return totalCost;
}

bool Partition::extends(Partition second) const
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
bool partitionCompare(const Partition& left, const Partition& right)
{
	int compare = memcmp((const char*)left.assignments.data.data(), (const char*)right.assignments.data.data(), left.assignments.actualSize*left.assignments.log2k/8);
	if (compare < 0)
	{
		return true;
	}
	if (compare > 0)
	{
		return false;
	}
	for (size_t i = (left.assignments.actualSize*left.assignments.log2k/8)*8/left.assignments.log2k; i < left.assignments.actualSize; i++)
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
std::vector<std::vector<size_t>> findExtensions(const std::vector<Partition>& lastRow, const std::vector<Partition>& newRow)
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
	std::sort(newLastRow.begin(), newLastRow.end(), [](const std::pair<size_t, Partition>& left, const std::pair<size_t, Partition>& right) { return partitionCompare(left.second, right.second); });
	std::sort(newNewRow.begin(), newNewRow.end(), [](const std::pair<size_t, Partition>& left, const std::pair<size_t, Partition>& right) { return partitionCompare(left.second, right.second); });

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

bool verifyExtensions(const std::vector<Partition>& lastRow, const std::vector<Partition>& newRow, const std::vector<std::vector<size_t>>& extensions)
{
	for (size_t i = 0; i < extensions.size(); i++)
	{
		for (auto y : extensions[i])
		{
			if (!newRow[i].extends(lastRow[y]))
			{
				return false;
			}
		}
	}
	return true;
}

Partition Partition::filter(size_t newMinRow, size_t newMaxRow) const
{
	newMinRow = std::max(newMinRow, minRow);
	newMaxRow = std::min(newMaxRow, maxRow);
	Partition ret {assignments.begin()+(newMinRow-minRow), assignments.end()-(maxRow-newMaxRow), k};
	ret.minRow = newMinRow;
	ret.maxRow = newMaxRow;
	return ret;
}

void Partition::unpermutate()
{
	std::vector<size_t> mapping;
	mapping.resize(k, -1);
	size_t nextNum = 0;
	for (size_t i = 0; i < assignments.size(); i++)
	{
		assert(assignments[i] < k);
		if (mapping[assignments[i]] == -1)
		{
			mapping[assignments[i]] = nextNum;
			nextNum++;
		}
		assignments[i] = mapping[assignments[i]];
	}
}

template <typename Iterator>
Partition::Partition(Iterator start, Iterator end, size_t k) :
	assignments(),
	minRow(),
	maxRow(),
	k(k)
{
	assignments.setk(k);
	size_t pos = 0;
	while (start != end)
	{
		assignments.push_back((size_t)*start);
		pos++;
		start++;
	}
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
	PartitionContainer optimalPartitions { k };
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

			assert(verifyExtensions(oldRowPartitions, newRowPartitions, extensions));

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
