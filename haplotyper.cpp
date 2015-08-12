#include <algorithm>
#include <iostream>
#include <cassert>
#include <chrono>
#include <tuple>
#include <cmath>
#include <cstring>
#include <set>

#include "variant_utils.h"
#include "haplotyper.h"

template <typename T>
class AssignOnce
{
public:
	AssignOnce() : assigned(false) {};
	operator=(const T& v)
	{
		assert(!assigned);
		value = v;
		assigned = true;
	}
	operator T()
	{
		assert(assigned);
		return value;
	}
private:
	T value;
	bool assigned;
};

AssignOnce<size_t> k;
AssignOnce<size_t> log2k;

//preserve numbering of overlapping haplotypes
//eg {0, 1, 2}
//      {0, 1, 0}
//=>
//   {0, 1, 2, 1}
std::vector<size_t> getNumbering(const std::vector<size_t>& left, const std::vector<size_t>& right, size_t k)
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

bool setEqual(const std::set<size_t>& left, const std::set<size_t>& right)
{
	if (left.size() != right.size())
	{
		return false;
	}
	return std::equal(left.begin(), left.end(), right.begin());
}

std::set<size_t> subsetIndices(const std::set<size_t>& subset, const std::set<size_t>& set)
{
	auto subsetIter = subset.begin();
	auto setIter = set.begin();
	size_t index = 0;
	std::set<size_t> ret;
	while (setIter != set.end() && subsetIter != subset.end())
	{
		if (*setIter == *subsetIter)
		{
			ret.insert(index);
			subsetIter++;
		}
		setIter++;
		index++;
	}
	assert(subsetIter == subset.end());
	return ret;
}

std::set<size_t> setIntersection(const std::set<size_t>& left, const std::set<size_t>& right)
{
	std::set<size_t> ret;
	std::set_intersection(left.begin(), left.end(), right.begin(), right.end(), std::inserter(ret, ret.begin()));
	return ret;
}

std::set<size_t> setUnion(const std::set<size_t>& left, const std::set<size_t>& right)
{
	std::set<size_t> ret;
	std::set_union(left.begin(), left.end(), right.begin(), right.end(), std::inserter(ret, ret.begin()));
	return ret;
}

int fact(int n)
{
	int ret = 1;
	for (int i = 2; i < n; i++)
	{
		ret *= i;
	}
	return ret;
}

TinyVectorMemoryAllocator::~TinyVectorMemoryAllocator()
{
	if (memory != nullptr)
	{
		delete [] memory;
	}
}

TinyVectorMemoryAllocator::TinyVectorMemoryAllocator(TinyVectorMemoryAllocator&& second) :
	bytesPerVector(second.bytesPerVector),
	nextEmpty(second.nextEmpty),
	memory(second.memory),
	size(second.size)
{
	second.memory = nullptr;
}

TinyVectorMemoryAllocator& TinyVectorMemoryAllocator::operator=(TinyVectorMemoryAllocator&& second)
{
	if (memory != nullptr)
	{
		delete [] memory;
	}
	memory = second.memory;
	bytesPerVector = second.bytesPerVector;
	nextEmpty = second.nextEmpty;
	size = second.size;
	second.memory = nullptr;
}

size_t getApproxNumberOfPartitions(size_t coverage, size_t k)
{
	return ceil((double)pow(k, coverage)/(double)fact(k))+k;
}

TinyVectorMemoryAllocator::TinyVectorMemoryAllocator(size_t coverage, size_t length, size_t k) :
	bytesPerVector(0),
	nextEmpty(0),
	memory(nullptr),
	size(0)
{
	size_t numOfVectors = getApproxNumberOfPartitions(coverage, k);
	bytesPerVector = (length*ceil(log2(k))+7)/8;
	size = numOfVectors*bytesPerVector;
	memory = new unsigned char[size];
	memset(memory, 0, size);
//	memory.resize(size, 0);
	nextEmpty = 0;
}

unsigned char* TinyVectorMemoryAllocator::allocate(size_t allocateSize)
{
	assert(memory != nullptr);
//	assert(size == bytesPerVector);
	if (nextEmpty+allocateSize > size)
	{
		std::ofstream file{"temp.temp"};
		int a = 1;
		int b = allocateSize+a;
		auto c = memory+b;
	}
	assert(nextEmpty+allocateSize <= size);
	unsigned char* ret = memory+nextEmpty;
	nextEmpty += allocateSize;
	return ret;
}

void TinyVectorMemoryAllocator::empty()
{
	delete [] memory;
	memory = nullptr;
//	(std::vector<unsigned char>{}).swap(memory);
}

SparsePartition::SparsePartition() :
	inner()
{
}

SparsePartition::SparsePartition(SolidPartition inner) :
	inner(inner)
{
}

size_t SparsePartition::getAssignment(size_t loc, const std::set<size_t>& actives) const
{
	assert(actives.size() == inner.assignments.size());
	size_t index = 0;
	auto iter = actives.begin();
	while (*iter != loc)
	{
		assert(iter != actives.end());
		assert(index < inner.assignments.size());
		iter++;
		index++;
		assert(iter != actives.end());
		assert(index < inner.assignments.size());
	}
	assert(index < inner.assignments.size());
	return inner.assignments[index];
}

std::vector<SparsePartition> SparsePartition::getAllPartitions(std::set<size_t> actives, TinyVectorMemoryAllocator& allocator)
{
	std::vector<SolidPartition> solids = SolidPartition::getAllPartitions(0, actives.size(), allocator);
	std::vector<SparsePartition> ret;
	for (size_t i = 0; i < solids.size()/2; i++)
	{
		ret.emplace_back(solids.back());
		solids.pop_back();
	}
	solids.shrink_to_fit();
	for (auto x : solids)
	{
		ret.emplace_back(x);
		assert(x.assignments.size() == actives.size());
	}
	return ret;
}

SolidPartition SparsePartition::getSolidFromIndices(const std::set<size_t>& pickThese, TinyVectorMemoryAllocator& allocator) const
{
	assert(pickThese.size() > 0);
	SolidPartition ret;
	ret.assignments.reserve(pickThese.size(), allocator);
#ifndef NDEBUG
	ret.assignments.actualSize = pickThese.size();
#endif
	size_t num = 0;
	for (auto x : pickThese)
	{
		assert(num < pickThese.size());
		ret.assignments[num] = inner.assignments[x];
		num++;
	}
	assert(num == pickThese.size());
	return ret;
}

SolidPartition SparsePartition::getSolid(const std::set<size_t>& actives, TinyVectorMemoryAllocator& allocator) const
{
	size_t min = *actives.lower_bound(0);
	size_t max = *actives.upper_bound(-1);
	for (size_t i = *actives.lower_bound(0); i < *actives.upper_bound(-1); i++)
	{
		assert(actives.count(i) == 1);
	}
	SolidPartition subset = getSubset(actives, actives, allocator);
	return subset;
}

double SparsePartition::deltaCost(const Column& col, const std::set<size_t>& actives) const
{
	std::vector<std::array<double, 4>> costs;
	std::vector<double> costSum;
	costs.resize(k, {0, 0, 0, 0});
	costSum.resize(k, 0);
	auto iter = actives.begin();
	size_t index = 0;
	while (iter != actives.end())
	{
		size_t assignment = inner.assignments[index];
		assert(assignment < costs.size());
		assert(assignment < costSum.size());
		size_t columnIndex = *iter;
		assert(columnIndex < col.costs.size());
		switch(col.variants[columnIndex])
		{
			case 'A':
				costs[assignment][0] += col.costs[columnIndex];
				costSum[assignment] += col.costs[columnIndex];
				break;
			case 'T':
				costs[assignment][1] += col.costs[columnIndex];
				costSum[assignment] += col.costs[columnIndex];
				break;
			case 'C':
				costs[assignment][2] += col.costs[columnIndex];
				costSum[assignment] += col.costs[columnIndex];
				break;
			case 'G':
				costs[assignment][3] += col.costs[columnIndex];
				costSum[assignment] += col.costs[columnIndex];
				break;
			default:
				break;
		}
		index++;
		iter++;
	}
	size_t totalCost = 0;
	for (size_t i = 0; i < k; i++)
	{
		totalCost += costSum[i]-std::max(std::max(costs[i][0], costs[i][1]), std::max(costs[i][2], costs[i][3]));
	}
	return totalCost;
}

SolidPartition SparsePartition::getSubset(const std::set<size_t>& subset, const std::set<size_t>& actives, TinyVectorMemoryAllocator& allocator) const
{
	assert(subset.size() > 0);
	SolidPartition ret;
	ret.assignments.reserve(subset.size(), allocator);
#ifndef NDEBUG
	ret.assignments.actualSize = subset.size();
#endif
	size_t num = 0;
	for (auto x : subset)
	{
		assert(num < subset.size());
		assert(num < ret.assignments.size());
		assert(getAssignment(x, actives) < k);
		ret.assignments[num] = getAssignment(x, actives);
		num++;
	}
	assert(num == subset.size());
	return ret;
}


SolidPartition::SolidPartition() : assignments() {}

PartitionAssignments::PartitionAssignmentElement::PartitionAssignmentElement(PartitionAssignments& container, size_t pos) :
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

size_t PartitionAssignments::PartitionAssignmentElement::getValue() const
{
	size_t bytePos = pos*log2k/8;
	size_t bitPos = (pos*log2k)%8;
	assert(container.data.size() > bytePos);
	size_t ret = 0;
	unsigned char mask;
	mask = 255-((1<<(bitPos))-1);
	ret = container.data[bytePos] & mask;
	ret >>= bitPos;
	if (bitPos+log2k <= 8)
	{
		ret &= (1<<(log2k))-1;
		return ret;
	}
	assert(bitPos < 8);
	size_t bitsEaten = 8-bitPos;
	for (; bitsEaten+8 <= log2k; bitsEaten += 8)
	{
		bytePos++;
		assert(bytePos < container.data.size());
		ret += ((size_t)container.data[bytePos]) << bitsEaten;
	}
	if (bitsEaten < log2k)
	{
		bytePos++;
		assert(bytePos < container.data.size());
		ret += ((size_t)(container.data[bytePos]&((1<<(log2k-bitsEaten))-1))) << bitsEaten;
	}
	return ret;
}

PartitionAssignments::PartitionAssignmentElement::operator size_t() const
{
	return getValue();
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
	assert(k > 0);
	assert(log2k > 0);
	size_t bytePos = ((pos)*log2k)/8;
	size_t bitPos = ((pos)*log2k)%8;
	if (container.data.size() <= bytePos)
	{
		int a = 1;
	}
	assert(container.data.size() > bytePos);
	for (size_t i = 0; i < log2k; i++)
	{
		assert(bytePos < container.data.size());
		if (value & 1)
		{
			container.data[bytePos] |= 1<<bitPos;
		}
		else
		{
			container.data[bytePos] &= (-1) ^ (1<<bitPos);
		}
		value >>= 1;
		if (i == log2k-1)
		{
			break;
		}
		bitPos++;
		if (bitPos >= 8)
		{
			bitPos = 0;
			bytePos++;
		}
	}
	if (container.data.size() <= bytePos)
	{
		int a = 1;
	}
	assert(container.data.size() > bytePos);
	return *this;
}

PartitionAssignments::PartitionAssignments() : 
#ifndef NDEBUG
	actualSize(0), 
#endif
	data()
{
}

PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElement> PartitionAssignments::begin()
{
	assert(log2k > 0);
	assert(k > 0);
	return PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElement>(*this, 0);
}

PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElement> PartitionAssignments::end(size_t size)
{
#ifndef NDEBUG
	assert(size == actualSize);
#endif
	assert(log2k > 0);
	assert(k > 0);
	return PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElement>(*this, size);
}

PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElementConst> PartitionAssignments::begin() const
{
	assert(log2k > 0);
	assert(k > 0);
	return PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElementConst>(const_cast<PartitionAssignments&>(*this), 0);
}

PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElementConst> PartitionAssignments::end(size_t size) const
{
#ifndef NDEBUG
	assert(size == actualSize);
#endif
	assert(log2k > 0);
	assert(k > 0);
	return PartitionAssignments::iterator<PartitionAssignments::PartitionAssignmentElementConst>(const_cast<PartitionAssignments&>(*this), size);
}

void PartitionAssignments::extendCapacity(size_t newCapacity, size_t defaultValue, TinyVectorMemoryAllocator& allocator)
{
	assert(log2k > 0);
	assert(k > 0);
	assert(defaultValue == 0);
	assert(newCapacity > 0);
	size_t newDataSize = (newCapacity*log2k+7)/8;
	data.resize(newDataSize, defaultValue, allocator);
	assert(capacity() >= newCapacity);
}

#ifndef NDEBUG

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

#endif

void PartitionAssignments::reserve(size_t newSize, TinyVectorMemoryAllocator& allocator)
{
	extendCapacity(newSize, 0, allocator);
}

PartitionAssignments::PartitionAssignmentElement PartitionAssignments::operator[](size_t pos)
{
	assert(pos < size());
	return PartitionAssignmentElement(*this, pos);
}

PartitionAssignments::PartitionAssignmentElementConst PartitionAssignments::operator[](size_t pos) const
{
	assert(pos < size());
	return PartitionAssignmentElementConst(*this, pos);
}

class SparsePartitionContainer
{
public:
	SparsePartitionContainer(const std::vector<std::set<size_t>>& readsPerSNP, size_t k);
	size_t insertPartition(const SparsePartition& partition, size_t SNPnum, size_t size);
	size_t extendPartition(size_t partitionNum, const SparsePartition& extension, size_t SNPnum, size_t maxSNP, const std::set<size_t>& oldActives, const std::set<size_t>& partitionActives, const std::set<size_t>& intersection);
	SparsePartition getPartition(size_t partitionNum, size_t maxSNP, TinyVectorMemoryAllocator& allocator) const;
	void clearUnused(std::vector<size_t> used);
	void reserveMore(size_t moreIndices);
private:
	std::vector<size_t> getPartitionAssignments(size_t partitionNum, const std::set<size_t>& indexes) const;
	size_t extend(size_t index, const std::vector<size_t>& extension);
	std::vector<size_t> unusedIndices;
	//assignment, position, previous index
	std::vector<std::tuple<unsigned char, size_t, size_t>> partitionsLinkedList;
	std::vector<size_t> readOrdering;
	std::vector<size_t> inverseReadOrdering;
	std::vector<size_t> SNPstarts;
	size_t k;
};

void SparsePartitionContainer::reserveMore(size_t moreIndices)
{
	if (moreIndices <= unusedIndices.size())
	{
		return;
	}
	partitionsLinkedList.reserve(partitionsLinkedList.size()+moreIndices-unusedIndices.size());
}

SparsePartitionContainer::SparsePartitionContainer(const std::vector<std::set<size_t>>& readsPerSNP, size_t k) :
	unusedIndices(),
	partitionsLinkedList(),
	readOrdering(),
	inverseReadOrdering(),
	k(k)
{
	std::set<size_t> usedReads;
	for (auto x : readsPerSNP)
	{
		SNPstarts.push_back(readOrdering.size());
		for (auto y : x)
		{
			if (usedReads.count(y) == 0)
			{
				usedReads.insert(y);
				readOrdering.push_back(y);
			}
		}
	}
	for (size_t i = 0; i < usedReads.size(); i++)
	{
		assert(usedReads.count(i) > 0);
	}
	assert(readOrdering.size() == (*std::max_element(readOrdering.begin(), readOrdering.end()))+1);
	assert(usedReads.size() == (*usedReads.upper_bound(-1)));
	SNPstarts.push_back(readOrdering.size());
	inverseReadOrdering.resize(readOrdering.size(), -1);
	for (size_t i = 0; i < readOrdering.size(); i++)
	{
		assert(readOrdering[i] < inverseReadOrdering.size());
		inverseReadOrdering[readOrdering[i]] = i;
	}
}

size_t SparsePartitionContainer::insertPartition(const SparsePartition& partition, size_t SNPnum, size_t size)
{
	assert(SNPstarts[SNPnum] == 0);
#ifndef NDEBUG
	assert(size == partition.inner.assignments.size());
#endif
	std::vector<size_t> extension;
	for (size_t i = 0; i < size; i++)
	{
		extension.push_back(partition.inner.assignments[i]);
	}
	return extend(-1, extension);
}

size_t SparsePartitionContainer::extendPartition(size_t partitionNum, const SparsePartition& partition, size_t SNPnum, size_t maxSNP, const std::set<size_t>& oldActives, const std::set<size_t>& partitionActives, const std::set<size_t>& intersection)
{
	size_t numNews = SNPstarts[SNPnum+1]-SNPstarts[SNPnum];
	assert(numNews == partitionActives.size()-intersection.size());
	if (numNews == 0)
	{
		return partitionNum;
	}
	assert(partitionNum < partitionsLinkedList.size());
	std::vector<size_t> rightNumbering;
	rightNumbering.reserve(intersection.size());
	for (auto x : intersection)
	{
		rightNumbering.push_back(partition.getAssignment(x, partitionActives));
	}
	std::vector<size_t> leftNumbering = getPartitionAssignments(partitionNum, intersection);
	std::vector<size_t> numbering = getNumbering(leftNumbering, rightNumbering, k);
	std::vector<size_t> extension;
	extension.reserve(SNPstarts[SNPnum+1]-SNPstarts[SNPnum]);
	for (size_t i = SNPstarts[SNPnum]; i < SNPstarts[SNPnum+1]; i++)
	{
		extension.push_back(numbering[partition.getAssignment(readOrdering[i], partitionActives)]);
	}
	assert(SNPstarts[SNPnum] == oldActives.size());
	assert(oldActives.size() + extension.size() == setUnion(oldActives, partitionActives).size());
	size_t result = extend(partitionNum, extension);
	return result;
}

size_t SparsePartitionContainer::extend(size_t index, const std::vector<size_t>& extension)
{
	size_t position = 0;
	if (index != -1)
	{
		position = std::get<1>(partitionsLinkedList[index])+1;
	}
	for (size_t i = 0; i < extension.size(); i++)
	{
		size_t putHere;
		if (unusedIndices.size() > 0)
		{
			putHere = unusedIndices.back();
			unusedIndices.pop_back();
		}
		else
		{
			putHere = partitionsLinkedList.size();
			partitionsLinkedList.emplace_back(-1, -1, -1);
		}
		assert(putHere < partitionsLinkedList.size());
		std::get<0>(partitionsLinkedList[putHere]) = extension[i];
		std::get<1>(partitionsLinkedList[putHere]) = position;
		std::get<2>(partitionsLinkedList[putHere]) = index;
		position++;
		index = putHere;
	}
	return index;
}

std::vector<size_t> SparsePartitionContainer::getPartitionAssignments(size_t partitionNum, const std::set<size_t>& indexes) const
{
	std::vector<size_t> ret;
	ret.resize(indexes.size(), -1);
	std::set<size_t> remaining = indexes;
	size_t index = partitionNum;
	size_t position;
	while (index != -1 && remaining.size() > 0)
	{
		position = std::get<1>(partitionsLinkedList[index]);
		if (remaining.count(readOrdering[position]) > 0)
		{
			remaining.erase(readOrdering[position]);
			size_t pos = 0;
			for (auto iter = indexes.begin(); iter != indexes.end(); iter++)
			{
				if (*iter == readOrdering[position])
				{
					ret[pos] = std::get<0>(partitionsLinkedList[index]);
					break;
				}
				pos++;
			}
		}
		index = std::get<2>(partitionsLinkedList[index]);
	}
	assert(remaining.size() == 0);
	assert(std::none_of(ret.begin(), ret.end(), [](size_t x) { return x == -1;}));

	return ret;
}

SparsePartition SparsePartitionContainer::getPartition(size_t partitionNum, size_t maxSNP, TinyVectorMemoryAllocator& allocator) const
{
	std::vector<size_t> assignments;
	assert(partitionNum < partitionsLinkedList.size());
	size_t index = partitionNum;
	size_t position = std::get<1>(partitionsLinkedList[index]);
	assignments.resize(position+1, -1);
	while (index != -1)
	{
		assert(index == partitionNum ^ std::get<1>(partitionsLinkedList[index])+1 == position);
		position = std::get<1>(partitionsLinkedList[index]);
		assignments[position] = std::get<0>(partitionsLinkedList[index]);
		index = std::get<2>(partitionsLinkedList[index]);
	}
	assert(position == 0);

	SparsePartition ret;
	ret.inner.assignments.reserve(assignments.size(), allocator);
#ifndef NDEBUG
	ret.inner.assignments.actualSize = assignments.size();
#endif
	size_t num = 0;
	for (size_t i = 0; i < inverseReadOrdering.size(); i++)
	{
		if (inverseReadOrdering[i] < assignments.size())
		{
			assert(assignments[inverseReadOrdering[i]] != -1);
			assert(num < assignments.size());
			assert(num < ret.inner.assignments.size());
			ret.inner.assignments[num] = assignments[inverseReadOrdering[i]];
			num++;
		}
	}
	assert(num == assignments.size());
	return ret;
}

void SparsePartitionContainer::clearUnused(std::vector<size_t> used)
{
	std::vector<bool> isUsed;
	isUsed.resize(partitionsLinkedList.size(), false);
	for (auto x : used)
	{
		assert(x < partitionsLinkedList.size());
		isUsed[x] = true;
		while (std::get<2>(partitionsLinkedList[x]) != -1 && !isUsed[std::get<2>(partitionsLinkedList[x])])
		{
			x = std::get<2>(partitionsLinkedList[x]);
			assert(x < partitionsLinkedList.size());
			isUsed[x] = true;
		}
	}
	for (auto x : unusedIndices)
	{
		assert(x < partitionsLinkedList.size());
		isUsed[x] = true;
/*		while (std::get<2>(partitionsLinkedList[x]) != -1 && !isUsed[std::get<2>(partitionsLinkedList[x])])
		{
			x = std::get<2>(partitionsLinkedList[x]);
			assert(x < partitionsLinkedList.size());
			isUsed[x] = true;
		}*/
	}
	while (partitionsLinkedList.size() > 0 && !isUsed[partitionsLinkedList.size()-1])
	{
		partitionsLinkedList.pop_back();
	}
	for (size_t i = 0; i < partitionsLinkedList.size(); i++)
	{
		if (!isUsed[i])
		{
			unusedIndices.push_back(i);
		}
	}
	unusedIndices.shrink_to_fit();
	partitionsLinkedList.shrink_to_fit();
}

//must not return permutations, otherwise will produce about k! times more partitions than necessary
//[start, end)
std::vector<SolidPartition> SolidPartition::getAllPartitions(size_t start, size_t end, TinyVectorMemoryAllocator& allocator)
{
	assert(end > start);
	std::vector<SolidPartition> ret;
	ret.reserve(getApproxNumberOfPartitions(end-start, k));
	std::vector<size_t> newPartition;
	newPartition.resize(end-start+5, 0);
	std::vector<size_t> rowsOccupied;
	rowsOccupied.resize(k+5, 0);
	rowsOccupied[0] = end-start;
	if (end == start+1)
	{
		ret.emplace_back(newPartition.begin(), newPartition.begin()+(end-start), end-start, allocator);
		return ret;
	}
	size_t loc = end-start-1;
	bool repeat = true;
	while (repeat)
	{
		ret.emplace_back(newPartition.begin(), newPartition.begin()+(end-start), end-start, allocator);
		assert(loc < end-start);
		assert(newPartition[loc] < k);
		assert(rowsOccupied[newPartition[loc]] > 0);
		assert(rowsOccupied[newPartition[loc]] <= end-start);
		rowsOccupied[newPartition[loc]]--;
		newPartition[loc]++;
		rowsOccupied[newPartition[loc]]++;
		assert(loc < end-start);
		assert(newPartition[loc] <= k);
		assert(rowsOccupied[newPartition[loc]] > 0);
		assert(rowsOccupied[newPartition[loc]] <= end-start);
		if (rowsOccupied[newPartition[loc]-1] == 0 || newPartition[loc] == k)
		{
			while (rowsOccupied[newPartition[loc]-1] == 0 || newPartition[loc] == k)
			{
				assert(loc < end-start);
				assert(newPartition[loc] <= k);
				assert(rowsOccupied[newPartition[loc]] > 0);
				assert(rowsOccupied[newPartition[loc]] <= end-start);
				rowsOccupied[newPartition[loc]]--;
				newPartition[loc] = 0;
				rowsOccupied[newPartition[loc]]++;
				assert(loc < end-start);
				assert(newPartition[loc] <= k);
				assert(rowsOccupied[newPartition[loc]] > 0);
				assert(rowsOccupied[newPartition[loc]] <= end-start);
				assert(loc < end-start);
				loc--;
				if (loc == 0)
				{
					repeat = false;
					break;
				}
				assert(loc < end-start);
				assert(newPartition[loc] <= k);
				assert(rowsOccupied[newPartition[loc]] > 0);
				assert(rowsOccupied[newPartition[loc]] <= end-start);
				rowsOccupied[newPartition[loc]]--;
				newPartition[loc]++;
				rowsOccupied[newPartition[loc]]++;
				assert(loc < end-start);
				assert(newPartition[loc] <= k);
				assert(rowsOccupied[newPartition[loc]] > 0);
				assert(rowsOccupied[newPartition[loc]] <= end-start);
			}
			loc = end-start-1;
		}
	}
	return ret;
}

//basically <
bool partitionCompare(const SolidPartition& left, const SolidPartition& right, size_t size)
{
#ifndef NDEBUG
	assert(size == left.assignments.actualSize);
	assert(size == right.assignments.actualSize);
#endif
	int compare = memcmp((const char*)left.assignments.data.data(), (const char*)right.assignments.data.data(), size*log2k/8);
	if (compare < 0)
	{
		return true;
	}
	if (compare > 0)
	{
		return false;
	}
	for (size_t i = (size*log2k/8)*8/log2k; i < size; i++)
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
std::vector<std::vector<size_t>> findExtensions(const std::vector<std::pair<size_t, SolidPartition>>& newLastRow, const std::vector<std::pair<size_t, SolidPartition>>& newNewRow, size_t size)
{
	std::vector<std::vector<size_t>> ret;
	ret.resize(newNewRow.size());
	size_t lastIndex = 0;
	size_t newIndex = 0;
	while (lastIndex < newLastRow.size() && newIndex < newNewRow.size())
	{
		if (partitionCompare(newLastRow[lastIndex].second, newNewRow[newIndex].second, size))
		{
			lastIndex++;
		}
		else if (partitionCompare(newNewRow[newIndex].second, newLastRow[lastIndex].second, size))
		{
			newIndex++;
		}
		else
		{
			size_t newIndexEnd = newIndex+1;
			size_t lastIndexEnd = lastIndex+1;
			while (newIndexEnd < newNewRow.size() && !partitionCompare(newNewRow[newIndexEnd-1].second, newNewRow[newIndexEnd].second, size))
			{
				newIndexEnd++;
			}
			while (lastIndexEnd < newLastRow.size() && !partitionCompare(newLastRow[lastIndexEnd-1].second, newLastRow[lastIndexEnd].second, size))
			{
				lastIndexEnd++;
			}
			for (size_t i = newIndex; i < newIndexEnd; i++)
			{
				for (size_t j = lastIndex; j < lastIndexEnd; j++)
				{
					assert(i < newNewRow.size());
					assert(j < newLastRow.size());
					ret[newNewRow[i].first].push_back(newLastRow[j].first);
				}
			}
			newIndex = newIndexEnd;
			lastIndex = lastIndexEnd;
		}
	}
	return ret;
}

std::vector<size_t> findExtensions(const std::vector<std::pair<size_t, SolidPartition>>& newLastRow, const std::vector<std::pair<size_t, SolidPartition>>& newNewRow, size_t size, const std::vector<double>& oldCosts)
{
	std::vector<size_t> ret;
	ret.resize(newNewRow.size(), -1);
	size_t lastIndex = 0;
	size_t newIndex = 0;
	while (lastIndex < newLastRow.size() && newIndex < newNewRow.size())
	{
		if (partitionCompare(newLastRow[lastIndex].second, newNewRow[newIndex].second, size))
		{
			lastIndex++;
		}
		else if (partitionCompare(newNewRow[newIndex].second, newLastRow[lastIndex].second, size))
		{
			newIndex++;
		}
		else
		{
			size_t newIndexEnd = newIndex+1;
			size_t lastIndexEnd = lastIndex+1;
			double minCost = oldCosts[newLastRow[lastIndex].first];
			size_t minCostIndex = newLastRow[lastIndex].first;
			while (lastIndexEnd < newLastRow.size() && !partitionCompare(newLastRow[lastIndexEnd-1].second, newLastRow[lastIndexEnd].second, size))
			{
				if (oldCosts[newLastRow[lastIndexEnd].first] < minCost)
				{
					minCost = oldCosts[newLastRow[lastIndexEnd].first];
					minCostIndex = newLastRow[lastIndexEnd].first;
				}
				lastIndexEnd++;
			}
			assert(newNewRow[newIndex].first < newNewRow.size());
			assert(ret[newNewRow[newIndex].first] == -1);
			ret[newNewRow[newIndex].first] = minCostIndex;
			while (newIndexEnd < newNewRow.size() && !partitionCompare(newNewRow[newIndexEnd-1].second, newNewRow[newIndexEnd].second, size))
			{
				assert(newNewRow[newIndexEnd].first < newNewRow.size());
				assert(ret[newNewRow[newIndexEnd].first] == -1);
				ret[newNewRow[newIndexEnd].first] = minCostIndex;
				newIndexEnd++;
			}
			newIndex = newIndexEnd;
			lastIndex = lastIndexEnd;
		}
	}
	assert(std::none_of(ret.begin(), ret.end(), [](size_t x) { return x == -1; }));
	return ret;
}

void SolidPartition::unpermutate(size_t size)
{
#ifndef NDEBUG
	assert(size == assignments.size());
#endif
	std::vector<size_t> mapping;
	mapping.resize(k, -1);
	size_t nextNum = 0;
	for (size_t i = 0; i < size; i++)
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
SolidPartition::SolidPartition(Iterator start, Iterator end, size_t numAssignments, TinyVectorMemoryAllocator& allocator) :
	assignments()
{
	assert(numAssignments > 0);
	assignments.reserve(numAssignments, allocator);
#ifndef NDEBUG
	assignments.actualSize = numAssignments;
#endif
	size_t pos = 0;
	while (start != end)
	{
		assert(pos < numAssignments);
		assignments[pos] = (size_t)*start;
		pos++;
		start++;
	}
	assert(pos == numAssignments);
}

std::vector<std::set<size_t>> getActiveRows(std::vector<SNPSupport> supports)
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

	std::vector<std::pair<size_t, size_t>> rowExtent;
	rowExtent.resize(maxRead, std::pair<size_t, size_t>{-1, 0});
	for (auto x : supports)
	{
		rowExtent[x.readNum].first = std::min(x.SNPnum, rowExtent[x.readNum].first);
		rowExtent[x.readNum].second = std::max(x.SNPnum, rowExtent[x.readNum].second);
	}

	std::vector<std::set<size_t>> ret;
	ret.resize(maxSNP);
	for (size_t i = 0; i < maxSNP; i++)
	{
		for (size_t a = 0; a < maxRead; a++)
		{
			if (rowExtent[a].first <= i && rowExtent[a].second >= i)
			{
				ret[i].insert(a);
			}
		}
	}
	return ret;
}

std::vector<size_t> findOptimalExtensions(const std::vector<std::vector<size_t>>& extensions, const std::vector<double>& oldCosts)
{
	std::vector<size_t> ret;
	ret.reserve(extensions.size());
	for (size_t i = 0; i < extensions.size(); i++)
	{
		assert(extensions[i].size() > 0);
		double minCost = oldCosts[extensions[i][0]];
		size_t index = extensions[i][0];
		for (size_t j = 1; j < extensions[i].size(); j++)
		{
			if (oldCosts[extensions[i][j]] < minCost)
			{
				minCost = oldCosts[extensions[i][j]];
				index = extensions[i][j];
			}
		}
		ret.push_back(index);
	}
	return ret;
}

std::vector<size_t> findOptimalExtensionsWithNoOverlap(size_t num, const std::vector<double>& costs)
{
	size_t index = 0;
	for (size_t i = 1; i < costs.size(); i++)
	{
		if (costs[i] < costs[index])
		{
			index = i;
		}
	}
	std::vector<size_t> ret;
	for (size_t i = 0; i < num; i++)
	{
		ret.push_back(index);
	}
	return ret;
}

std::vector<std::pair<size_t, SolidPartition>> splitIntersection(const std::vector<SparsePartition>& partitions, std::set<size_t> intersection, std::set<size_t> actives, TinyVectorMemoryAllocator& allocator)
{
	std::vector<std::pair<size_t, SolidPartition>> result;
	result.reserve(partitions.size());
	std::set<size_t> pickThese = subsetIndices(intersection, actives);
	size_t size = pickThese.size();
	for (size_t i = 0; i < partitions.size(); i++)
	{
		SolidPartition insertion = partitions[i].getSolidFromIndices(pickThese, allocator);
		insertion.unpermutate(size);
		result.emplace_back(i, insertion);
	}
	std::sort(result.begin(), result.end(), [size](const std::pair<size_t, SolidPartition>& left, const std::pair<size_t, SolidPartition>& right) { return partitionCompare(left.second, right.second, size); });
	return result;
}

template <typename T>
void clearVector(std::vector<T>& o)
{
	(std::vector<T>{}).swap(o);
}

//returns optimal partition and its score
std::tuple<std::vector<size_t>, double> haplotype(std::vector<SNPSupport> supports, size_t inK)
{
	k = inK;
	log2k = ceil(log2(k));
	std::cerr << "extents\n";
	size_t maxSNP = 0;
	size_t maxRead = 0;
	for (auto x : supports)
	{
		maxSNP = std::max(maxSNP, x.SNPnum);
		maxRead = std::max(maxRead, x.readNum);
	}
	maxSNP++;
	maxRead++;

	std::cerr << "split supports per SNP\n";
	std::vector<std::vector<SNPSupport>> supportsPerSNP;
	supportsPerSNP.resize(maxSNP);
	for (auto x : supports)
	{
		supportsPerSNP[x.SNPnum].push_back(x);
	}
	size_t firstSNP = 0;
	while (supportsPerSNP[firstSNP].size() == 0)
	{
		firstSNP++;
		assert(firstSNP < supportsPerSNP.size());
	}

	std::cerr << "active rows\n";
	std::vector<std::set<size_t>> activeRowsPerColumn = getActiveRows(supports);
	SparsePartitionContainer optimalPartitions { activeRowsPerColumn, k };

	std::cerr << "column " << firstSNP << " (" << activeRowsPerColumn[firstSNP].size() << ")";
	auto lastColumnTime = std::chrono::system_clock::now();
	Column oldColumn { supportsPerSNP[firstSNP].begin(), supportsPerSNP[firstSNP].end(), firstSNP, 0, maxRead };

	TinyVectorMemoryAllocator oldRowMemoryAllocator {activeRowsPerColumn[firstSNP].size(), supportsPerSNP[firstSNP].size(), k};

	std::vector<SparsePartition> oldRowPartitions = SparsePartition::getAllPartitions(activeRowsPerColumn[firstSNP], oldRowMemoryAllocator);
	std::vector<double> oldRowCosts;
	std::vector<size_t> oldOptimalPartitions;
	for (size_t i = 0; i < oldRowPartitions.size(); i++)
	{
		oldOptimalPartitions.push_back(optimalPartitions.insertPartition(oldRowPartitions[i], firstSNP, activeRowsPerColumn[firstSNP].size()));
	}
	for (auto x : oldRowPartitions)
	{
		oldRowCosts.push_back(x.deltaCost(oldColumn, activeRowsPerColumn[firstSNP]));
	}

	std::set<size_t> all = activeRowsPerColumn[firstSNP];

	for (size_t snp = firstSNP+1; snp < maxSNP; snp++)
	{
		auto newColumnTime = std::chrono::system_clock::now();
		auto diff = std::chrono::duration_cast<std::chrono::duration<int,std::milli>>(newColumnTime-lastColumnTime);
		std::set<size_t> intersect = setIntersection(activeRowsPerColumn[snp], activeRowsPerColumn[snp-1]);
		lastColumnTime = newColumnTime;
		std::cerr << " " << diff.count() << "ms\n";
		std::cerr << "column " << snp << " (" << activeRowsPerColumn[snp].size() << ", " << intersect.size() << ")";
		if (setEqual(activeRowsPerColumn[snp], activeRowsPerColumn[snp-1]))
		{
			Column col { supportsPerSNP[snp].begin(), supportsPerSNP[snp].end(), snp, 0, maxRead };
			for (size_t i = 0; i < oldRowCosts.size(); i++)
			{
				oldRowCosts[i] += oldRowPartitions[i].deltaCost(col, activeRowsPerColumn[snp]);
			}
			all = setUnion(all, activeRowsPerColumn[snp]);
			continue;
		}
		TinyVectorMemoryAllocator newRowMemoryAllocator { activeRowsPerColumn[snp].size(), activeRowsPerColumn[snp].size(), k };
		std::vector<SparsePartition> newRowPartitions = SparsePartition::getAllPartitions(activeRowsPerColumn[snp], newRowMemoryAllocator);
		std::cerr << " (" << newRowPartitions.size() << " partitions)";
		std::vector<size_t> optimalExtensions;
		if (intersect.size() == 0)
		{
			optimalExtensions = findOptimalExtensionsWithNoOverlap(newRowPartitions.size(), oldRowCosts);
		}
		else
		{
			TinyVectorMemoryAllocator tempOldRowMemoryAllocator { activeRowsPerColumn[snp-1].size(), intersect.size(), k };
			auto tempOldRowPartitions = splitIntersection(oldRowPartitions, intersect, activeRowsPerColumn[snp-1], tempOldRowMemoryAllocator);
			clearVector(oldRowPartitions);
			oldRowMemoryAllocator.empty();
			TinyVectorMemoryAllocator tempNewRowMemoryAllocator { activeRowsPerColumn[snp].size(), intersect.size(), k };
			auto tempNewRowPartitions = splitIntersection(newRowPartitions, intersect, activeRowsPerColumn[snp], tempNewRowMemoryAllocator);
			optimalExtensions = findExtensions(tempOldRowPartitions, tempNewRowPartitions, intersect.size(), oldRowCosts);
//			auto extensions = findExtensions(tempOldRowPartitions, tempNewRowPartitions, intersect.size());
//			auto optimalExtensions2 = findOptimalExtensions(extensions, oldRowCosts);
//			assert(std::equal(optimalExtensions.begin(), optimalExtensions.end(), optimalExtensions2.begin()));
		}
		Column col { supportsPerSNP[snp].begin(), supportsPerSNP[snp].end(), snp, 0, maxRead };
		std::vector<double> newRowCosts;
		std::vector<size_t> newOptimalPartitions;
		newRowCosts.reserve(newRowPartitions.size());
		newOptimalPartitions.reserve(newRowPartitions.size());
		optimalPartitions.reserveMore(optimalExtensions.size()*(activeRowsPerColumn[snp].size()-intersect.size()));
		for (size_t j = 0; j < optimalExtensions.size(); j++)
		{
			assert(optimalExtensions[j] < oldOptimalPartitions.size());
			assert(optimalExtensions[j] < oldRowCosts.size());
			newRowCosts.push_back(oldRowCosts[optimalExtensions[j]]+newRowPartitions[j].deltaCost(col, activeRowsPerColumn[snp]));
			newOptimalPartitions.push_back(optimalPartitions.extendPartition(oldOptimalPartitions[optimalExtensions[j]], newRowPartitions[j], snp, maxSNP, all, activeRowsPerColumn[snp], intersect));
		}
		clearVector(optimalExtensions);
		optimalPartitions.clearUnused(newOptimalPartitions);
		oldRowPartitions = std::move(newRowPartitions);
		oldRowCosts = std::move(newRowCosts);
		oldOptimalPartitions = std::move(newOptimalPartitions);
		oldRowMemoryAllocator = std::move(newRowMemoryAllocator);
		all = setUnion(all, activeRowsPerColumn[snp]);
	}

	size_t optimalResultIndex = 0;
	for (size_t i = 1; i < oldRowCosts.size(); i++)
	{
		if (oldRowCosts[i] < oldRowCosts[optimalResultIndex])
		{
			optimalResultIndex = i;
		}
	}

	TinyVectorMemoryAllocator allocator { 2, maxSNP, k };
	SolidPartition partition = optimalPartitions.getPartition(oldOptimalPartitions[optimalResultIndex], maxSNP, allocator).getSolid(all, allocator);
	double score = oldRowCosts[optimalResultIndex];

	std::vector<size_t> result;
	for (auto iter = partition.assignments.begin(); iter != partition.assignments.end(maxRead); iter++)
	{
		result.push_back((size_t)*iter);
	}
	return std::tuple<std::vector<size_t>, double> { result, score };
}
