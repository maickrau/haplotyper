#ifndef haplotyper_h
#define haplotyper_h

#include <cstdlib>
#include <vector>
#include <set>
#include <cassert>

#include "variant_utils.h"

class Column
{
public:
	template <typename Iterator>
	Column(Iterator start, Iterator end, size_t column, size_t minRow, size_t maxRow) : variants(), minRow(minRow), maxRow(maxRow)
	{
		static_assert(std::is_constructible<SNPSupport, decltype(*start)>::value, "");
		variants.resize(maxRow+1, 0);
		costs.resize(maxRow+1, 0);
		while (start != end)
		{
			SNPSupport s = *start;
			if (s.SNPnum == column)
			{
				assert(s.readNum <= maxRow);
				assert(s.variant == 'A' || s.variant == 'T' || s.variant == 'C' || s.variant == 'G');
				variants[s.readNum] = s.variant;
				costs[s.readNum] = s.support;
			}
			start++;
		}
	};
	std::vector<char> variants;
	std::vector<double> costs;
	size_t minRow;
	size_t maxRow;
};

class SolidPartition;

class PartitionAssignments
{
public:
	class PartitionAssignmentElement
	{
	public:
		operator size_t() const;
		PartitionAssignmentElement& operator=(size_t value);
		PartitionAssignmentElement& operator=(const PartitionAssignmentElement) = delete;
	private:
		PartitionAssignmentElement(PartitionAssignments& container, size_t pos);
		PartitionAssignments& container;
		size_t pos;
		friend class PartitionAssignments;
	};
	class PartitionAssignmentElementConst
	{
	public:
		operator size_t() const;
	private:
		PartitionAssignmentElementConst(const PartitionAssignments& container, size_t pos);
		PartitionAssignmentElementConst(PartitionAssignmentElement el);
		PartitionAssignmentElement el;
		friend class PartitionAssignments;
	};
	template <typename ElementType>
	class iterator
	{
	public:
		typedef size_t difference_type;
		typedef ElementType value_type;
		typedef ElementType& reference;
		typedef ElementType* pointer;
		typedef std::bidirectional_iterator_tag iterator_category;
		ElementType operator*() const
		{
			return container[pos];
		}
		iterator(const iterator& second) = default;
		iterator& operator=(const iterator& second) = default;
		iterator operator++()
		{
			pos++;
			return *this;
		}
		iterator operator++(int)
		{
			PartitionAssignments::iterator<ElementType> ret = *this;
			pos++;
			return ret;
		}
		iterator operator--()
		{
			pos--;
			return *this;
		}
		iterator operator--(int)
		{
			PartitionAssignments::iterator<ElementType> ret = *this;
			pos--;
			return ret;
		}
		iterator operator+(size_t add)
		{
			pos += add;
			return *this;
		}
		iterator operator-(size_t deduct)
		{
			pos -= deduct;
			return *this;
		}
		bool operator==(const iterator& second) const
		{
			return &container == &second.container && pos == second.pos;
		}
		bool operator!=(const iterator& second) const
		{
			return !(*this == second);
		}
	private:
		iterator(PartitionAssignments& container, size_t pos) :
			container(container),
			pos(pos)
		{
		}

		PartitionAssignments& container;
		size_t pos;
		friend class PartitionAssignments;
	};
	PartitionAssignments(size_t k);
	PartitionAssignmentElement operator[](size_t pos);
	PartitionAssignmentElementConst operator[](size_t pos) const;
	void push_back(size_t value);
	size_t size() const;
	size_t capacity() const;
	iterator<PartitionAssignmentElement> begin();
	iterator<PartitionAssignmentElement> end();
	iterator<PartitionAssignmentElementConst> begin() const;
	iterator<PartitionAssignmentElementConst> end() const;
private:
	void extendCapacity(size_t newCapacity, size_t defaultValue);
	size_t k;
	size_t log2k;
	size_t actualSize;
	std::vector<unsigned char> data;
	friend class PartitionAssignments::PartitionAssignmentElement;
	friend bool partitionCompare(const SolidPartition& left, const SolidPartition& right);
};

class SolidPartition
{
public:
	SolidPartition(size_t k);
	static std::vector<SolidPartition> getAllPartitions(size_t start, size_t end, size_t k);
	template <typename Iterator>
	SolidPartition(Iterator start, Iterator end, size_t k);
	bool extends(const SolidPartition& second) const;
	void unpermutate();
	size_t getk() const;
	size_t getValue(size_t index) const;

	PartitionAssignments assignments;
	size_t minRow;
	size_t maxRow;
private:
	size_t k;
};

class SparsePartition
{
public:
	SparsePartition(size_t k);
	SparsePartition(SolidPartition inner, std::set<size_t> actives);
	static std::vector<SparsePartition> getAllPartitions(std::set<size_t> actives, size_t k);
	SparsePartition merge(const SparsePartition& second) const;
	SolidPartition getComparableIntersection(const SparsePartition& second) const;
	SolidPartition getComparableIntersection(const std::set<size_t>& newActives) const;
	SolidPartition getSolid() const;
	bool extends(const SparsePartition& second) const;
	double deltaCost(const Column& col) const;
	size_t getk() const;

	SolidPartition inner;
	std::set<size_t> actives;
private:
	size_t getAssignment(size_t loc) const;
	SolidPartition getSubset(const std::set<size_t>& subset) const;
	size_t k;
};

std::pair<SolidPartition, double> haplotype(std::vector<SNPSupport> supports, size_t k);

#endif