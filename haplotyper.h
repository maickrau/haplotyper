#ifndef haplotyper_h
#define haplotyper_h

#include <tuple>
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

class TinyVectorMemoryAllocator
{
public:
	TinyVectorMemoryAllocator(size_t coverage, size_t length, size_t k);
	~TinyVectorMemoryAllocator();
	TinyVectorMemoryAllocator(const TinyVectorMemoryAllocator& second) = delete;
	TinyVectorMemoryAllocator& operator=(const TinyVectorMemoryAllocator& second) = delete;
	TinyVectorMemoryAllocator(TinyVectorMemoryAllocator&& second);
	TinyVectorMemoryAllocator& operator=(TinyVectorMemoryAllocator&& second);
	unsigned char* allocate(size_t size);
	void empty();
private:
	unsigned char* memory;
//	std::vector<unsigned char> memory;
	size_t bytesPerVector;
	size_t nextEmpty;
	size_t size;
};

template <typename T>
class TinyVector
{
public:
	TinyVector() : contents(nullptr), contentsSize(0) {};
	~TinyVector() {}
	TinyVector(const TinyVector& second) : contents(nullptr), contentsSize(0)
	{
		*this = second;
	}
	TinyVector& operator=(const TinyVector& second)
	{
		if (this == &second)
		{
			return *this;
		}
/*		if (contents != nullptr)
		{
			delete [] contents;
		}
		contents = new T[second.contentsSize];
		for (int i = 0; i < second.contentsSize; i++)
		{
			contents[i] = second.contents[i];
		}*/
		contents = second.contents;
		contentsSize = second.contentsSize;
		return *this;
	}
	T* begin()
	{
		assert(contents != nullptr);
		return contents;
	}
	T* end()
	{
		assert(contents != nullptr);
		return contents+contentsSize;
	}
	const T* begin() const
	{
		assert(contents != nullptr);
		return contents;
	}
	const T* end() const
	{
		assert(contents != nullptr);
		return contents+contentsSize;
	}
	size_t size() const
	{
		assert(contents != nullptr);
		return contentsSize;
	}
	size_t capacity() const
	{
		assert(contents != nullptr);
		return contentsSize;
	}
	void resize(unsigned char newCapacity, const T& defaultValue, TinyVectorMemoryAllocator& allocator)
	{
//		T* newData = new T[newCapacity];
		T* newData = allocator.allocate(newCapacity);
		if (contents != nullptr)
		{
			for (int i = 0; i < contentsSize && i < newCapacity; i++)
			{
				newData[i] = contents[i];
			}
			for (int i = contentsSize; i < newCapacity; i++)
			{
				newData[i] = defaultValue;
			}
//			delete [] contents;
		}
		else
		{
			for (int i = 0; i < newCapacity; i++)
			{
				newData[i] = defaultValue;
			}
		}
		contents = newData;
		contentsSize = newCapacity;
	}
	T* data()
	{
		assert(contents != nullptr);
		return contents;
	}
	const T* data() const
	{
		assert(contents != nullptr);
		return contents;
	}
	T& operator[](size_t index)
	{
		assert(contents != nullptr);
		assert(index < contentsSize);
		return contents[index];
	}
	const T& operator[](size_t index) const
	{
		assert(contents != nullptr);
		assert(index < contentsSize);
		return contents[index];
	}
private:
	T* contents;
	unsigned char contentsSize;
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
		size_t getValue() const;
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
			assert(pos < container.size());
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
	PartitionAssignments();
	PartitionAssignmentElement operator[](size_t pos);
	PartitionAssignmentElementConst operator[](size_t pos) const;
#ifndef NDEBUG
	size_t size() const;
	size_t capacity() const;
	size_t actualSize;
#endif
	void reserve(size_t numAssignments, TinyVectorMemoryAllocator& allocator);
	iterator<PartitionAssignmentElement> begin();
	iterator<PartitionAssignmentElement> end(size_t size);
	iterator<PartitionAssignmentElementConst> begin() const;
	iterator<PartitionAssignmentElementConst> end(size_t size) const;

private:
	void extendCapacity(size_t newCapacity, size_t defaultValue, TinyVectorMemoryAllocator& allocator);
	TinyVector<unsigned char> data;
	friend class PartitionAssignments::PartitionAssignmentElement;
	friend bool partitionCompare(const SolidPartition& left, const SolidPartition& right, size_t size);
};

class SolidPartition
{
public:
	SolidPartition();
	static std::vector<SolidPartition> getAllPartitions(size_t start, size_t end, TinyVectorMemoryAllocator& allocator);
	template <typename Iterator>
	SolidPartition(Iterator start, Iterator end, size_t numAssignments, TinyVectorMemoryAllocator& allocator);
	void unpermutate(size_t size);
	size_t getk() const;

	PartitionAssignments assignments;
private:
};

class SparsePartition
{
public:
	SparsePartition();
	SparsePartition(SolidPartition inner);
	static std::vector<SparsePartition> getAllPartitions(std::set<size_t> actives, TinyVectorMemoryAllocator& allocator);
	SolidPartition getSolid(const std::set<size_t>& actives, TinyVectorMemoryAllocator& allocator) const;
	SolidPartition getSolidFromIndices(const std::set<size_t>& pickThese, TinyVectorMemoryAllocator& allocator) const;
	double deltaCost(const Column& col, const std::set<size_t>& actives) const;
	size_t getk() const;
	size_t getAssignment(size_t loc, const std::set<size_t>& actives) const;

	SolidPartition inner;
private:
	SolidPartition getSubset(const std::set<size_t>& subset, const std::set<size_t>& actives, TinyVectorMemoryAllocator& allocator) const;
};

std::tuple<std::vector<size_t>, double> haplotype(std::vector<SNPSupport> supports, size_t k);

#endif