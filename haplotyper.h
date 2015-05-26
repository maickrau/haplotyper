#ifndef haplotyper_h
#define haplotyper_h

#include <cstdlib>
#include <vector>

#include "variant_utils.h"

class Column
{
public:
	template <typename Iterator>
	Column(Iterator start, Iterator end, size_t column, size_t minRow, size_t maxRow) : variants(), minRow(minRow), maxRow(maxRow)
	{
		static_assert(std::is_constructible<SNPSupport, decltype(*start)>::value, "");
		variants.resize(maxRow+1, 0);
		costs.resize(maxRow+1);
		while (start != end)
		{
			SNPSupport s = *start;
			if (s.SNPnum == column)
			{
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

class Partition
{
public:
	Partition();
	static std::vector<Partition> getAllPartitions(size_t start, size_t end, size_t k);
	template <typename Iterator>
	Partition(Iterator start, Iterator end);
	Partition filter(size_t newMinRow, size_t newMaxRow);
	bool extends(Partition second);
	void unpermutate();
	double deltaCost(Column col);
	Partition merge(Partition second);
	std::vector<size_t> assignments;
	size_t minRow;
	size_t maxRow;
	size_t k;
private:
	double wCost(Column col, char variant, size_t haplotype);
};

std::pair<Partition, double> haplotype(std::vector<SNPSupport> supports, size_t k);

#endif