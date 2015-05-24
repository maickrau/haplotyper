#ifndef haplotyper_h
#define haplotyper_h

#include <cstdlib>
#include <vector>

#include "variant_utils.h"

class Column
{
public:
	template <typename Iterator>
	Column(Iterator start, Iterator end, size_t column) : variants(), minRow(-1), maxRow(0)
	{
		static_assert(std::is_constructible<SNPSupport, decltype(*start)>::value, "");
		while (start != end)
		{
			SNPSupport s = *start;
			if (s.SNPnum == column)
			{
				if (variants.size() <= s.readNum)
				{
					variants.resize(s.readNum+1, 0);
				}
				variants[s.readNum] = s.variant;
				if (costs.size() <= s.readNum)
				{
					costs.resize(s.readNum+1, 0);
				}
				costs[s.readNum] = s.support;
				minRow = std::min(minRow, s.readNum);
				maxRow = std::max(maxRow, s.readNum+1);
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
	static std::vector<Partition> getAllPartitions(size_t start, size_t end, size_t k);
	template <typename Iterator>
	Partition(Iterator start, Iterator end);
	Partition filter(size_t newMinRow, size_t newMaxRow);
	bool extends(Partition second);
	void unpermutate();
	std::vector<size_t> assignments;
	double deltaCost(Column col);
	Partition merge(Partition second);
	size_t minRow;
	size_t maxRow;
	size_t k;
private:
	double wCost(Column col, char variant, size_t haplotype);
};

std::pair<Partition, double> haplotype(std::vector<SNPSupport> supports, size_t k);

#endif