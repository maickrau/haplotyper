#ifndef variant_utils_h
#define variant_utils_h

#include <cassert>
#include <string>
#include <set>

#include "fasta_utils.h"

class SNPPosition
{
public:
	SNPPosition(Genome gen1, Genome gen2);
	std::string leftFlank;
	std::string rightFlank;
	std::set<char> variants;
};

class SupportRenumbering
{
public:
	static SupportRenumbering identity(size_t maxRead, size_t maxSNP);
	SupportRenumbering();
	SupportRenumbering reverse();
	SupportRenumbering merge(SupportRenumbering second);
	void overwriteReadRenumbering(size_t oldRead, size_t newRead);
	void addReadRenumbering(size_t oldRead, size_t newRead);
	void addSNPRenumbering(size_t oldSNP, size_t newSNP);
	void swapRows(size_t firstRow, size_t secondRow);
	void swapColumns(size_t firstColumn, size_t secondColumn);
	size_t getReadRenumbering(size_t oldRead) const;
	size_t getSNPRenumbering(size_t oldSNP) const;
	bool hasReadRenumbering(size_t oldRead) const;
	bool hasSNPRenumbering(size_t oldSNP) const;
	size_t readSize() const;
	size_t SNPSize() const;
	bool checkValidity() const;
private:
	std::vector<size_t> readRenumbering;
	std::vector<size_t> SNPRenumbering;
	friend SupportRenumbering loadRenumbering(std::string fileName);
	friend void writeRenumbering(SupportRenumbering renumbering, std::string fileName);
};

class SNPSupport
{
public:
	SNPSupport(size_t readNum, size_t SNPnum, char variant);
	SNPSupport(size_t readNum, size_t SNPnum, char variant, double support);
	size_t readNum;
	size_t SNPnum;
	double support;
	char variant;
};

std::vector<SNPSupport> loadSupports(std::string fileName);
void writeSupports(std::vector<SNPSupport> supports, std::string fileName);
std::vector<SNPSupport> renumberSupports(std::vector<SNPSupport> supports, SupportRenumbering renumbering);
void writeRenumbering(SupportRenumbering renumbering, std::string fileName);
SupportRenumbering loadRenumbering(std::string fileName);


class SNPLine
{
public:
	template <typename Iterator>
	SNPLine(Iterator supportsStart, Iterator supportsEnd, size_t readNum) : variantsAtLocations(), supportsAtLocations(), readNum(readNum)
	{
		static_assert(std::is_constructible<SNPSupport, decltype(*supportsStart)>::value, "");
		size_t lastSNPnum = 0;
		while (supportsStart != supportsEnd)
		{
			SNPSupport x = *supportsStart;
			if (x.readNum == readNum)
			{
				assert(x.SNPnum >= lastSNPnum);
				variantsAtLocations.emplace_back(x.SNPnum, x.variant);
				supportsAtLocations.push_back(x.support);
				lastSNPnum = x.SNPnum;
			}
			supportsStart++;
		}
	}
	SNPLine();
	std::vector<std::pair<size_t, char>> variantsAtLocations;
	std::vector<double> supportsAtLocations;
	void merge(SNPLine second);
	size_t readNum;
	bool operator==(const SNPLine& second) const;
	bool operator!=(const SNPLine& second) const;
	bool contains(const SNPLine& second) const;
	char variantAt(size_t loc) const;
	void mergeSubset(SNPLine second);
	std::vector<SNPSupport> toSupports() const;
};

std::vector<SNPLine> makeLines(std::vector<SNPSupport> supports);
std::vector<SNPSupport> mergeRows(std::vector<SNPSupport> oldSupports, size_t row1, size_t row2);

#endif