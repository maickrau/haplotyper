#ifndef variant_utils_h
#define variant_utils_h

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
	SupportRenumbering();
	SupportRenumbering reverse();
	SupportRenumbering merge(SupportRenumbering second);
	void addReadRenumbering(size_t oldRead, size_t newRead);
	void addSNPRenumbering(size_t oldSNP, size_t newSNP);
	void swapRows(size_t firstRow, size_t secondRow);
	void swapColumns(size_t firstColumn, size_t secondColumn);
	size_t getReadRenumbering(size_t oldRead) const;
	size_t getSNPRenumbering(size_t oldSNP) const;
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

#endif