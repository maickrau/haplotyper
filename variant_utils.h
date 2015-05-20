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

#endif