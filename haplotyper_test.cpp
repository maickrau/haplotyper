//g++ haplotyper_test.cpp haplotyper.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o haplotyper_test.exe
//./haplotyperTest.exe

#include <cassert>

#include "haplotyper.h"

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports {
		{1, 1, 'A', 1},
		{1, 2, 'A', 1},
		{1, 3, 'G', 1},
		{2, 1, 'T', 1},
		{2, 2, 'T', 1},
		{2, 3, 'T', 1},
		{3, 2, 'C', 1},
		{3, 3, 'A', 1},
		{3, 4, 'A', 1},
		{4, 3, 'T', 1},
		{4, 4, 'T', 1}
	};
	assert(haplotype(supports, 1).second == 6);
	assert(haplotype(supports, 2).second == 2);
	assert(haplotype(supports, 3).second == 0);
	supports[2].variant = 'A';
	assert(haplotype(supports, 1).second == 6);
	assert(haplotype(supports, 2).second == 1);
	assert(haplotype(supports, 3).second == 0);
	supports[7].variant = 'G';
	assert(haplotype(supports, 1).second == 6);
	assert(haplotype(supports, 2).second == 2);
	assert(haplotype(supports, 3).second == 0);
	supports[6].variant = 'A';
	assert(haplotype(supports, 1).second == 5);
	assert(haplotype(supports, 2).second == 1);
	assert(haplotype(supports, 3).second == 0);
}