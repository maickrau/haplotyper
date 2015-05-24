//g++ reposition_reads.cpp fasta_utils.cpp -o reposition_reads.exe -std=c++11

#include <sstream>
#include <algorithm>

#include "fasta_utils.h"

size_t getReadPos(Genome read)
{
	std::istringstream s { read.name };
	std::string unused;
	size_t pos;
	s >> unused >> pos;
	return pos;
}

int main(int argc, char** argv)
{
	std::vector<Genome> reads = loadFastas(argv[1]);
	std::sort(reads.begin(), reads.end(), [](Genome left, Genome right) { return getReadPos(left) < getReadPos(right); });
	writeFasta(reads.begin(), reads.end(), argv[2]);
}