//./make_reads.exe genome.fasta reads.fasta numOfReads lengthOfReads
//g++ make_reads.cpp fasta_utils.cpp -o make_reads.exe -std=c++11

#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <chrono>

#include "fasta_utils.h"

std::vector<Genome> createReads(Genome genome, int num, int length)
{
    std::mt19937 mt(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<size_t> dist(0, genome.bases.size()-length);
    std::vector<Genome> ret;
    for (int i = 0; i < num; i++)
    {
    	size_t startPos = dist(mt);
    	Genome read;
    	read.bases = genome.bases.substr(startPos, length);
    	read.name = genome.name+" "+std::to_string(startPos)+" "+std::to_string(i);
    	ret.push_back(read);
    }
    return ret;
}

int main(int argc, char** argv)
{
	Genome genome = loadFasta(argv[1]);
	int numOfReads = std::stoi(argv[3]);
	int lengthOfReads = std::stoi(argv[4]);
	auto reads = createReads(genome, numOfReads, lengthOfReads);
	writeFasta(reads.begin(), reads.end(), argv[2]);
}