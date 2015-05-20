//./brute_force_mapper.exe allReadsFile snpsFile outputFile
//g++ brute_force_mapper.cpp fasta_utils.cpp variant_utils.cpp -std=c++11 -o brute_force_mapper.exe

#include <cstring>
#include <set>
#include <algorithm>
#include <vector>
#include <iostream>

#include "fasta_utils.h"
#include "variant_utils.h"

std::vector<SNPPosition> parseSNPs(std::vector<Genome> SNPcontigs)
{
	std::vector<SNPPosition> ret;
	for (size_t i = 0; i < SNPcontigs.size(); i += 2)
	{
		ret.emplace_back(SNPcontigs[i], SNPcontigs[i+1]);
	}
	return ret;
}

bool stringMatch(std::string left, std::string right)
{
	return strncmp(left.c_str(), right.c_str(), std::min(left.size(), right.size())) == 0;
}

//Since the reads are error-free, the reads must match the SNP contigs completely.
//The SNP/read pair is interesting only if the read contains the SNP.
//So, find pairs where the read contains the SNP, and both flanks match exactly until the SNP contig or read ends.
//returns: either the variant in the read, or 0 for no match
char overlaps(Genome read, SNPPosition SNP)
{
	for (size_t a = 0; a < read.bases.size(); a++)
	{
		if (SNP.variants.count(read.bases[a]) == 1)
		{
			if (stringMatch(SNP.rightFlank, read.bases.substr(a+1)))
			{
				std::string reverseReadLeftpart = read.bases.substr(0, a);
				std::reverse(reverseReadLeftpart.begin(), reverseReadLeftpart.end());
				if (stringMatch(SNP.leftFlank, reverseReadLeftpart))
				{
					return read.bases[a];
				}
			}
		}
	}
	return 0;
}

std::vector<SNPSupport> matchReadsToSNPs(std::vector<Genome> reads, std::vector<SNPPosition> SNPs)
{
	std::vector<SNPSupport> ret;
	for (size_t SNPIndex = 0; SNPIndex < SNPs.size(); SNPIndex++)
	{
		for (size_t readIndex = 0; readIndex < reads.size(); readIndex++)
		{
			char overlap = overlaps(reads[readIndex], SNPs[SNPIndex]);
			if (overlap != 0)
			{
				ret.emplace_back(readIndex, SNPIndex, overlap);
			}
		}
	}
	return ret;
}

std::vector<Genome> toUpperCase(std::vector<Genome> in)
{
	std::vector<Genome> out;
	for (auto x : in)
	{
		std::string result;
		result.resize(x.bases.size());
		std::transform(x.bases.begin(), x.bases.end(), result.begin(), toupper);
		out.emplace_back(x.name, result);
	}
	return out;
}

int main(int argc, char** argv)
{
	std::cerr << "reads ";
	std::vector<Genome> reads = loadFastas(argv[1]);
	std::cerr << reads.size() << "\n";
	std::cerr << "contigs ";
	std::vector<Genome> SNPcontigs = loadFastas(argv[2]);
	std::cerr << SNPcontigs.size() << "\n";
	std::cerr << "uppercase ";
	SNPcontigs = toUpperCase(SNPcontigs);
	std::cerr << SNPcontigs.size() << "\n";
	std::cerr << "parse ";
	std::vector<SNPPosition> SNPs = parseSNPs(SNPcontigs);
	std::cerr << SNPs.size() << "\n";
	std::cerr << "match ";
	std::vector<SNPSupport> readSNPMatching = matchReadsToSNPs(reads, SNPs);
	std::cerr << readSNPMatching.size() << "\n";
	std::cerr << "sort ";
	std::sort(readSNPMatching.begin(), readSNPMatching.end(), [](SNPSupport left, SNPSupport right) { return left.readNum < right.readNum || (left.readNum == right.readNum && left.SNPnum < right.SNPnum); });
	std::cerr << readSNPMatching.size() << "\n";
	std::cerr << "print ";
	writeSupports(readSNPMatching, argv[3]);
}