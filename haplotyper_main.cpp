//g++ haplotyper_main.cpp haplotyper.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o haplotyper_main.exe
//./haplotyper_main.exe supportsFile k

#include <iostream>

#include "haplotyper.h"

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);
	auto result = haplotype(supports, std::stoi(argv[2]));
	for (auto x = std::get<0>(result).begin(); x != std::get<0>(result).end(); x++)
	{
		std::cout << *x << " ";
	}
	std::cout << "\n";
	std::cout << std::get<1>(result);
	std::cout << "\n";
}