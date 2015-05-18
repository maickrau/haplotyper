//./find_snps.exe genomesFile resultsFile
//g++ find_snps.cpp -std=c++11 -o find_snps.exe

#include <fstream>
#include <set>

#include "fasta_utils.cpp"

int main(int argc, char** argv)
{
	std::ofstream result { argv[2] };
	std::vector<Genome> genomes = loadFastas(argv[1]);
	for (size_t i = 0; i < genomes[0].bases.size(); i++)
	{
		std::set<char> variants;
		for (size_t a = 1; a < genomes.size(); a++)
		{
			if (genomes[a].bases[i] != genomes[a-1].bases[i])
			{
				variants.insert(genomes[a].bases[i]);
				variants.insert(genomes[a-1].bases[i]);
			}
		}
		if (variants.size() > 0)
		{
			result << i;
			for (char x : variants)
			{
				result << " " << x;
			}
			result << "\n";
		}
	}
}