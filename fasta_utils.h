#ifndef fasta_utils_h
#define fasta_utils_h

#include <string>
#include <fstream>
#include <vector>

class Genome
{
public:
	Genome(std::string name, std::string bases) : name(name), bases(bases) {};
	Genome() : name(), bases() {};
	std::string name;
	std::string bases;
};

Genome loadFasta(std::string fileName);
std::vector<Genome> loadFastas(std::string fileName);
void writeFasta(Genome genome, std::string fileName);

template <typename GenomeIterator>
void writeFasta(GenomeIterator genomeStart, GenomeIterator genomeEnd, std::string fileName)
{
	static_assert(std::is_constructible<Genome, decltype(*genomeStart)>::value, "");
	std::ofstream file {fileName};
	while (genomeStart != genomeEnd)
	{
		file << ">" << genomeStart->name << "\n";
		size_t loc = 0;
		std::string genome = genomeStart->bases;
		while (loc+80 < genome.size())
		{
			file << genome.substr(loc, 80) << "\n";
			loc += 80;
		}
		file << genome.substr(loc) << "\n";

		genomeStart++;
	}
}

#endif