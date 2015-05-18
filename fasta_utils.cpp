#include "fasta_utils.h"

Genome loadFasta(std::string fileName)
{
	return loadFastas(fileName)[0];
}

std::vector<Genome> loadFastas(std::string fileName)
{
	std::ifstream file { fileName };
	std::vector<Genome> ret;
	char a = file.get();
	while (file.good())
	{
		if (a == '>')
		{
			ret.emplace_back();
			std::getline(file, ret.back().name);
			a = file.get();
		}
		if (a != '\n')
		{
			ret.back().bases.push_back(a);
		}
		a = file.get();
	}
	return ret;
}

void writeFasta(Genome genome, std::string fileName)
{
	std::ofstream file { fileName };
	file << ">" << genome.name << "\n";
	size_t loc = 0;
	while (loc+80 < genome.bases.size())
	{
		file << genome.bases.substr(loc, 80) << "\n";
		loc += 80;
	}
	file << genome.bases.substr(loc) << "\n";
}
