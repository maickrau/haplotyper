//g++ renumberer.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o renumberer.exe
//./renumberer.exe outputResultFile inputResultFile inputRenumberingFile1 inputRenumberingFile2 ...

#include "variant_utils.h"

std::pair<std::vector<size_t>, size_t> loadResult(std::string fileName)
{
	std::ifstream file { fileName };
	std::pair<std::vector<size_t>, size_t> result;
	while (file.good())
	{
		size_t read;
		file >> read;
		if (file.good())
		{
			result.first.push_back(read);
		}
	}
	result.second = result.first.back();
	result.first.pop_back();
	return result;
}

void writeResult(std::vector<size_t> assignments, size_t error, std::string fileName)
{
	std::ofstream file { fileName };
	for (size_t i = 0; i < assignments.size(); i++)
	{
		file << assignments[i] << " ";
	}
	file << "\n" << error << "\n";
}

int main(int argc, char** argv)
{
	SupportRenumbering renumbering = loadRenumbering(argv[3]);
	for (int i = 4; i < argc; i++)
	{
		SupportRenumbering add = loadRenumbering(argv[i]);
		renumbering = renumbering.merge(add);
	}
	std::pair<std::vector<size_t>, size_t> result = loadResult(argv[2]);
	std::vector<size_t> actualResult;
	for (size_t i = 0; i < renumbering.readSize(); i++)
	{
		if (!renumbering.hasReadRenumbering(i))
		{
			actualResult.push_back(-1);
		}
		else
		{
			actualResult.push_back(result.first[renumbering.getReadRenumbering(i)]);
		}
	}
	writeResult(actualResult, result.second, argv[1]);
}
