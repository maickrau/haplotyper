//g++ switch_distance_scorer.cpp fasta_utils.cpp -std=c++11 -o switch_distance_scorer.exe
//./switch_distance_scorer.exe haplotypingResultsFile allReadsFile allGenomesFile maxSubstitutions k

#include <algorithm>
#include <tuple>
#include <iostream>
#include <sstream>
#include <map>
#include <cassert>

#include "fasta_utils.h"

std::vector<size_t> findSNPIndices(std::vector<Genome> genomes)
{
	std::vector<size_t> ret;
	for (size_t i = 0; i < genomes[0].bases.size(); i++)
	{
		for (size_t a = 1; a < genomes.size(); a++)
		{
			if (genomes[a].bases[i] != genomes[a-1].bases[i])
			{
				ret.push_back(i);
				continue;
			}
		}
	}
	return ret;
}

std::vector<char> getNucleotidesAtSNPPositions(Genome genome, std::vector<size_t> SNPindices)
{
	std::vector<char> ret;
	for (size_t i = 0; i < SNPindices.size(); i++)
	{
		ret.push_back(genome.bases[SNPindices[i]]);
	}
	return ret;
}

Genome getConsensusFromReads(const std::vector<Genome>& reads, size_t actualSize)
{
	std::vector<std::map<char, size_t>> consensusSupport;
	consensusSupport.resize(actualSize);
	for (auto x : reads)
	{
		std::istringstream s { x.name };
		std::string originalGenome;
		s >> originalGenome;
		size_t location;
		s >> location;
		assert(location >= 0);
		assert(location+x.bases.size() <= actualSize);
		for (size_t i = 0; i < x.bases.size(); i++)
		{
			consensusSupport[location+i][x.bases[i]]++;
		}
	}
	Genome ret;
	ret.bases.resize(actualSize);
	for (size_t i = 0; i < actualSize; i++)
	{
		size_t maxValue = 0;
		for (auto x : consensusSupport[i])
		{
			maxValue = std::max(maxValue, x.second);
		}
		for (auto x : consensusSupport[i])
		{
			if (maxValue == x.second)
			{
				ret.bases[i] = x.first;
			}
		}
	}
	return ret;
}

std::vector<int> loadAssignments(std::string fileName, int k)
{
	std::vector<int> ret;
	std::ifstream file {fileName};
	while (file.good())
	{
		char a = file.peek();
		if (a == '\n')
		{
			break;
		}
		while (a == ' ')
		{
			file.get();
			a = file.peek();
		}
		if (a == '-')
		{
			file.get();
			ret.push_back(-1);
		}
		else
		{
			int assignment;
			file >> assignment;
			ret.push_back(assignment);
		}
	}
	ret.pop_back();
	return ret;
}

std::vector<std::vector<Genome>> splitReadsByAssignments(const std::vector<Genome>& reads, std::vector<int> assignments, int k)
{
	std::vector<std::vector<Genome>> ret;
	ret.resize(k);
	for (int i = 0; i < assignments.size(); i++)
	{
		if (assignments[i] != -1)
		{
			assert(assignments[i] >= 0);
			assert(assignments[i] < k);
			ret[assignments[i]].push_back(reads[i]);
		}
	}
	return ret;
}

std::vector<std::vector<int>> getAllPossibleAssignments(int k)
{
	std::vector<std::vector<int>> ret;
	std::vector<int> thisAssignment;
	for (int i = 0; i < k; i++)
	{
		thisAssignment.push_back(0);
	}
	while (ret.size() < pow(k, k))
	{
		ret.push_back(thisAssignment);
		for (int i = k-1; i >= 0; i--)
		{
			thisAssignment[i]++;
			if (thisAssignment[i] == k)
			{
				thisAssignment[i] = 0;
			}
			else
			{
				break;
			}
		}
	}
	return ret;
}

std::vector<std::vector<int>> getPermutations(int k)
{
	std::vector<std::vector<int>> ret;
	std::vector<int> thisPermutation;
	for (int i = 0; i < k; i++)
	{
		thisPermutation.push_back(i);
	}
	bool more = true;
	do
	{
		ret.push_back(thisPermutation);
		more = std::next_permutation(thisPermutation.begin(), thisPermutation.end());
	} while (more);
	return ret;
}

int permutationDifference(const std::vector<int>& p1, const std::vector<int>& p2)
{
	int ret = 0;
	for (size_t i = 0; i < p1.size(); i++)
	{
		if (p1[i] != p2[i])
		{
			ret += 1;
		}
	}
	return ret;
}

int substitutions(const std::vector<std::vector<char>>& genomes, const std::vector<std::vector<char>>& consensuses, const std::vector<int>& permutation, size_t index)
{
	int result = 0;
	for (size_t i = 0; i < permutation.size(); i++)
	{
		if (consensuses[i][index] > 0 && genomes[permutation[i]][index] != consensuses[i][index])
		{
			result++;
		}
	}
	return result;
}

//switch distance, substitutions
std::tuple<int, int> getScore(std::vector<std::vector<char>> genomes, std::vector<std::vector<char>> consensuses, int maxSubstitutions, size_t k)
{
//	std::vector<std::vector<int>> permutations = getPermutations(k);
	std::vector<std::vector<int>> permutations = getAllPossibleAssignments(k);

	//permutation index, substitutions, score
	std::vector<std::tuple<size_t, int, int>> score;
	size_t length = genomes[0].size();
	for (size_t i = 0; i < permutations.size(); i++)
	{
		score.emplace_back(i, 0, 0);
	}
	for (size_t i = 0; i < length; i++)
	{
		std::cerr << i << "/" << length << " (" << score.size() << ")\n";
		std::vector<std::tuple<size_t, int, int>> newScores;
		for (size_t j = 0; j < permutations.size(); j++)
		{
			std::vector<int> bestScores;
			bestScores.resize(maxSubstitutions, 0);
			std::vector<bool> hasScores;
			hasScores.resize(maxSubstitutions, false);
			for (size_t k = 0; k < score.size(); k++)
			{
				int newSubstitutions = std::get<1>(score[k])+substitutions(genomes, consensuses, permutations[std::get<0>(score[k])], i);
				if (newSubstitutions < maxSubstitutions)
				{
					int newScore = std::get<2>(score[k])+permutationDifference(permutations[j], permutations[std::get<0>(score[k])]);
					if (!hasScores[newSubstitutions] || newScore < bestScores[newSubstitutions])
					{
						bestScores[newSubstitutions] = newScore;
						hasScores[newSubstitutions] = true;
					}
				}
			}
			for (int k = 0; k < maxSubstitutions; k++)
			{
				if (hasScores[k])
				{
					newScores.emplace_back(j, k, bestScores[k]);
				}
			}
		}
		score = newScores;
	}
	auto bestScore = *std::min_element(score.begin(), score.end(), [](std::tuple<size_t, int, int> left, std::tuple<size_t, int, int> right) {return std::get<2>(left) < std::get<2>(right);});

	return std::make_tuple(std::get<2>(bestScore), std::get<1>(bestScore));
}

int main(int argc, char** argv)
{
	int maxSubstitutions = std::stoi(argv[4]);
	int k = std::stoi(argv[5]);
	std::cerr << "genomes\n";
	std::vector<Genome> genomes = loadFastas(argv[3]);
	std::cerr << "reads\n";
	std::vector<Genome> reads = loadFastas(argv[2]);
	std::cerr << "assignments\n";
	std::vector<int> assignments = loadAssignments(argv[1], k);
	std::cerr << "split reads\n";
	std::vector<std::vector<Genome>> splitReads = splitReadsByAssignments(reads, assignments, k);
	std::cerr << "snp positions\n";
	std::vector<size_t> SNPPositions = findSNPIndices(genomes);
	std::cerr << "genome nucleotides\n";
	std::vector<std::vector<char>> genomeNucleotides;
	for (auto x : genomes)
	{
		genomeNucleotides.push_back(getNucleotidesAtSNPPositions(x, SNPPositions));
	}
	std::cerr << "consensus nucleotides\n";
	std::vector<std::vector<char>> consensusNucleotides;
	for (auto x : splitReads)
	{
		Genome consensus = getConsensusFromReads(x, genomes[0].bases.size());
		consensusNucleotides.push_back(getNucleotidesAtSNPPositions(consensus, SNPPositions));
	}
	for (auto x : genomeNucleotides)
	{
		for (auto y : x)
		{
			std::cout << y;
		}
		std::cout << "\n";
	}
	std::cout << "\n";
	for (auto x : consensusNucleotides)
	{
		for (auto y : x)
		{
			std::cout << y;
		}
		std::cout << "\n";
	}
	std::cerr << "scoring\n";
	auto score = getScore(genomeNucleotides, consensusNucleotides, maxSubstitutions, k);
	std::cout << std::get<0>(score) << " " << std::get<1>(score) << "\n";
}