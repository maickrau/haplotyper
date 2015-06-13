//g++ matrix_bander.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o matrix_bander.exe
//./matrix_bander.exe inputSupportsFile barycentricIterations outputSupportsFile renumberingsFile temperature temperatureMultiplier annealingIterationr

#include <iostream>
#include <map>
#include <algorithm>
#include <random>
#include <chrono>
#include <cmath>
#include <cassert>

#include "variant_utils.h"

class MovedSupport
{
public:
	MovedSupport(const MovedSupport& second) = default;
	MovedSupport(SNPSupport s) : newRowNum(s.readNum), newColNum(s.SNPnum), support(s) {};
	size_t newRowNum;
	size_t newColNum;
	SNPSupport support;
};

std::vector<MovedSupport> transpose(const std::vector<MovedSupport>& old)
{
	std::vector<MovedSupport> result;
	for (auto x : old)
	{
		result.emplace_back(x);
		result.back().newRowNum = x.newColNum;
		result.back().newColNum = x.newRowNum;
	}
	return result;
}

std::vector<MovedSupport> barycentricSort(const std::vector<MovedSupport>& old)
{
	std::vector<std::pair<double, std::vector<MovedSupport>>> rows;

	for (auto x : old)
	{
		if (rows.size() <= x.newRowNum)
		{
			rows.resize(x.newRowNum+1);
		}
		rows[x.newRowNum].second.emplace_back(x);
	}

	for (auto& x : rows)
	{
		double sum = 0;
		for (auto y : x.second)
		{
			sum += y.newColNum;
		}
		sum /= (double)x.second.size();
		x.first = sum;
	}

	std::sort(rows.begin(), rows.end(), [](const std::pair<double, std::vector<MovedSupport>>& left, const std::pair<double, std::vector<MovedSupport>>& right) { return left.first < right.first; });

	std::vector<MovedSupport> ret;
	for (size_t i = 0; i < rows.size(); i++)
	{
		for (auto x : rows[i].second)
		{
			ret.emplace_back(x);
			ret.back().newRowNum = i;
		}
	}

	return ret;
}

std::vector<SNPSupport> getMovedSupports(const std::vector<MovedSupport>& supports)
{
	std::vector<SNPSupport> ret;
	for (auto x : supports)
	{
		ret.push_back(x.support);
		ret.back().readNum = x.newRowNum;
		ret.back().SNPnum = x.newColNum;
	}
	return ret;
}

double getEnergy(std::vector<SNPSupport> supports)
{
	std::vector<size_t> minColumn;
	std::vector<size_t> maxColumn;
	for (auto x : supports)
	{
		if (minColumn.size() <= x.readNum)
		{
			minColumn.resize(x.readNum+1, -1);
		}
		if (maxColumn.size() <= x.readNum)
		{
			maxColumn.resize(x.readNum+1, 0);
		}
		minColumn[x.readNum] = std::min(minColumn[x.readNum], x.SNPnum);
		maxColumn[x.readNum] = std::max(maxColumn[x.readNum], x.SNPnum);
	}
	for (size_t i = 1; i < maxColumn.size(); i++)
	{
		maxColumn[i] = std::max(maxColumn[i], maxColumn[i-1]);
	}
	for (size_t i = minColumn.size()-2; i < minColumn.size(); i--)
	{
		minColumn[i] = std::min(minColumn[i], minColumn[i+1]);
	}
	double ret = 0;
	for (size_t i = 0; i < minColumn.size(); i++)
	{
		ret += pow(2, maxColumn[i]-minColumn[i]);
	}
	return ret;
}

double getEnergy(const std::vector<MovedSupport>& supports)
{
	return getEnergy(getMovedSupports(supports));
}

SupportRenumbering getSupportRenumbering(const std::vector<MovedSupport>& supports)
{
	SupportRenumbering ret;
	for (auto x : supports)
	{
		ret.addReadRenumbering(x.support.readNum, x.newRowNum);
		ret.addSNPRenumbering(x.support.SNPnum, x.newColNum);
	}
	return ret;
}

SupportRenumbering makeBandedBarycentric(const std::vector<SNPSupport>& supports, size_t iterations)
{
	std::vector<MovedSupport> locations;
	for (auto x : supports)
	{
		locations.emplace_back(x);
	}

	std::vector<MovedSupport> best { locations };
	double bestEnergy = getEnergy(best);

	for (size_t i = 0; i < iterations; i++)
	{
		//sort/transpose twice so result won't be transposed
		locations = barycentricSort(locations);
		locations = transpose(locations);
		locations = barycentricSort(locations);
		locations = transpose(locations);

		double newEnergy = getEnergy(locations);
		if (newEnergy < bestEnergy)
		{
			std::cerr << "iteration " << i << " new best " << newEnergy << "\n";
			best = locations;
			bestEnergy = newEnergy;
			assert(getSupportRenumbering(best).checkValidity());
		}
	}
	return getSupportRenumbering(best);
}

SupportRenumbering transpose(const SupportRenumbering& renumbering)
{
	assert(renumbering.checkValidity());
	SupportRenumbering ret;
	for (size_t i = 0; i < renumbering.readSize(); i++)
	{
		ret.addSNPRenumbering(i, renumbering.getReadRenumbering(i));
	}
	for (size_t i = 0; i < renumbering.SNPSize(); i++)
	{
		ret.addReadRenumbering(i, renumbering.getSNPRenumbering(i));
	}
	assert(ret.checkValidity());
	return ret;
}

SupportRenumbering swapK(const SupportRenumbering& numbering, int k)
{
	SupportRenumbering ret { numbering };

	std::mt19937 mt {(size_t)std::chrono::system_clock::now().time_since_epoch().count()};
	std::uniform_int_distribution<size_t> rowSelector {0, numbering.readSize()-1};
	std::uniform_int_distribution<size_t> columnSelector {0, numbering.SNPSize()-1};
	for (int i = 0; i < k; i++)
	{
		size_t firstRow = rowSelector(mt);
		size_t secondRow = rowSelector(mt);
		ret.swapRows(firstRow, secondRow);

		size_t firstColumn = columnSelector(mt);
		size_t secondColumn = columnSelector(mt);
		ret.swapColumns(firstColumn, secondColumn);
	}

	return ret;
}

SupportRenumbering swapOneRow(const SupportRenumbering& numbering)
{
	SupportRenumbering ret { numbering };

	std::mt19937 mt {(size_t)std::chrono::system_clock::now().time_since_epoch().count()};
	std::uniform_int_distribution<size_t> rowSelector {0, numbering.readSize()-1};
	size_t firstRow = rowSelector(mt);
	size_t secondRow = rowSelector(mt);
	ret.swapRows(firstRow, secondRow);

	return ret;
}

SupportRenumbering swapOneColumn(const SupportRenumbering& numbering)
{
	SupportRenumbering ret { numbering };

	std::mt19937 mt {(size_t)std::chrono::system_clock::now().time_since_epoch().count()};
	std::uniform_int_distribution<size_t> columnSelector {0, numbering.SNPSize()-1};
	size_t firstColumn = columnSelector(mt);
	size_t secondColumn = columnSelector(mt);
	ret.swapColumns(firstColumn, secondColumn);

	return ret;
}

SupportRenumbering adjSwapK(const SupportRenumbering& numbering, int k)
{
	SupportRenumbering ret { numbering };

	std::mt19937 mt {(size_t)std::chrono::system_clock::now().time_since_epoch().count()};
	std::uniform_int_distribution<size_t> rowSelector {0, numbering.readSize()-1};
	std::uniform_int_distribution<size_t> columnSelector {0, numbering.SNPSize()-1};
	for (int i = 0; i < k; i++)
	{
		size_t firstRow = rowSelector(mt);
		ret.swapRows(firstRow, (firstRow+1) % numbering.readSize());

		size_t firstColumn = columnSelector(mt);
		ret.swapColumns(firstColumn, (firstColumn+1) % numbering.SNPSize());
	}

	return ret;
}

SupportRenumbering adjSwapOneRow(const SupportRenumbering& numbering)
{
	SupportRenumbering ret { numbering };

	std::mt19937 mt {(size_t)std::chrono::system_clock::now().time_since_epoch().count()};
	std::uniform_int_distribution<size_t> rowSelector {0, numbering.readSize()-1};
	size_t firstRow = rowSelector(mt);
	ret.swapRows(firstRow, (firstRow+1) % numbering.readSize());

	return ret;
}

SupportRenumbering adjSwapOneColumn(const SupportRenumbering& numbering)
{
	SupportRenumbering ret { numbering };

	std::mt19937 mt {(size_t)std::chrono::system_clock::now().time_since_epoch().count()};
	std::uniform_int_distribution<size_t> columnSelector {0, numbering.SNPSize()-1};
	size_t firstColumn = columnSelector(mt);
	ret.swapColumns(firstColumn, (firstColumn+1) % numbering.SNPSize());

	return ret;
}

SupportRenumbering reverseRows(const SupportRenumbering& numbering)
{
	std::vector<size_t> rows;
	for (size_t i = 0; i < numbering.readSize(); i++)
	{
		rows.push_back(numbering.getReadRenumbering(i));
	}

	std::mt19937 mt {(size_t)std::chrono::system_clock::now().time_since_epoch().count()};
	std::uniform_int_distribution<size_t> rowSelector {0, numbering.readSize()-1};

	size_t firstRow = rowSelector(mt);
	size_t secondRow;
	do
	{
		secondRow = rowSelector(mt);
	} while (secondRow == firstRow);
	size_t lengthDiff;
	if (firstRow < secondRow)
	{
		lengthDiff = secondRow-firstRow;
	}
	else
	{
		lengthDiff = numbering.readSize()+secondRow-firstRow;
	}

	for (size_t i = 0; i < lengthDiff/2; i++)
	{
		std::swap(rows[(firstRow+i) % rows.size()], rows[(secondRow+rows.size()-i) % rows.size()]);
	}

	SupportRenumbering ret;
	for (size_t i = 0; i < rows.size(); i++)
	{
		ret.addReadRenumbering(i, rows[i]);
	}
	for (size_t i = 0; i < numbering.SNPSize(); i++)
	{
		ret.addSNPRenumbering(i, numbering.getSNPRenumbering(i));
	}
	return ret;
}

SupportRenumbering reverse(const SupportRenumbering& numbering)
{
	SupportRenumbering ret = reverseRows(numbering);
	assert(ret.checkValidity());
	ret = transpose(ret);
	ret = reverseRows(ret);
	assert(ret.checkValidity());
	ret = transpose(ret);
	return ret;
}

SupportRenumbering relocateRows(const SupportRenumbering& numbering)
{
	std::vector<size_t> rows;
	for (size_t i = 0; i < numbering.readSize(); i++)
	{
		rows.push_back(numbering.getReadRenumbering(i));
	}

	std::mt19937 mt {(size_t)std::chrono::system_clock::now().time_since_epoch().count()};
	std::uniform_int_distribution<size_t> rowSelector {0, numbering.readSize()-1};

	size_t firstRow = rowSelector(mt);
	size_t secondRow;
	do
	{
		secondRow = rowSelector(mt);
	} while (secondRow == firstRow);
	size_t k;
	do
	{
		k = rowSelector(mt);
	} while (k == 0);

	std::vector<size_t> tempK;
	for (size_t i = 0; i < k; i++)
	{
		tempK.push_back(rows[(secondRow+i) % rows.size()]);
	}
	for (size_t i = secondRow+rows.size()+k-1; i % rows.size() != (firstRow+k)%rows.size(); i--)
	{
		rows[i % rows.size()] = rows[(i+rows.size()-k) % rows.size()];
	}
	rows[(firstRow+k) % rows.size()] = rows[firstRow];
	for (size_t i = firstRow; i < firstRow + k; i++)
	{
		rows[i % rows.size()] = tempK[i-firstRow];
	}

	SupportRenumbering ret;
	for (size_t i = 0; i < rows.size(); i++)
	{
		ret.addReadRenumbering(i, rows[i]);
	}
	for (size_t i = 0; i < numbering.SNPSize(); i++)
	{
		ret.addSNPRenumbering(i, numbering.getSNPRenumbering(i));
	}
	return ret;
}

SupportRenumbering relocate(const SupportRenumbering& numbering)
{
	SupportRenumbering ret = relocateRows(numbering);
	assert(ret.checkValidity());
	ret = transpose(ret);
	ret = relocateRows(ret);
	assert(ret.checkValidity());
	ret = transpose(ret);
	return ret;
}

SupportRenumbering getNeighbor(const SupportRenumbering& numbering)
{
    std::mt19937 mt {(size_t)std::chrono::system_clock::now().time_since_epoch().count()};
    std::uniform_int_distribution<int> method{0, 3};
    int chosen = method(mt);
    SupportRenumbering ret {numbering};

    switch (chosen)
    {
	case 0:
		ret = swapOneRow(numbering);
		break;
	case 1:
		ret = swapOneColumn(numbering);
		break;
	case 2:
		ret = adjSwapOneRow(numbering);
		break;
	case 3:
		ret = adjSwapOneColumn(numbering);
		break;
/*	case 0:
		ret = swapK(numbering, 1);
		assert(ret.checkValidity());
		break;
	case 1:
		ret = adjSwapK(numbering, 1);
		assert(ret.checkValidity());
		break;
	case 2:
		ret = reverse(numbering);
		assert(ret.checkValidity());
		break;
	case 3:
		ret = relocate(numbering);
		assert(ret.checkValidity());
		break;*/
    }
    
    return ret;
}

SupportRenumbering getRandomRenumbering(const std::vector<SNPSupport>& supports)
{
	size_t maxRow = supports[0].readNum;
	size_t maxSNP = supports[0].SNPnum;
	for (auto x : supports)
	{
		maxRow = std::max(maxRow, x.readNum);
		maxSNP = std::max(maxSNP, x.SNPnum);
	}
	std::vector<size_t> rowPermutation;
	std::vector<size_t> SNPPermutation;
	for (size_t i = 0; i <= maxRow; i++)
	{
		rowPermutation.push_back(i);
	}
	for (size_t i = 0; i <= maxSNP; i++)
	{
		SNPPermutation.push_back(i);
	}
    std::mt19937 mt {(size_t)std::chrono::system_clock::now().time_since_epoch().count()};

	std::shuffle(rowPermutation.begin(), rowPermutation.end(), mt);
	std::shuffle(SNPPermutation.begin(), SNPPermutation.end(), mt);

	SupportRenumbering ret;
	for (size_t i = 0; i < rowPermutation.size(); i++)
	{
		ret.addReadRenumbering(i, rowPermutation[i]);
	}
	for (size_t i = 0; i < SNPPermutation.size(); i++)
	{
		ret.addSNPRenumbering(i, SNPPermutation[i]);
	}
	return ret;
}

SupportRenumbering makeBandedSimulatedAnnealing(const std::vector<SNPSupport>& supports, SupportRenumbering start, int iterations, double temperature, double temperatureMultiplier)
{
	SupportRenumbering best = start;
	double bestEnergy = getEnergy(renumberSupports(supports, best));
	SupportRenumbering current = best;
	double currentEnergy = bestEnergy;

    std::mt19937 mt {(size_t)std::chrono::system_clock::now().time_since_epoch().count()};
	std::uniform_real_distribution<double> changeCurrent {0, 1};

	for (int i = 0; i < iterations; i++)
	{
		SupportRenumbering newRenumbering = getNeighbor(current);
		double newEnergy = getEnergy(renumberSupports(supports, newRenumbering));
		if (newEnergy < bestEnergy)
		{
			best = newRenumbering;
			bestEnergy = newEnergy;
			std::cerr << "iteration " << i << " new best " << bestEnergy << "\n";
		}
		if (changeCurrent(mt) < std::min(1.0, exp((currentEnergy-newEnergy)/temperature)))
		{
			current = newRenumbering;
			currentEnergy = newEnergy;
		}
		temperature *= temperatureMultiplier;
	}
	return best;
}

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);

	std::cerr << "barycentric banding\n";
	SupportRenumbering firstRenumbering = makeBandedBarycentric(supports, std::stol(argv[2]));
	assert(firstRenumbering.checkValidity());

	std::cerr << "annealing banding\n";
	SupportRenumbering secondRenumbering = makeBandedSimulatedAnnealing(supports, firstRenumbering, std::stoi(argv[7]), std::stod(argv[5]), std::stod(argv[6]));
	assert(secondRenumbering.checkValidity());

	std::cerr << "writing output\n";
	SupportRenumbering renumbering = secondRenumbering;
	std::vector<SNPSupport> result = renumberSupports(supports, renumbering);
	assert(renumbering.checkValidity());
	writeSupports(result, argv[3]);
	writeRenumbering(renumbering, argv[4]);
}