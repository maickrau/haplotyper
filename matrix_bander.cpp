//g++ matrix_bander.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o matrix_bander.exe
//./matrix_bander.exe inputSupportsFile barycentricIterations outputSupportsFile renumberingsFile temperature temperatureMultiplier annealingIterationr rowDistanceIterations rowDistanceCutoff

#include <iostream>
#include <map>
#include <algorithm>
#include <random>
#include <chrono>
#include <cmath>
#include <cassert>
#include <functional>

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
	std::vector<std::pair<size_t, size_t>> rowExtents;
	size_t maxRead = 0;
	size_t maxSNP = 0;
	for (auto x : supports)
	{
		maxRead = std::max(maxRead, x.readNum);
		maxSNP = std::max(maxSNP, x.SNPnum);
	}
	maxRead++;
	maxSNP++;
	rowExtents.resize(maxRead, {-1, 0});
	for (auto x : supports)
	{
		rowExtents[x.readNum].first = std::min(x.SNPnum, rowExtents[x.readNum].first);
		rowExtents[x.readNum].second = std::max(x.SNPnum, rowExtents[x.readNum].second);
	}
	double total = 0;
	for (size_t i = 0; i < maxSNP; i++)
	{
		size_t numRows = 0;
		for (size_t j = 0; j < rowExtents.size(); j++)
		{
			if (rowExtents[j].first <= i && rowExtents[j].second >= i)
			{
				numRows++;
			}
		}
		total += pow(2, numRows);
	}
	return total;
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

SupportRenumbering getIdentityRenumbering(const std::vector<SNPSupport>& supports)
{
	SupportRenumbering ret;
	for (auto x : supports)
	{
		ret.addReadRenumbering(x.readNum, x.readNum);
		ret.addSNPRenumbering(x.SNPnum, x.SNPnum);
	}
	return ret;
}

template <typename RowFunction>
SupportRenumbering makeBandedAlternating(const std::vector<SNPSupport>& supports, SupportRenumbering initialRenumbering, size_t iterations, RowFunction rowOrderer)
{
	std::vector<MovedSupport> locations;
	
	{
		std::vector<SNPSupport> renumbered = renumberSupports(supports, initialRenumbering);
		for (auto x : renumbered)
		{
			locations.emplace_back(x);
		}
	}

	std::vector<MovedSupport> best { locations };
	double bestEnergy = getEnergy(best);

	for (size_t i = 0; i < iterations; i++)
	{
		//sort/transpose twice so result won't be transposed
		locations = rowOrderer(locations);
		locations = transpose(locations);
		locations = rowOrderer(locations);
		locations = transpose(locations);

		assert(getSupportRenumbering(locations).checkValidity());
		double newEnergy = getEnergy(locations);
		if (newEnergy < bestEnergy)
		{
			std::cerr << "iteration " << i << " new best " << newEnergy << "\n";
			best = locations;
			bestEnergy = newEnergy;
			assert(getSupportRenumbering(best).checkValidity());
		}
	}
	SupportRenumbering result = getSupportRenumbering(best);
	return initialRenumbering.merge(result);
}

SupportRenumbering makeBandedBarycentric(const std::vector<SNPSupport>& supports, SupportRenumbering numbering, size_t iterations)
{
	return makeBandedAlternating(supports, numbering, iterations, barycentricSort);
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

SupportRenumbering shrinkBiggestColumn(const SupportRenumbering& numbering, const std::vector<SNPSupport>& supports)
{
	std::vector<size_t> minRow;
	std::vector<size_t> maxRow;
	minRow.resize(numbering.SNPSize(), -1);
	maxRow.resize(numbering.SNPSize(), 0);
	for (auto x : supports)
	{
		size_t readNum = numbering.getReadRenumbering(x.readNum);
		size_t SNPnum = numbering.getSNPRenumbering(x.SNPnum);
		minRow[SNPnum] = std::min(minRow[SNPnum], readNum);
		maxRow[SNPnum] = std::max(maxRow[SNPnum], readNum);
	}
	size_t biggestSize = maxRow[0]-minRow[0];
	size_t biggestIndex = 0;
	for (size_t i = 1; i < maxRow.size(); i++)
	{
		if (maxRow[i]-minRow[i] > biggestSize)
		{
			biggestSize = maxRow[i]-minRow[i];
			biggestIndex = i;
		}
	}
	SupportRenumbering ret { numbering };
	ret.swapRows(minRow[biggestIndex], minRow[biggestIndex]+1);
	ret.swapRows(maxRow[biggestIndex], maxRow[biggestIndex]-1);
	return ret;
}

SupportRenumbering permutateColumns(const SupportRenumbering& numbering, std::set<size_t> permutableIndices)
{
	std::mt19937 mt {(size_t)std::chrono::system_clock::now().time_since_epoch().count()};
	std::uniform_real_distribution<double> continueSelector {0, 1.0};
	std::uniform_int_distribution<size_t> rowSelector {0, numbering.SNPSize()-1};

	while (permutableIndices.size() < 2 || continueSelector(mt) < 0.8)
	{
		permutableIndices.insert(rowSelector(mt));
	}

	std::vector<size_t> indices;
	for (auto x : permutableIndices)
	{
		indices.push_back(x);
	}
	std::shuffle(indices.begin(), indices.end(), mt);

	std::vector<size_t> columnIndices;
	columnIndices.reserve(numbering.SNPSize());
	size_t num = 0;
	for (size_t i = 0; i < numbering.SNPSize(); i++)
	{
		if (permutableIndices.count(i) == 0)
		{
			columnIndices.push_back(numbering.getSNPRenumbering(i));
		}
		else
		{
			assert(num < indices.size());
			columnIndices.push_back(numbering.getSNPRenumbering(indices[num]));
			num++;
		}
	}
	assert(num == indices.size());

	SupportRenumbering ret;
	for (size_t i = 0; i < numbering.SNPSize(); i++)
	{
		ret.addSNPRenumbering(i, columnIndices[i]);
	}
	for (size_t i = 0; i < numbering.readSize(); i++)
	{
		ret.addReadRenumbering(i, numbering.getReadRenumbering(i));
	}

	return ret;
}

SupportRenumbering permutateColumns(const SupportRenumbering& numbering)
{
	return permutateColumns(numbering, std::set<size_t>{});
}

SupportRenumbering permutateColumnsGuided(const SupportRenumbering& numbering, const std::vector<SNPSupport>& supports)
{
	std::vector<std::pair<size_t, size_t>> rowExtents;
	rowExtents.resize(numbering.readSize(), {-1, 0});
	for (auto x : supports)
	{
		rowExtents[numbering.getReadRenumbering(x.readNum)].first = std::max(rowExtents[numbering.getReadRenumbering(x.readNum)].first, numbering.getSNPRenumbering(x.SNPnum));
		rowExtents[numbering.getReadRenumbering(x.readNum)].second = std::min(rowExtents[numbering.getReadRenumbering(x.readNum)].second, numbering.getSNPRenumbering(x.SNPnum));
	}
	size_t maxSNPIndex = 0;
	size_t maxSNPCoverage = 0;
	for (size_t i = 0; i < numbering.SNPSize(); i++)
	{
		size_t coverage = 0;
		for (size_t j = 0; j < numbering.readSize(); j++)
		{
			if (rowExtents[j].first <= i && rowExtents[j].second >= i)
			{
				coverage++;
			}
		}
		if (coverage > maxSNPCoverage)
		{
			maxSNPIndex = i;
			maxSNPCoverage = coverage;
		}
	}
	return permutateColumns(numbering, std::set<size_t> { maxSNPIndex });
}

SupportRenumbering getNeighbor(const SupportRenumbering& numbering, const std::vector<SNPSupport>& supports)
{
    std::mt19937 mt {(size_t)std::chrono::system_clock::now().time_since_epoch().count()};
    std::uniform_int_distribution<int> method{0, 1};
    int chosen = method(mt);
    SupportRenumbering ret;
    switch (chosen)
    {
	case 0:
		ret = permutateColumns(numbering);
		break;
	case 1:
		ret = permutateColumnsGuided(numbering, supports);
		break;
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
		SupportRenumbering newRenumbering = getNeighbor(current, supports);
		assert(newRenumbering.checkValidity());
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

size_t rowDistance(const std::vector<MovedSupport>& supports, size_t firstRowNum, size_t secondRowNum)
{
	std::vector<bool> firstRow;
	std::vector<bool> secondRow;
	for (auto x : supports)
	{
		if (x.newRowNum == firstRowNum || x.newRowNum == secondRowNum)
		{
			if (firstRow.size() <= x.newColNum)
			{
				firstRow.resize(x.newColNum+1, false);
				secondRow.resize(x.newColNum+1, false);
			}
		}
		if (x.newRowNum == firstRowNum)
		{
			firstRow[x.newColNum] = true;
		}
		if (x.newRowNum == secondRowNum)
		{
			secondRow[x.newColNum] = true;
		}
	}
	size_t result = 0;
	for (size_t i = 0; i < firstRow.size(); i++)
	{
		if (firstRow[i] != secondRow[i])
		{
			result++;
		}
	}
	return result;
}

size_t rowLeftness(const std::vector<MovedSupport>& supports, size_t rowIndex)
{
	size_t ret = -1;
	for (auto x : supports)
	{
		if (x.newRowNum == rowIndex)
		{
			ret = std::min(ret, x.newColNum);
		}
	}
	return ret;
}

std::vector<MovedSupport> greedyRowDistanceSorter(const std::vector<MovedSupport>& old, int breakpointDistance)
{
	size_t rows = old[0].newRowNum;
	for (auto x : old)
	{
		rows = std::max(rows, x.newRowNum);
	}
	rows += 1;

	std::vector<size_t> breakpoints;
	breakpoints.push_back(0);
	for (size_t a = 1; a < rows; a++)
	{
		if (rowDistance(old, a-1, a) >= breakpointDistance)
		{
			breakpoints.push_back(a);
		}
	}
	if (breakpoints.size() < 2)
	{
		return old;
	}
	std::vector<bool> breakpointUsed;
	std::vector<size_t> breakpointOrdering;

	size_t leftmostBlockIndex = 0;
	size_t leftmostBlockLeftness = rowLeftness(old, 0);
	for (size_t i = 0; i < breakpoints.size(); i++)
	{
		size_t newBlockLeftness = rowLeftness(old, breakpoints[i]);
		if (newBlockLeftness < leftmostBlockLeftness)
		{
			leftmostBlockLeftness = newBlockLeftness;
			leftmostBlockIndex = i;
		}
	}
	breakpointOrdering.push_back(leftmostBlockIndex);

	breakpoints.push_back(rows);
	breakpointUsed.resize(breakpoints.size(), false);
	breakpointUsed[breakpointOrdering[0]] = true;
	for (size_t iteration = 0; iteration < breakpoints.size()-2; iteration++)
	{
		size_t bestDistance = 0;
		size_t bestDistanceIndex = 0;
		for (size_t i = 0; i < breakpoints.size()-1; i++)
		{
			if (!breakpointUsed[i])
			{
				size_t distance = rowDistance(old, breakpoints[breakpointOrdering.back()+1]-1, breakpoints[i]);
				if (distance < bestDistance || bestDistanceIndex == 0)
				{
					bestDistanceIndex = i;
					bestDistance = distance;
				}
			}
		}
		breakpointUsed[bestDistanceIndex] = true;
		breakpointOrdering.push_back(bestDistanceIndex);
	}

	std::vector<MovedSupport> ret;
	size_t currentRow = 0;
	for (size_t i = 0; i < breakpointOrdering.size(); i++)
	{
		for (auto x : old)
		{
			if (x.newRowNum >= breakpoints[breakpointOrdering[i]] && x.newRowNum < breakpoints[breakpointOrdering[i]+1])
			{
				ret.emplace_back(x);
				ret.back().newRowNum = currentRow+x.newRowNum-breakpoints[breakpointOrdering[i]];
			}
		}
		currentRow += breakpoints[breakpointOrdering[i]+1]-breakpoints[breakpointOrdering[i]];
	}
	return ret;
}

SupportRenumbering makeBandedRowDistance(const std::vector<SNPSupport>& supports, SupportRenumbering numbering, size_t iterations, size_t rowDistance)
{
	return makeBandedAlternating(supports, numbering, iterations, std::bind(greedyRowDistanceSorter, std::placeholders::_1, rowDistance));
}

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);

	std::cerr << "starting score " << getEnergy(supports) << "\n";

	SupportRenumbering numbering = getIdentityRenumbering(supports);

	std::cerr << "row distance banding\n";
	numbering = makeBandedRowDistance(supports, numbering, std::stol(argv[8]), std::stol(argv[9]));
	assert(numbering.checkValidity());

	std::cerr << "barycentric banding\n";
	numbering = makeBandedBarycentric(supports, numbering, std::stol(argv[2]));
	assert(numbering.checkValidity());

	std::cerr << "annealing banding\n";
	numbering = makeBandedSimulatedAnnealing(supports, numbering, std::stoi(argv[7]), std::stod(argv[5]), std::stod(argv[6]));
	assert(numbering.checkValidity());

	std::cerr << "writing output\n";
	SupportRenumbering renumbering = numbering;
	std::vector<SNPSupport> result = renumberSupports(supports, renumbering);
	assert(renumbering.checkValidity());
	writeSupports(result, argv[3]);
	writeRenumbering(renumbering, argv[4]);
}