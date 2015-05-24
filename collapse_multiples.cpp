//./collapse_multiples.exe inputFile outputFile renumberingFile
//g++ collapse_multiples.cpp variant_utils.cpp -std=c++11 -o collapse_multiples.exe

#include <cassert>
#include <algorithm>
#include <iostream>
#include <unordered_map>

#include "variant_utils.h"

class SNPLine
{
public:
	template <typename Iterator>
	SNPLine(Iterator supportsStart, Iterator supportsEnd, size_t readNum) : variantsAtLocations(), supportsAtLocations(), readNum(readNum)
	{
		static_assert(std::is_constructible<SNPSupport, decltype(*supportsStart)>::value, "");
		while (supportsStart != supportsEnd)
		{
			SNPSupport x = *supportsStart;
			if (x.readNum == readNum)
			{
				variantsAtLocations.emplace_back(x.SNPnum, x.variant);
				supportsAtLocations.push_back(x.support);
			}
			supportsStart++;
		}
	}
	SNPLine() {};
	std::vector<std::pair<size_t, char>> variantsAtLocations;
	std::vector<double> supportsAtLocations;
	void merge(SNPLine second);
	size_t readNum;
	bool operator==(const SNPLine& second)
	{
		return variantsAtLocations == second.variantsAtLocations;
	}
	bool operator!=(const SNPLine& second)
	{
		return !(*this == second);
	}
	std::vector<SNPSupport> toSupports();
};

std::vector<SNPSupport> SNPLine::toSupports()
{
	std::vector<SNPSupport> ret;
	for (size_t i = 0; i < variantsAtLocations.size(); i++)
	{
		ret.emplace_back(readNum, variantsAtLocations[i].first, variantsAtLocations[i].second, supportsAtLocations[i]);
	}
	return ret;
}

void SNPLine::merge(SNPLine second)
{
	assert(supportsAtLocations.size() == second.supportsAtLocations.size());
	for (size_t i = 0; i < supportsAtLocations.size(); i++)
	{
		supportsAtLocations[i] += second.supportsAtLocations[i];
	}
}

std::vector<SNPSupport> mergeSupports(std::vector<SNPSupport> supports)
{
	std::cout << "start\n";
	std::sort(supports.begin(), supports.end(), [](SNPSupport left, SNPSupport right) { return left.SNPnum < right.SNPnum; });
	std::stable_sort(supports.begin(), supports.end(), [](SNPSupport left, SNPSupport right) { return left.readNum < right.readNum; });
	std::vector<SNPLine> rows;
	size_t maxReadNum = 0;
	for (auto x : supports)
	{
		maxReadNum = std::max(maxReadNum, x.readNum);
	}
	rows.resize(maxReadNum+1);
	std::cout << "make lines\n";
	size_t lastReadStart = 0;
	for (size_t i = 1; i < supports.size(); i++)
	{
		if (supports[i].readNum != supports[i-1].readNum)
		{
			rows[supports[lastReadStart].readNum] = SNPLine {supports.begin()+lastReadStart, supports.begin()+i, supports[lastReadStart].readNum};
			lastReadStart = i;
		}
	}
	rows[supports[lastReadStart].readNum] = SNPLine {supports.begin()+lastReadStart, supports.end(), supports[lastReadStart].readNum};
	std::cout << "made lines\n";
	size_t lastRead = 0;
	std::vector<SNPLine> merged;
	merged.push_back(rows[0]);
	double currentSupport = 1;
	std::cout << "merge\n";
	for (size_t i = 1; i < rows.size(); i++)
	{
		bool exists = false;
		for (size_t a = 0; a < merged.size(); a++)
		{
			if (rows[i] == merged[a])
			{
				exists = true;
				merged[a].merge(rows[i]);
			}
		}
		if (!exists)
		{
			merged.push_back(rows[i]);
		}
	}
	std::cout << "merged from " << rows.size() << " rows to " << merged.size() << " rows\n";
	std::cout << "get snpsupports\n";
	std::vector<SNPSupport> ret;
	for (auto x : merged)
	{
		std::vector<SNPSupport> newSupports = x.toSupports();
		ret.insert(ret.end(), newSupports.begin(), newSupports.end());
	}
	std::stable_sort(ret.begin(), ret.end(), [](SNPSupport left, SNPSupport right) { return left.SNPnum < right.SNPnum; });
	std::cout << "return\n";
	return ret;
}

std::pair<std::unordered_map<size_t, size_t>, std::vector<SNPSupport>> renumber(std::vector<SNPSupport> supports)
{
	std::cout << "renumber\n";
	std::unordered_map<size_t, size_t> newNumber;
	size_t nextNum = 1;
	for (auto& x : supports)
	{
		if (newNumber.count(x.readNum) == 0)
		{
			newNumber[x.readNum] = nextNum;
			nextNum++;
		}
		x.readNum = newNumber[x.readNum];
	}
	return std::pair<std::unordered_map<size_t, size_t>, std::vector<SNPSupport>>(newNumber, supports);
}

void writeRenumbering(std::unordered_map<size_t, size_t> numbering, std::string fileName)
{
	std::ofstream file { fileName };
	for (auto x : numbering)
	{
		file << x.first << " " << x.second << "\n";
	}
}

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);
	std::vector<SNPSupport> merged = mergeSupports(supports);
	auto x = renumber(merged);
	writeSupports(x.second, argv[2]);
	writeRenumbering(x.first, argv[3]);
}