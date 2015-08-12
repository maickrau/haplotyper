#include <algorithm>
#include <cassert>
#include <unordered_map>

#include "variant_utils.h"

SNPPosition::SNPPosition(Genome gen1, Genome gen2) : leftFlank(), rightFlank(), variants()
{
	size_t loc = 0;
	while (gen1.bases[loc] == gen2.bases[loc])
	{
		loc++;
	}
	leftFlank = gen1.bases.substr(0, loc);
	rightFlank = gen1.bases.substr(loc+1);
	variants.insert(gen1.bases[loc]);
	variants.insert(gen2.bases[loc]);
	std::reverse(leftFlank.begin(), leftFlank.end());
}

SNPSupport::SNPSupport(size_t readNum, size_t SNPnum, char variant) : readNum(readNum), SNPnum(SNPnum), support(1), variant(variant) 
{
}

SNPSupport::SNPSupport(size_t readNum, size_t SNPnum, char variant, double support) : readNum(readNum), SNPnum(SNPnum), support(support), variant(variant) 
{
}

std::vector<SNPSupport> loadSupports(std::string fileName)
{
	std::vector<SNPSupport> result;
	std::ifstream file { fileName };
	do
	{
		size_t readNum, SNPnum;
		char variant;
		double support;
		file >> readNum >> SNPnum >> variant >> support;
		if (file.good())
		{
			result.emplace_back(readNum, SNPnum, variant, support);
		}
	} while (file.good());
	return result;
}

void writeSupports(std::vector<SNPSupport> supports, std::string fileName)
{
	std::ofstream file { fileName };
	for (auto x : supports)
	{
		file << x.readNum << " " << x.SNPnum << " " << x.variant << " " << x.support << "\n";
	}
}

SupportRenumbering SupportRenumbering::identity(size_t maxRead, size_t maxSNP)
{
	SupportRenumbering ret;
	for (size_t i = 0; i < maxRead; i++)
	{
		ret.addReadRenumbering(i, i);
	}
	for (size_t i = 0; i < maxSNP; i++)
	{
		ret.addSNPRenumbering(i, i);
	}
	return ret;
}

SupportRenumbering::SupportRenumbering()
{
}

bool SupportRenumbering::checkValidity() const
{
	if (readRenumbering.size() == 0)
	{
		return false;
	}
	if (SNPRenumbering.size() == 0)
	{
		return false;
	}
	std::vector<bool> usedReads;
	usedReads.resize(readRenumbering.size(), false);
	for (size_t i = 0; i < readRenumbering.size(); i++)
	{
		if (readRenumbering[i] >= readRenumbering.size())
		{
			return false;
		}
		if (usedReads[readRenumbering[i]])
		{
			return false;
		}
		usedReads[readRenumbering[i]] = true;
	}

	std::vector<bool> usedSNPs;
	usedSNPs.resize(SNPRenumbering.size(), false);
	for (size_t i = 0; i < SNPRenumbering.size(); i++)
	{
		if (SNPRenumbering[i] >= SNPRenumbering.size())
		{
			return false;
		}
		if (usedSNPs[SNPRenumbering[i]])
		{
			return false;
		}
		usedSNPs[SNPRenumbering[i]] = true;
	}

	for (size_t i = 0; i < readRenumbering.size(); i++)
	{
		if (!usedReads[i])
		{
			return false;
		}
	}
	for (size_t i = 0; i < SNPRenumbering.size(); i++)
	{
		if (!usedSNPs[i])
		{
			return false;
		}
	}
	return true;
}

std::vector<SNPSupport> renumberSupports(std::vector<SNPSupport> supports, SupportRenumbering renumbering)
{
	std::vector<SNPSupport> ret;
	for (auto x : supports)
	{
		ret.push_back(x);
		ret.back().readNum = renumbering.getReadRenumbering(ret.back().readNum);
		ret.back().SNPnum = renumbering.getSNPRenumbering(ret.back().SNPnum);
	}
	return ret;
}

SupportRenumbering SupportRenumbering::merge(SupportRenumbering second)
{
	SupportRenumbering ret;
	for (size_t i = 0; i < readRenumbering.size(); i++)
	{
		if (readRenumbering[i] == -1)
		{
			ret.addReadRenumbering(i, -1);
		}
		else
		{
			assert(second.readRenumbering.size() > readRenumbering[i]);
			if (second.hasReadRenumbering(readRenumbering[i]))
			{
				ret.addReadRenumbering(i, second.getReadRenumbering(readRenumbering[i]));
			}
			else
			{
				ret.addReadRenumbering(i, -1);
			}
		}
	}
	for (size_t i = 0; i < SNPRenumbering.size(); i++)
	{
		if (SNPRenumbering[i] == -1)
		{
			ret.addSNPRenumbering(i, -1);
		}
		else
		{
			assert(second.SNPRenumbering.size() > SNPRenumbering[i]);
			if (second.hasSNPRenumbering(SNPRenumbering[i]))
			{
				ret.addSNPRenumbering(i, second.getSNPRenumbering(SNPRenumbering[i]));
			}
			else
			{
				ret.addSNPRenumbering(i, -1);
			}
		}
	}
	return ret;
}

SupportRenumbering SupportRenumbering::reverse()
{
	SupportRenumbering ret;
	for (size_t i = 0; i < readRenumbering.size(); i++)
	{
		ret.addReadRenumbering(readRenumbering[i], i);
	}
	for (size_t i = 0; i < SNPRenumbering.size(); i++)
	{
		ret.addSNPRenumbering(SNPRenumbering[i], i);
	}
	return ret;
}

void SupportRenumbering::overwriteReadRenumbering(size_t oldRead, size_t newRead)
{
	assert(oldRead < readRenumbering.size());
	readRenumbering[oldRead] = newRead;
}

void SupportRenumbering::addReadRenumbering(size_t oldRead, size_t newRead)
{
	if (readRenumbering.size() <= oldRead)
	{
		readRenumbering.resize(oldRead+1, -1);
	}
	if (readRenumbering[oldRead] != -1)
	{
		assert(readRenumbering[oldRead] == newRead);
	}
	readRenumbering[oldRead] = newRead;
}

void SupportRenumbering::addSNPRenumbering(size_t oldSNP, size_t newSNP)
{
	if (SNPRenumbering.size() <= oldSNP)
	{
		SNPRenumbering.resize(oldSNP+1, -1);
	}
	if (SNPRenumbering[oldSNP] != -1)
	{
		assert(SNPRenumbering[oldSNP] == newSNP);
	}
	SNPRenumbering[oldSNP] = newSNP;
}

bool SupportRenumbering::hasReadRenumbering(size_t oldRead) const
{
	return readRenumbering[oldRead] != -1;
}

bool SupportRenumbering::hasSNPRenumbering(size_t oldSNP) const
{
	return SNPRenumbering[oldSNP] != -1;
}

size_t SupportRenumbering::getReadRenumbering(size_t oldRead) const
{
	assert(readRenumbering.size() > oldRead);
	assert(readRenumbering[oldRead] != -1);
	return readRenumbering[oldRead];
}

size_t SupportRenumbering::getSNPRenumbering(size_t oldSNP) const
{
	assert(SNPRenumbering.size() > oldSNP);
	assert(SNPRenumbering[oldSNP] != -1);
	return SNPRenumbering[oldSNP];
}

size_t SupportRenumbering::readSize() const
{
	return readRenumbering.size();
}

size_t SupportRenumbering::SNPSize() const
{
	return SNPRenumbering.size();
}

void SupportRenumbering::swapRows(size_t firstRow, size_t secondRow)
{
	assert(firstRow < readRenumbering.size());
	assert(secondRow < readRenumbering.size());
	assert(readRenumbering[firstRow] != -1);
	assert(readRenumbering[secondRow] != -1);
	std::swap(readRenumbering[firstRow], readRenumbering[secondRow]);
}

void SupportRenumbering::swapColumns(size_t firstColumn, size_t secondColumn)
{
	assert(firstColumn < SNPRenumbering.size());
	assert(secondColumn < SNPRenumbering.size());
	assert(SNPRenumbering[firstColumn] != -1);
	assert(SNPRenumbering[secondColumn] != -1);
	std::swap(SNPRenumbering[firstColumn], SNPRenumbering[secondColumn]);
}

void writeRenumbering(SupportRenumbering renumbering, std::string fileName)
{
	std::ofstream file { fileName };
	file << renumbering.readRenumbering.size();
	for (auto x : renumbering.readRenumbering)
	{
		file << " " << x;
	}
	file << "\n" << renumbering.SNPRenumbering.size();
	for (auto x : renumbering.SNPRenumbering)
	{
		file << " " << x;
	}
	file << "\n";
}

SupportRenumbering loadRenumbering(std::string fileName)
{
	SupportRenumbering ret;
	std::ifstream file { fileName };
	size_t reads;
	file >> reads;
	for (size_t i = 0; i < reads; i++)
	{
		size_t readPos;
		file >> readPos;
		ret.addReadRenumbering(i, readPos);
	}
	size_t SNPs;
	file >> SNPs;
	for (size_t i = 0; i < SNPs; i++)
	{
		size_t SNPpos;
		file >> SNPpos;
		ret.addSNPRenumbering(i, SNPpos);
	}
	assert(file.good());
	return ret;
}

SNPLine::SNPLine() {};

bool SNPLine::operator==(const SNPLine& second) const
{
	return variantsAtLocations == second.variantsAtLocations;
}

bool SNPLine::operator!=(const SNPLine& second) const
{
	return !(*this == second);
}

bool SNPLine::contains(const SNPLine& second) const
{
	if (variantsAtLocations.size() < second.variantsAtLocations.size())
	{
		return false;
	}
	size_t index = 0;
	size_t secondIndex = 0;
	while (index < variantsAtLocations.size() && secondIndex < second.variantsAtLocations.size())
	{
		if (variantsAtLocations[index].first < second.variantsAtLocations[secondIndex].first)
		{
			index++;
		}
		else if (variantsAtLocations[index].first == second.variantsAtLocations[secondIndex].first)
		{
			if (variantsAtLocations[index].second != second.variantsAtLocations[secondIndex].second)
			{
				return false;
			}
			index++;
			secondIndex++;
		}
		else if (variantsAtLocations[index].first > second.variantsAtLocations[secondIndex].first)
		{
			return false;
		}
	}
	if (index == variantsAtLocations.size() && secondIndex < second.variantsAtLocations.size())
	{
		return false;
	}
	return true;
}

void SNPLine::mergeSubset(SNPLine second)
{
	size_t index = 0;
	size_t secondIndex = 0;
	while (index < variantsAtLocations.size() && secondIndex < second.variantsAtLocations.size())
	{
		if (variantsAtLocations[index].first < second.variantsAtLocations[secondIndex].first)
		{
			index++;
		}
		else if (variantsAtLocations[index].first == second.variantsAtLocations[secondIndex].first)
		{
			assert(variantsAtLocations[index].second == second.variantsAtLocations[secondIndex].second);
			supportsAtLocations[index] += second.supportsAtLocations[secondIndex];
			index++;
			secondIndex++;
		}
		else
		{
			assert(false);
		}
	}
}

std::vector<SNPSupport> SNPLine::toSupports() const
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

char SNPLine::variantAt(size_t loc) const
{
	for (auto x : variantsAtLocations)
	{
		if (x.first == loc)
		{
			return x.second;
		}
	}
	assert(false);
}

std::vector<SNPLine> makeLines(std::vector<SNPSupport> supports)
{
	std::sort(supports.begin(), supports.end(), [](SNPSupport left, SNPSupport right) { return left.SNPnum < right.SNPnum; });
	std::set<size_t> reads;
	for (auto x : supports)
	{
		reads.insert(x.readNum);
	}
	std::vector<SNPLine> result;
	for (auto x : reads)
	{
		result.emplace_back(supports.begin(), supports.end(), x);
	}
	return result;
}

std::vector<SNPSupport> mergeRows(std::vector<SNPSupport> oldSupports, size_t row1, size_t row2)
{
	std::vector<SNPSupport> result;
	std::unordered_map<size_t, SNPSupport> newRow;
	for (auto x : oldSupports)
	{
		if (x.readNum != row1 && x.readNum != row2)
		{
			result.push_back(x);
		}
		else
		{
			if (newRow.find(x.SNPnum) == newRow.end())
			{
				newRow.emplace(x.SNPnum, x);
				newRow.at(x.SNPnum).readNum = row1;
			}
			else
			{
				assert(newRow.at(x.SNPnum).variant == x.variant);
				newRow.at(x.SNPnum).support += x.support;
			}
		}
	}
	for (auto x : newRow)
	{
		result.push_back(x.second);
	}
	return result;
}
