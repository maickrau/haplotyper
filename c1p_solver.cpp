//g++ c1p_solver.cpp variant_utils.cpp fasta_utils.cpp -std=c++11 -o c1p_solver.exe -I../eigen
//./c1p_solver.exe inputSupportsFile renumberingFile outputSupportsFile

//http://eigen.tuxfamily.org/index.php?title=Main_Page
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <iostream>

#include "variant_utils.h"

double correlationSimilarity(const std::vector<bool>& left, const std::vector<bool>& right)
{
	double leftAverage = 0;
	double rightAverage = 0;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left[i])
		{
			leftAverage += 1;
		}
		if (right[i])
		{
			rightAverage += 1;
		}
	}
	leftAverage /= left.size();
	rightAverage /= right.size();
	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;
	for (size_t i = 0; i < left.size(); i++)
	{
		sum1 += (left[i]-leftAverage)*(right[i]-rightAverage);
		sum2 += (left[i]-leftAverage)*(left[i]-leftAverage);
		sum3 += (right[i]-rightAverage)*(right[i]-rightAverage);
	}
	return (1.0+(sum1/(sqrt(sum2)*sqrt(sum3))))/2.0;
}

double dotProductSimilarity(const std::vector<bool>& left, const std::vector<bool>& right)
{
	size_t result = 0;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left[i] && right[i])
		{
			result++;
		}
	}
	return result;
}

double jaccardSimilarity(const std::vector<bool>& left, const std::vector<bool>& right)
{
	size_t intersectSize = 0;
	size_t unionSize = 0;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left[i] == right[i])
		{
			intersectSize++;
		}
		if (left[i] || right[i])
		{
			unionSize++;
		}
	}
	assert(unionSize > 0);
	return (double)intersectSize/(double)unionSize;
}

double simpleSimilarity(const std::vector<bool>& left, const std::vector<bool>& right)
{
	double result = 0;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left[i] == right[i])
		{
			result += 1;
		}
	}
	return result;
}

template <typename F>
Eigen::MatrixXd laplacianMatrix(const std::vector<SNPSupport>& supports, F similarity)
{
	size_t maxSNP = 0;
	size_t maxRead = 0;
	for (auto x : supports)
	{
		maxSNP = std::max(maxSNP, x.SNPnum);
		maxRead = std::max(maxRead, x.readNum);
	}
	maxSNP++;
	maxRead++;
	std::vector<std::vector<bool>> binaryMatrix;
	binaryMatrix.resize(maxSNP);
	for (size_t i = 0; i < maxSNP; i++)
	{
		binaryMatrix[i].resize(maxRead, false);
	}
	for (auto x : supports)
	{
		binaryMatrix[x.SNPnum][x.readNum] = true;
	}
	Eigen::MatrixXd ret(maxSNP, maxSNP);
	for (size_t i = 0; i < maxSNP; i++)
	{
		for (size_t j = 0; j < maxSNP; j++)
		{
			ret(i, j) = similarity(binaryMatrix[i], binaryMatrix[j]);
		}
	}
	std::vector<double> diagonal;
	diagonal.resize(maxSNP, 0);
	for (size_t i = 0; i < maxSNP; i++)
	{
		for (size_t j = 0; j < maxSNP; j++)
		{
			diagonal[i] += ret(i, j);
		}
	}
	for (size_t i = 0; i < maxSNP; i++)
	{
		ret(i, i) -= diagonal[i];
	}
	for (size_t i = 0; i < maxSNP; i++)
	{
		for (size_t j = 0; j < maxSNP; j++)
		{
			ret(i, j) *= -1;
		}
	}

	return ret;
}

template <typename F>
std::vector<double> getFiedlerVector(const std::vector<SNPSupport>& supports, F similarity)
{
	Eigen::MatrixXd matrix = laplacianMatrix(supports, similarity);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
	auto fiedlerVector = solver.eigenvectors().col(1);
	std::vector<double> values;
	for (size_t i = 0; i < matrix.cols(); i++)
	{
		values.push_back(fiedlerVector(i));
	}
	return values;
}

template <typename F>
SupportRenumbering getSpectralOrdering(const std::vector<SNPSupport>& supports, F similarity)
{
	std::vector<double> fiedlerVector = getFiedlerVector(supports, similarity);
	std::vector<std::pair<size_t, double>> values;
	values.resize(fiedlerVector.size());
	for (size_t i = 0; i < fiedlerVector.size(); i++)
	{
		values[i] = std::pair<size_t, double> { i, fiedlerVector[i] };
	}
	std::sort(values.begin(), values.end(), [](std::pair<size_t, double> left, std::pair<size_t, double> right) { return left.second < right.second; });
	SupportRenumbering result;
	size_t maxSNP = 0;
	size_t maxRead = 0;
	for (auto x : supports)
	{
		maxSNP = std::max(maxSNP, x.SNPnum);
		maxRead = std::max(maxRead, x.readNum);
	}
	maxSNP++;
	maxRead++;
	for (size_t i = 0; i < values.size(); i++)
	{
		result.addSNPRenumbering(values[i].first, i);
	}
	for (size_t i = 0; i < maxRead; i++)
	{
		result.addReadRenumbering(i, i);
	}
	return result;
}

double getScore(const std::vector<SNPSupport>& supports)
{
	size_t maxSNP = 0;
	size_t maxRead = 0;
	for (auto x : supports)
	{
		maxSNP = std::max(maxSNP, x.SNPnum);
		maxRead = std::max(maxRead, x.readNum);
	}
	maxRead++;
	maxSNP++;
	std::vector<std::pair<size_t, size_t>> rowExtents;
	rowExtents.resize(maxRead, {-1, 0});
	for (auto x : supports)
	{
		rowExtents[x.readNum].first = std::min(rowExtents[x.readNum].first, x.SNPnum);
		rowExtents[x.readNum].second = std::max(rowExtents[x.readNum].second, x.SNPnum);
	}
	double result = 0;
	for (size_t i = 0; i < maxSNP; i++)
	{
		double coverage = 0;
		for (size_t j = 0; j < maxRead; j++)
		{
			if (rowExtents[j].first <= i && rowExtents[j].second >= i)
			{
				coverage += 1;
			}
		}
		result += pow(2, coverage);
	}
	return result;
}

SupportRenumbering sortRows(const SupportRenumbering& renumbering, const std::vector<SNPSupport>& supports)
{
	SupportRenumbering ret;
	for (size_t i = 0; i < renumbering.SNPSize(); i++)
	{
		ret.addSNPRenumbering(i, renumbering.getSNPRenumbering(i));
	}
	std::vector<std::pair<size_t, size_t>> minSNP;
	for (size_t i = 0; i < renumbering.readSize(); i++)
	{
		minSNP.emplace_back(i, -1);
	}
	for (auto x : supports)
	{
		minSNP[x.readNum].second = std::min(minSNP[x.readNum].second, ret.getSNPRenumbering(x.SNPnum));
	}
	std::sort(minSNP.begin(), minSNP.end(), [](std::pair<size_t, size_t> left, std::pair<size_t, size_t> right) { return left.second < right.second; });
	for (size_t i = 0; i < minSNP.size(); i++)
	{
		ret.addReadRenumbering(minSNP[i].first, i);
	}
	return ret;
}

template <typename F>
std::tuple<double, SupportRenumbering, std::vector<SNPSupport>> getRenumbered(const std::vector<SNPSupport>& supports, F similarity)
{
	SupportRenumbering renumbering = getSpectralOrdering(supports, similarity);
	renumbering = sortRows(renumbering, supports);
	std::vector<SNPSupport> result = renumberSupports(supports, renumbering);
	double score = getScore(result);
	return std::tuple<double, SupportRenumbering, std::vector<SNPSupport>> { score, renumbering, result };
}

SupportRenumbering identityRenumbering(std::vector<SNPSupport> supports)
{
	size_t maxSNP = 0;
	size_t maxRead = 0;
	for (auto x : supports)
	{
		maxSNP = std::max(maxSNP, x.SNPnum);
		maxRead = std::max(maxRead, x.readNum);
	}
	maxSNP++;
	maxRead++;
	SupportRenumbering ret;
	for (size_t i = 0; i < maxSNP; i++)
	{
		ret.addSNPRenumbering(i, i);
	}
	for (size_t i = 0; i < maxRead; i++)
	{
		ret.addReadRenumbering(i, i);
	}
	return ret;
}

int main(int argc, char** argv)
{
	std::vector<SNPSupport> supports = loadSupports(argv[1]);
	std::cerr << "original score: " << getScore(supports) << "\n";
	auto simple = getRenumbered(supports, simpleSimilarity);
	std::cerr << "simple score: " << std::get<0>(simple) << "\n";
	auto jaccard = getRenumbered(supports, jaccardSimilarity);
	std::cerr << "jaccard score: " << std::get<0>(jaccard) << "\n";
	auto dotProduct = getRenumbered(supports, dotProductSimilarity);
	std::cerr << "dot product score: " << std::get<0>(dotProduct) << "\n";
	auto correlation = getRenumbered(supports, correlationSimilarity);
	std::cerr << "correlation score: " << std::get<0>(correlation) << "\n";

	std::tuple<double, SupportRenumbering, std::vector<SNPSupport>> best;
	std::get<0>(best) = getScore(supports);
	std::get<1>(best) = identityRenumbering(supports);
	std::get<2>(best) = supports;

	if (std::get<0>(simple) < std::get<0>(best))
	{
		best = simple;
	}
	if (std::get<0>(jaccard) < std::get<0>(best))
	{
		best = jaccard;
	}
	if (std::get<0>(dotProduct) < std::get<0>(best))
	{
		best = dotProduct;
	}
	if (std::get<0>(correlation) < std::get<0>(best))
	{
		best = correlation;
	}

	writeRenumbering(std::get<1>(best), argv[2]);
	writeSupports(std::get<2>(best), argv[3]);
}