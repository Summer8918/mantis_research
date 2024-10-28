#ifndef _VALIDATE_MANTIS_UTILS_H_
#define _VALIDATE_MANTIS_UTILS_H_

#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>

double calculateMean(const std::vector<uint64_t> &data);

double calculateMedian(std::vector<uint64_t> data);

double calculateStandardDeviation(const std::vector<uint64_t> &data, double mean);

std::vector<int> getRanks(const std::vector<uint64_t> &data);

#endif