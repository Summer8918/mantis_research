#include "validateMantisUtils.h"
#include <iostream>

void test_mean_median_stddev() {
	std::vector<uint64_t> data = {2, 1, 3, 4, 5};
	double mean = calculateMean(data);
	double median = calculateMedian(data);
	double stdDev = calculateStandardDeviation(data, mean);
    std::cout << "mean:" << mean << " median:" << median << " stdDev:" << stdDev << std::endl;
}

void test_get_rank() {
	std::vector<uint64_t> data = {2, 1, 3, 4, 5};
	std::vector<int> ranks = getRanks(data);
	std::cout << "Ranks:" << std::endl;
	for (int i = 0; i < data.size(); ++i) {
		std::cout << "data[i]:" << data[i] << " rank:" << ranks[i] << std::endl;
	}
}

int main() {
    test_mean_median_stddev();
	test_get_rank();
    return 0;
}