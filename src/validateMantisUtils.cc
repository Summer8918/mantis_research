#include "validateMantisUtils.h"


double calculateMean(const std::vector<uint64_t> &data) {
	int sum = accumulate(data.begin(), data.end(), 0);
	return static_cast<double>(sum) / data.size();
}

double calculateMedian(std::vector<uint64_t> data) {
	sort(data.begin(), data.end());
	int n = data.size();
	if (n % 2 == 0) {
		return (data[n / 2 - 1] + data[n / 2]) / 2.0;
	} else {
		return data[n / 2];
	}
}

double calculateStandardDeviation(const std::vector<uint64_t> &data, double mean) {
	double sum = 0.0;
	for (auto val : data) {
		sum += pow(val - mean, 2);
	}
	return sqrt(sum / data.size());
}

std::vector<int> getRanks(const std::vector<uint64_t> &data) {
    int len = data.size();
    std::vector<int> indices(len), ranks(len);
    for (int i = 0; i < len; ++i) {
        indices[i] = i;
    }

    // Sort indices based on values in A in descending order
    std::sort(indices.begin(), indices.end(), [&data](int i, int j) {
        return data[i] > data[j];
    });

    for (int r = 0; r < len; ++r) {
        ranks[indices[r]] = r + 1;
    }
    return ranks;
}