#include "sparsepp/spp.h"
#include "tsl/sparse_map.h"

#include "MantisFS.h"
#include "ProgOpts.h"
#include "coloreddbg.h"
#include "squeakrconfig.h"
#include "json.hpp"
#include "mantis_utils.hpp"
#include "mantisconfig.hpp"

#include<vector>
#include<iostream>

void hash128DebugHelper(__uint128_t vec_hash) {
    // Extracting the lower and higher 64 bits
    uint64_t lower_part = (uint64_t)vec_hash;  // Lower 64 bits
    uint64_t higher_part = (uint64_t)(vec_hash >> 64);  // Higher 64 bits

    // Printing the two parts
    std::cout << "vec_hash (lower 64 bits): " << lower_part << std::endl;
    std::cout << "vec_hash (higher 64 bits): " << higher_part << std::endl;
}



int main() {
    std::vector<int> vec1 = {0,0,0,1,0,1};
    std::vector<int> vec2 = {0,0,0,1,0,1};
    __uint128_t vec1_hash = MurmurHash128A((void*)vec1.data(),
											(((vec1.size() * 32 + 63) >> 6) << 6) / 8, 2038074743,
											2038074751);
    __uint128_t vec2_hash = MurmurHash128A((void*)vec2.data(),
											(((vec2.size() * 32 + 63) >> 6) << 6) / 8, 2038074743,
											2038074751);
    hash128DebugHelper(vec1_hash);
    hash128DebugHelper(vec2_hash);
	vec1_hash = MurmurHash128A((void*)vec1.data(),
											(((vec1.size() * 32 + 63) >> 6) << 6) / 8, 2038074743,
											2038074751);
	vec2_hash = MurmurHash128A((void*)vec2.data(),
											(((vec2.size() * 32 + 63) >> 6) << 6) / 8, 2038074743,
											2038074751);
	hash128DebugHelper(vec1_hash);
    hash128DebugHelper(vec2_hash);	
    vec1 = {1,0,0,1,0};
    vec1_hash = MurmurHash128A((void*)vec1.data(),
											(((vec1.size() * 32 + 63) >> 6) << 6) / 8, 2038074743,
											2038074751);	
    hash128DebugHelper(vec1_hash);
    return 0;
}