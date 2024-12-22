#include <vector>
#include <unordered_map>
#include <iostream>
#include <cassert>
#include <chrono>

#define MAX_LL_INTEGER (0xFFFFFFFFFFFFFFFF)
#define SECOND_MAX_LL_INTEGER (0xFFFFFFFFFFFFFFFE)
#define ROW_DIST_THRESHOLD 0x10

struct pair_hash {
    template <class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2>& p) const {
        auto hash1 = std::hash<T1>{}(p.first);
    	auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ (hash2 << 1); // Combine the two hash values
    }
};



void test_compress() {
    uint64_t num_samples = 50;
	uint64_t len = 236423;
    std::vector<uint64_t> iv_buffer(len, 0);
	auto start = std::chrono::high_resolution_clock::now();
    for (uint64_t i = 0; i < iv_buffer.size(); i++) {
        iv_buffer[i] = i + 10;
    }
	auto end = std::chrono::high_resolution_clock::now();
	// Calculate the duration
    std::chrono::duration<double> duration = end - start;
    std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;

    uint64_t rows = iv_buffer.size() / num_samples;
	// Define a custom hash function
	std::cout << "Start to compress iv buffer" << std::endl;
	std::unordered_map<std::pair<uint64_t, uint64_t>, int, pair_hash> row_vec_dists;

    for (uint64_t i = 0; i < rows; i++) {
		uint64_t i_idx = i * num_samples;
		if (i % 100 == 0) {
        	std::cout << "i:" << i << std::endl;
		}
		for (uint64_t j = i + 1; j < rows; j++) {
			uint64_t j_idx = j * num_samples;
			uint64_t ManhattanDist = 0;
			for (uint64_t k = 0; k < num_samples; k++) {
				if (iv_buffer[j_idx + k] > iv_buffer[i_idx + k]) {
					ManhattanDist += iv_buffer[j_idx + k] - iv_buffer[i_idx + k];
				} else {
					ManhattanDist += iv_buffer[i_idx + k] -  iv_buffer[j_idx + k];
				}
			}
			std::pair<uint64_t, uint64_t> p = {i, j};
			row_vec_dists[p] = ManhattanDist;
		}
	}
    std::cout << "Finishes cal manhattan distance for each two row vectors" << std::endl;

    // Form sets for compression. The distance of the any two rows vectors in the set is <= threshold
	//<key, value> : <eq class id : combined set id>
	std::vector<uint64_t> combined_set_id_map(rows, MAX_LL_INTEGER);
	// <key, value> : <combined set id : set> 
	std::unordered_map<int, std::vector<uint64_t>> combined_set;
	uint32_t combined_set_id = 0;

    for (auto it = row_vec_dists.begin(); it != row_vec_dists.end(); it++) {
		if (it->second <= ROW_DIST_THRESHOLD) {
			uint64_t eq_class_id1 = it->first.first;
			uint64_t eq_class_id2 = it->first.second;
			bool in_a_set_flag = true;
			if (combined_set_id_map[eq_class_id1] == MAX_LL_INTEGER && combined_set_id_map[eq_class_id2 == MAX_LL_INTEGER]) {
				combined_set_id_map[eq_class_id1] = combined_set_id;
                combined_set_id_map[eq_class_id2] = combined_set_id;
				combined_set[combined_set_id].push_back({eq_class_id1});
                combined_set[combined_set_id].push_back({eq_class_id2});
				combined_set_id++;
			} else if (combined_set_id_map[eq_class_id1] != MAX_LL_INTEGER && combined_set_id_map[eq_class_id2] == MAX_LL_INTEGER) {
				in_a_set_flag = true;
				for (auto &node : combined_set[combined_set_id_map[eq_class_id1]]) {
					std::pair<uint64_t, uint64_t> p = {node, eq_class_id2};
					if (node > eq_class_id2) {
						p = {eq_class_id2, node};
					}
					assert(row_vec_dists.find(p) != row_vec_dists.end());
					if (row_vec_dists[p] > ROW_DIST_THRESHOLD) {
						in_a_set_flag = false;
						break;
					}
				}
				if (in_a_set_flag) {
					combined_set[combined_set_id_map[eq_class_id1]].push_back(eq_class_id2);
					combined_set_id_map[eq_class_id2] = combined_set_id_map[eq_class_id1];
				}
			} else if (combined_set_id_map[eq_class_id2] != MAX_LL_INTEGER && combined_set_id_map[eq_class_id1] == MAX_LL_INTEGER) {
				in_a_set_flag = true;
				for (auto &node : combined_set[combined_set_id_map[eq_class_id2]]) {
					std::pair<uint64_t, uint64_t> p = {node, eq_class_id1};
					if (node > eq_class_id1) {
						p = {eq_class_id1, node};
					}
					assert(row_vec_dists.find(p) != row_vec_dists.end());
					if (row_vec_dists[p] > ROW_DIST_THRESHOLD) {
						in_a_set_flag = false;
						break;
					}
				}
				if (in_a_set_flag) {
					combined_set[combined_set_id_map[eq_class_id2]].push_back(eq_class_id1);
					combined_set_id_map[eq_class_id1] = combined_set_id_map[eq_class_id2];
				}
			} else {
				// When eq_class_id1 and eq_class_id2 belongs to different set, try to merge them.
				in_a_set_flag = true;
				for (auto &node1 : combined_set[combined_set_id_map[eq_class_id1]]) {
					for (auto &node2 : combined_set[combined_set_id_map[eq_class_id2]]) {
						std::pair<uint64_t, uint64_t> p = {node1, node2};
						if (node1 > node2) {
							p = {node2, node1};
						}
						if (row_vec_dists[p] > ROW_DIST_THRESHOLD) {
							in_a_set_flag = false;
							break;
						}
					}
					if (in_a_set_flag == false) {
						break;
					}
				}
				if (in_a_set_flag) {
					for (auto &node2 : combined_set[combined_set_id_map[eq_class_id2]]) {
						combined_set[combined_set_id_map[eq_class_id1]].push_back(node2);
					}
					combined_set[combined_set_id_map[eq_class_id2]] = {};
					combined_set_id_map[eq_class_id2] = combined_set_id_map[eq_class_id1];
				}
			}
		}
	}

	std::cout << "Finishes combine row vectors" << std::endl;
	// Generate compressed iv_buffer
	uint64_t new_eqclass_id = 0;
	std::vector<uint64_t> new_iv_buffer(iv_buffer.size(), 0);
	std::unordered_map<int, int> prev_eqclassid_to_new_eqclassid;
	for (int i = 0; i < combined_set_id_map.size(); i++) {
		// -2 is the visited row
		if (combined_set_id_map[i] != MAX_LL_INTEGER && combined_set_id_map[i] != SECOND_MAX_LL_INTEGER) {
			std::vector<uint64_t> row_set = combined_set[combined_set_id_map[i]];
			std::vector<uint64_t> new_row_vec(num_samples, 0);
			for (int j = 0; j < num_samples; j++) {
				for (auto &row : row_set) {
					if (iv_buffer[row * num_samples + j] > new_row_vec[i]) {
						new_row_vec[j] = iv_buffer[row * num_samples + i];
					}
				}
			}
			for (int j = 0; j < num_samples; j++) {
				new_iv_buffer[new_eqclass_id + j] = new_row_vec[j];
			}

			for (auto &row : row_set) {
				combined_set_id_map[row] = SECOND_MAX_LL_INTEGER;
				prev_eqclassid_to_new_eqclassid[row] = new_eqclass_id;
			}
			new_eqclass_id++;
		} else if (combined_set_id_map[i] == MAX_LL_INTEGER) {
			for (int j = 0; j < num_samples; j++) {
				new_iv_buffer[new_eqclass_id + j] = iv_buffer[i * num_samples + j];
			}
			combined_set_id_map[i] == SECOND_MAX_LL_INTEGER;
			new_eqclass_id++;
		}
	}
	std::cout << "Finishes generate compressed iv_buffer" << std::endl;
}

// Custom hash function for std::vector
template <typename T>
struct VectorHash {
    std::size_t operator()(const std::vector<T>& vec) const {
        std::size_t seed = 0;
        std::hash<T> hasher;
        for (const auto& elem : vec) {
            // Combine hashes of individual elements
            seed ^= hasher(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

void test_compress2() {
	uint64_t num_samples = 3;
	uint64_t rows = 3;
    std::vector<uint64_t> iv_buffer = {1, 2, 3, 100, 200, 300, 2, 3, 4};

	std::unordered_map<uint64_t, std::vector<uint64_t>> normalized_row_vecs_hash;
	for (uint64_t i = 0; i < rows; i++) {
		std::vector<uint64_t> normalized_counter(num_samples + 1, 0);
		uint64_t tmp_sum = 0;
		for (uint64_t j = 0; j < num_samples; j++) {
			normalized_counter[j] = iv_buffer[i * num_samples + j] / (ROW_DIST_THRESHOLD + 1);
			tmp_sum += iv_buffer[i * num_samples + j];
		}
		normalized_counter[num_samples] = tmp_sum / (ROW_DIST_THRESHOLD + 1);
		VectorHash<uint64_t> vectorHasher;
		uint64_t hash_val = vectorHasher(normalized_counter);
		normalized_row_vecs_hash[hash_val].push_back(i);
	}
	for (auto it = normalized_row_vecs_hash.begin(); it != normalized_row_vecs_hash.end(); it++) {
		for (auto &data : it->second) {
			std::cout << "data:" << data << " ";
		}
		std::cout << " end" <<  std::endl;
	}

	// Form sets for compression. The distance of the any two rows vectors in the set is <= threshold
	// <key, value> : <eq class id : combined set id>
	// MAX_LL_INTEGER means not assign set id, every set has a id
	std::vector<uint64_t> combined_set_id_map(rows, MAX_LL_INTEGER);
	// <key, value> : <set id : row set>
	std::unordered_map<uint64_t, std::vector<uint64_t>> combined_set;
	uint64_t combined_set_id = 0;
	uint64_t tmp_cnt = 0;
	for (auto it = normalized_row_vecs_hash.begin(); it != normalized_row_vecs_hash.end(); it++) {
		for (auto &data : it->second) {
			combined_set_id_map[data] = combined_set_id;
			combined_set[combined_set_id].push_back(data);
		}
		tmp_cnt += it->second.size();
		combined_set_id++;
		std::cout << "completed percent:" << 1.0 * tmp_cnt / rows << std::endl;
	}
	for (int i = 0; i < combined_set_id_map.size(); i++) {
		std::cout << "combined_set_id_map: " << i << " val:" << combined_set_id_map[i] << std::endl;
	}
	for (auto it = combined_set.begin(); it != combined_set.end(); it++) {
		for (auto &data : it->second) {
			std::cout << "combined_set data:" << data << " ";
		}
		std::cout << " end" <<  std::endl;
	}
}

/*

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::compress_iv_buffer_brutal_force() {

	// column = num_samples
	uint64_t rows = iv_buffer.size() / num_samples;    // 20,000,000
	if (get_num_eqclasses() % mantis::NUM_IV_BUFFER != 0) { // when iv buffer is not full
		rows = get_num_eqclasses() * num_samples;
	}
	
	std::cout << "Start to compress iv buffer" << std::endl;
	std::cout << "rows:" << rows << std::endl;
	// Use a custom hash function for pair hash
	std::unordered_map<std::pair<uint64_t, uint64_t>, int, pair_hash> row_vec_dists;
	// Calculate Manhattan Distance for each two row vectors
	for (uint64_t i = 0; i < rows; i++) {
		uint64_t i_idx = i * num_samples;   // index for row i
		std::cout << "i:" << i << std::endl;
		for (uint64_t j = i + 1; j < rows; j++) {
			uint64_t j_idx = j * num_samples;   // index for row j
			uint64_t ManhattanDist = 0;
			for (uint64_t k = 0; k < num_samples; k++) {
				if (iv_buffer[j_idx + k] > iv_buffer[i_idx + k]) {
					ManhattanDist += iv_buffer[j_idx + k] - iv_buffer[i_idx + k];
				} else {
					ManhattanDist += iv_buffer[i_idx + k] -  iv_buffer[j_idx + k];
				}
			}
			std::pair<uint64_t, uint64_t> p = {i, j};  // always i < j
			row_vec_dists[p] = ManhattanDist;
		}
	}
	std::cout << "Finishes cal manhattan distance for each two row vectors" << std::endl;
	// Form sets for compression. The distance of the any two rows vectors in the set is <= threshold
	//<key, value> : <eq class id : combined set id>
	// MAX_LL_INTEGER means not assign set id, every set has a id
	std::vector<uint64_t> combined_set_id_map(rows, MAX_LL_INTEGER);
	// <key, value> : <set id : row set> 
	std::unordered_map<int, std::vector<uint64_t>> combined_set;
	uint32_t combined_set_id = 0;
	for (auto it = row_vec_dists.begin(); it != row_vec_dists.end(); it++) {
		if (it->second <= ROW_DIST_THRESHOLD) {
			uint64_t eq_class_id1 = it->first.first;
			uint64_t eq_class_id2 = it->first.second;
			bool in_a_set_flag = true;
			// row1 and row2 do not belong to any set neither
			if (combined_set_id_map[eq_class_id1] == MAX_LL_INTEGER && combined_set_id_map[eq_class_id2 == MAX_LL_INTEGER]) {
				combined_set_id_map[eq_class_id1] = combined_set_id;
                combined_set_id_map[eq_class_id2] = combined_set_id;
				combined_set[combined_set_id].push_back({eq_class_id1});
                combined_set[combined_set_id].push_back({eq_class_id2});
				combined_set_id++;
			}
			// row1 belong to one set, row2 donot belong to any set, check if row2 is close to the set row1 belong
			else if (combined_set_id_map[eq_class_id1] != MAX_LL_INTEGER && combined_set_id_map[eq_class_id2] == MAX_LL_INTEGER) {
				in_a_set_flag = true;
				for (auto &node : combined_set[combined_set_id_map[eq_class_id1]]) {
					std::pair<uint64_t, uint64_t> p = {node, eq_class_id2};
					if (node > eq_class_id2) {
						p = {eq_class_id2, node};
					}
					assert(row_vec_dists.find(p) != row_vec_dists.end());
					if (row_vec_dists[p] > ROW_DIST_THRESHOLD) {
						in_a_set_flag = false;
						break;
					}
				}
				if (in_a_set_flag) {
					combined_set[combined_set_id_map[eq_class_id1]].push_back(eq_class_id2);
					combined_set_id_map[eq_class_id2] = combined_set_id_map[eq_class_id1];
				}
			}
			// row2 belong to one set, row1 donot belong to any set, check if row1 is close to the set row2 belong
			else if (combined_set_id_map[eq_class_id2] != MAX_LL_INTEGER && combined_set_id_map[eq_class_id1] == MAX_LL_INTEGER) {
				in_a_set_flag = true;
				for (auto &node : combined_set[combined_set_id_map[eq_class_id2]]) {
					std::pair<uint64_t, uint64_t> p = {node, eq_class_id1};
					if (node > eq_class_id1) {
						p = {eq_class_id1, node};
					}
					assert(row_vec_dists.find(p) != row_vec_dists.end());
					if (row_vec_dists[p] > ROW_DIST_THRESHOLD) {
						in_a_set_flag = false;
						break;
					}
				}
				if (in_a_set_flag) {
					combined_set[combined_set_id_map[eq_class_id2]].push_back(eq_class_id1);
					combined_set_id_map[eq_class_id1] = combined_set_id_map[eq_class_id2];
				}
			}
			// When eq_class_id1 and eq_class_id2 belongs to different set, try to merge them.
			else {
				in_a_set_flag = true;
				for (auto &node1 : combined_set[combined_set_id_map[eq_class_id1]]) {
					for (auto &node2 : combined_set[combined_set_id_map[eq_class_id2]]) {
						std::pair<uint64_t, uint64_t> p = {node1, node2};
						if (node1 > node2) {
							p = {node2, node1};
						}
						if (row_vec_dists[p] > ROW_DIST_THRESHOLD) {
							in_a_set_flag = false;
							break;
						}
					}
					if (in_a_set_flag == false) {
						break;
					}
				}
				if (in_a_set_flag) {
					for (auto &node2 : combined_set[combined_set_id_map[eq_class_id2]]) {
						combined_set[combined_set_id_map[eq_class_id1]].push_back(node2);
						combined_set_id_map[node2] = combined_set_id_map[eq_class_id1];
					}
					combined_set[combined_set_id_map[eq_class_id2]] = {};
				}
			}
		}
	}

	std::cout << "Finishes combine row vectors" << std::endl;
	// Generate compressed iv_buffer
	uint64_t new_eqclass_id = 0;
	IntVector new_iv_buffer(mantis::NUM_IV_BUFFER * num_samples, 0);
	std::unordered_map<uint64_t, uint64_t> prev_eqclassid_to_new_eqclassid;
	// SECOND_MAX_LL_INTEGER means the row is compressed, MAX_LL_INTEGER means the row form a set with one element
	for (uint64_t i = 0; i < combined_set_id_map.size(); i++) {
		if (combined_set_id_map[i] != MAX_LL_INTEGER && combined_set_id_map[i] != SECOND_MAX_LL_INTEGER) {
			std::vector<uint64_t> row_set = combined_set[combined_set_id_map[i]];
			std::vector<uint64_t> new_row_vec(num_samples, 0);
			// set the row vector of a set with more than to two element with largest count for every sample
			for (uint64_t j = 0; j < num_samples; j++) {
				for (auto &row : row_set) {
					if (iv_buffer[row * num_samples + j] > new_row_vec[i]) {
						new_row_vec[j] = iv_buffer[row * num_samples + i];
					}
				}
			}
			for (uint64_t j = 0; j < num_samples; j++) {
				new_iv_buffer[new_eqclass_id + j] = new_row_vec[j];
			}

			for (auto &row : row_set) {
				combined_set_id_map[row] = SECOND_MAX_LL_INTEGER;
				prev_eqclassid_to_new_eqclassid[row] = new_eqclass_id;
			}
			new_eqclass_id++;
		}
		// set only with one element, the row vector is same
		else if (combined_set_id_map[i] == MAX_LL_INTEGER) {
			for (uint64_t j = 0; j < num_samples; j++) {
				new_iv_buffer[new_eqclass_id + j] = iv_buffer[i * num_samples + j];
			}
			combined_set_id_map[i] == SECOND_MAX_LL_INTEGER;
			prev_eqclassid_to_new_eqclassid[i] = new_eqclass_id;
			new_eqclass_id++;
		}
	}
	std::cout << "Finishes generate compressed iv_buffer" << std::endl;
	iv_buffer = new_iv_buffer;

	// todo: update the previoud eqclass id to new eqclass id
}

*/

/*
template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::update_eqclassid_of_dbg(std::unordered_map<int, int> &prev_eqclassid_to_new_eqclassid) {
	struct Iterator {
		QFi qfi;
		typename key_obj::kmer_t kmer;
		uint64_t count;
		uint32_t id;
		bool do_madvice{false};
		Iterator(uint32_t id, const QF* cqf, bool flag): id(id), do_madvice(flag)
		{
			if (qf_iterator_from_position(cqf, &qfi, 0) != QFI_INVALID) {
				get_key();
        		if (do_madvice) {
          			qfi_initial_madvise(&qfi);
				}
      		}
		}
		bool next() {
			if (do_madvice) {
				if (qfi_next_madvise(&qfi) == QFI_INVALID) return false;
			} else {
				if (qfi_next(&qfi) == QFI_INVALID) return false;
			}
			get_key();
			return true;
		}
		bool end() const {
			return qfi_end(&qfi);
		}
		bool operator>(const Iterator& rhs) const {
			return key() > rhs.key();
		}
		const typename key_obj::kmer_t& key() const { return kmer; }
		uint64_t get_count() const {return count;}
		private:
		void get_key() {
			uint64_t value;
			qfi_get_hash(&qfi, &kmer, &value, &count);
		}
	};
}
*/

int main() {
    test_compress2();
    return 0;
}