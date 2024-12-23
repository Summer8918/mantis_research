/*
 * ============================================================================
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *                  Mike Ferdman (), mferdman@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * ============================================================================
 */

#ifndef _COLORED_DBG_H_
#define _COLORED_DBG_H_

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <chrono>

#include <iostream>
#include <vector>
#include <fstream>
#include <inttypes.h>

#include "sparsepp/spp.h"
#include "tsl/sparse_map.h"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/vectors.hpp"
#include "sdsl/coder_elias_gamma.hpp"
#include "sdsl/coder_comma.hpp"
#include "sdsl/coder_fibonacci.hpp"
#include "gqf_cpp.h"
#include "gqf/hashutil.h"
#include "common_types.h"
#include "mantisconfig.hpp"

#define MANTIS_DBG_IN_MEMORY (0x01)
#define MANTIS_DBG_ON_DISK (0x02)
#define INT_VECTOR_BIT_NUM (0x40)
#define ROW_DIST_THRESHOLD 0x8
#define MAX_LL_INTEGER (0xFFFFFFFFFFFFFFFF)
#define SECOND_MAX_LL_INTEGER (0xFFFFFFFFFFFFFFFE)

typedef sdsl::int_vector<INT_VECTOR_BIT_NUM> IntVector;
typedef sdsl::vlc_vector<sdsl::coder::fibonacci> CompressedIntVector;

struct hash128 {
	uint64_t operator()(const __uint128_t& val128) const
	{
		__uint128_t val = val128;
		// Using the same seed as we use in k-mer hashing.
		return MurmurHash64A((void*)&val, sizeof(__uint128_t),
												 2038074743);
	}
};

template <typename Key, typename Value>
using cdbg_bv_map_t = spp::sparse_hash_map<Key, Value, hash128>;

using default_cdbg_bv_map_t = cdbg_bv_map_t<__uint128_t,
			std::pair<uint64_t,uint64_t>>;

template <class qf_obj, class key_obj>
class ColoredDbg {
	public:
		ColoredDbg(std::string& cqf_file, std::vector<std::string>& eqclass_files,
							 std::string& sample_file, int flag);

		ColoredDbg(uint64_t qbits, uint64_t key_bits, enum qf_hashmode hashmode,
							 uint32_t seed, std::string& prefix, uint64_t nqf, int flag);

		void build_sampleid_map(qf_obj *incqfs);

		default_cdbg_bv_map_t&
			construct(qf_obj *incqfs, uint64_t num_kmers);

		void set_console(spdlog::logger* c) { console = c; }
		const CQF<key_obj> *get_cqf(void) const { return &dbg; }
		uint64_t get_num_bitvectors(void) const;
		uint64_t get_num_eqclasses(void) const { return eqclass_map.size(); }
		uint64_t get_num_samples(void) const { return num_samples; }
		std::string get_sample(uint32_t id) const;
		uint32_t seed(void) const { return dbg.seed(); }
		uint64_t range(void) const { return dbg.range(); }

		std::vector<uint64_t>
			find_samples(const mantis::QuerySet& kmers);

        std::unordered_map<uint64_t, std::vector<uint64_t>>
            find_samples(const std::unordered_map<mantis::KmerHash, uint64_t> &uniqueKmers);

		std::unordered_map<mantis::KmerHash, std::vector<uint64_t>> find_samples3(const mantis::QuerySet& kmers);

		uint64_t getEqclassid(uint64_t kmer);

		void serialize();
		
		void reinit(default_cdbg_bv_map_t& map);
		void set_flush_eqclass_dist(void) { flush_eqclass_dis = true; }

	private:
		// returns true if adding this k-mer increased the number of equivalence
		// classes
		// and false otherwise.
		bool add_kmer3(const typename key_obj::kmer_t& hash, const std::vector<uint64_t>&
									vector);
		void add_intvector(const std::vector<uint64_t> & vector, uint64_t eq_id);
		uint64_t get_next_available_id(void);
		//void cv_buffer_serialize();
		void iv_buffer_serialize();
		void compress_iv_buffer();
		void process_row_vec_group(std::vector<uint64_t> &group,
				uint64_t start,
				uint64_t end,
				std::vector<uint64_t> &combined_set_id_map,
				std::unordered_map<uint64_t, std::vector<uint64_t>> &combined_set,
				uint64_t &combined_set_id);
		void update_dbg();
		//void reshuffle_bit_vectors(cdbg_bv_map_t<__uint128_t, std::pair<uint64_t,
		//													 uint64_t>>& map);
		//void reshuffle_count_vectors(cdbg_bv_map_t<__uint128_t, std::pair<uint64_t,
		//													 uint64_t>>& map);
		void reshuffle_int_vectors(cdbg_bv_map_t<__uint128_t, std::pair<uint64_t,
															 uint64_t>>& map);
		std::unordered_map<uint64_t, std::string> sampleid_map;
		// bit_vector --> <eq_class_id, abundance>
		cdbg_bv_map_t<__uint128_t, std::pair<uint64_t, uint64_t>> eqclass_map;
		// previous eq class id -> new eq class id
		std::unordered_map<uint64_t, uint64_t> prev_eqclassid_to_new_eqclassid;
		// the previous largest eq id of lastest compressed row vectors
		uint64_t prevEqId = 0;
		CQF<key_obj> dbg;

		IntVector iv_buffer;
		std::vector<CompressedIntVector> eqclasses;

		std::string prefix;
		uint64_t num_samples;
		uint64_t num_serializations;
		int dbg_alloc_flag;
		bool flush_eqclass_dis{false};
		std::time_t start_time_;
		spdlog::logger* console;
		uint64_t new_eqclass_id = 1;
};

template <class T>
class SampleObject {
	public:
		SampleObject() : obj(), sample_id(), id(0) {}
		SampleObject(T o, std::string& s = std::string(),
								 uint32_t id = 0) : obj(o), sample_id(s), id(id) {}
		SampleObject(const SampleObject& o) : obj(o.obj),
		sample_id(o.sample_id), id(o.id) {}

		T obj;
		std::string sample_id;
		uint32_t id;
};

template <class T>
struct compare {
	bool operator()(const SampleObject<T>& lhs, const SampleObject<T>& rhs)
	{
		return lhs.obj.key > rhs.obj.key;
	}
};

template <class qf_obj, class key_obj>
inline uint64_t ColoredDbg<qf_obj, key_obj>::get_next_available_id(void) {
	return get_num_eqclasses() + 1;
}

template <class qf_obj, class key_obj>
std::string ColoredDbg<qf_obj, key_obj>::get_sample(uint32_t id) const {
	auto it = sampleid_map.find(id);
	if (it == sampleid_map.end())
		return std::string();
	else
		return it->second;
}

template <class qf_obj, class key_obj>
uint64_t ColoredDbg<qf_obj, key_obj>::get_num_bitvectors(void) const {
	uint64_t total = 0;
	for (uint32_t i = 0; i < num_serializations; i++)
		total += eqclasses[i].size();

	return total / num_samples;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj,key_obj>::reshuffle_int_vectors(cdbg_bv_map_t<__uint128_t, 
			std::pair<uint64_t, uint64_t>>& map) {
	IntVector new_iv_buffer(mantis::NUM_IV_BUFFER * num_samples, 0);
	for (auto &it_input : map) {
		auto it_local = eqclass_map.find(it_input.first);

		if (it_local == eqclass_map.end()) {
			console->error("Can't find the vector hash during shuffling");
			exit(1);
		}
		assert(it_local->second.first <= mantis::NUM_IV_BUFFER &&
							it_input.second.first <= mantis::NUM_IV_BUFFER);
		uint64_t src_idx = ((it_local->second.first - 1) * num_samples);
		uint64_t dest_idx = ((it_input.second.first - 1) * num_samples);
		for (uint32_t i = 0; i < num_samples; i++, src_idx++, dest_idx++) {
			if (iv_buffer[src_idx]) {
				new_iv_buffer[dest_idx] = iv_buffer[src_idx];
			}
		}
	}
	iv_buffer = new_iv_buffer;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::reinit(cdbg_bv_map_t<__uint128_t,
										std::pair<uint64_t, uint64_t>>& map) {
	// dbg.reset();
	uint64_t qbits = log2(dbg.numslots());
	uint64_t keybits = dbg.keybits();
	enum qf_hashmode hashmode = dbg.hash_mode();
	uint64_t seed = dbg.seed();
	dbg.delete_file();
	CQF<key_obj>cqf(qbits, keybits, hashmode, seed, prefix + mantis::CQF_FILE);
	dbg = cqf;

	reshuffle_int_vectors(map);
	// Check if the current bit vector buffer is full and needs to be serialized.
	// This happens when the sampling phase fills up the bv buffer.
	if (get_num_eqclasses() % mantis::NUM_IV_BUFFER == 0) {
		// The bit vector buffer is full.
		console->info("Serializing int vector with {} eq classes.",
									get_num_eqclasses());
		compress_iv_buffer();
		iv_buffer_serialize();
	}
	eqclass_map = map;
}

template <class qf_obj, class key_obj>
bool ColoredDbg<qf_obj, key_obj>::add_kmer3(const typename key_obj::kmer_t& key, const std::vector<uint64_t>&
									vector) {
    // A kmer (hash) is seen only once during the merge process.
	// So we insert every kmer in the dbg
	uint64_t eq_id;
	__uint128_t vec_hash = MurmurHash128A((void*)vector.data(),
											vector.size() * 8, 2038074743,
											2038074751);
	auto it = eqclass_map.find(vec_hash);

	bool added_eq_class{false};
	// Find if the eqclass of the kmer is already there.
	// If it is there then increment the abundance.
	// Else create a new eq class.
	if (it == eqclass_map.end()) {
		eq_id = get_next_available_id();
		eqclass_map.emplace(std::piecewise_construct,
												std::forward_as_tuple(vec_hash),
												std::forward_as_tuple(eq_id, 1));
		add_intvector(vector, eq_id - 1);
		added_eq_class = true;
	} else { // eq class is seen before so increment the abundance.
		eq_id = it->second.first;
		// If eq_id is a previous eq id, update the record in eqclass_map and delete the item in prev_eqclassid_to_new_eqclassid
		auto it2 = prev_eqclassid_to_new_eqclassid.find(eq_id);
		if (it2 != prev_eqclassid_to_new_eqclassid.end()) {
			prev_eqclassid_to_new_eqclassid.erase(eq_id);
			eq_id = it2->second;
			it->second.first = eq_id;
		}
		// with standard map
		it->second.second += 1; // update the abundance.
	}

	// check: the k-mer should not already be present.
	uint64_t count = dbg.query(KeyObject(key,0,eq_id), QF_NO_LOCK |
													   QF_KEY_IS_HASH);
	if (count > 0) {
		console->error("in add k-mer 2 K-mer was already present. kmer: {} eqid: {}", key, count);
		exit(1);
	}

	// we use the count to store the eqclass ids
	int ret = dbg.insert(KeyObject(key,0,eq_id), QF_NO_LOCK | QF_KEY_IS_HASH);
	if (ret == QF_NO_SPACE) {
		// This means that auto_resize failed. 
		console->error("The CQF is full and auto resize failed. Please rerun build with a bigger size.");
		exit(1);
	}
	return added_eq_class;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::add_intvector(const std::vector<uint64_t> & vector, uint64_t eq_id) {
	uint64_t start_idx = (eq_id % mantis::NUM_IV_BUFFER) * num_samples;
	//console->info("vector size:{}, num_samples:{}", vector.size(), num_samples);
	//console->info("eq_id:{}", eq_id);
	for (uint32_t i = 0; i < num_samples; i++) {
		iv_buffer[start_idx + i] += vector[i];
		// if (vector[i] != 0) {
        //      console->info("vector i:{}, val: {} ", i, vector[i]);
        // }
	}
}

struct pair_hash {
    template <class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2>& p) const {
        auto hash1 = std::hash<T1>{}(p.first);
		auto hash2 = std::hash<T2>{}(p.second);
		return hash1 ^ (hash2 << 1); // Combine the two hash values
    }
};

// <key, value> : <eq class id : combined set id>
template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::process_row_vec_group(std::vector<uint64_t> &group,
				uint64_t start,
				uint64_t end,
				std::vector<uint64_t> &combined_set_id_map,
				std::unordered_map<uint64_t, std::vector<uint64_t>> &combined_set,
				uint64_t &combined_set_id) {
	// Use a custom hash function for pair hash
	std::unordered_map<std::pair<uint64_t, uint64_t>, int, pair_hash> row_vec_dists;
	std::cout << "Hash group size: " << end - start << std::endl;
	for (int i = start; i < end; i++) {
		uint64_t i_idx = group[i] * num_samples;   // index for row i
		// if (end - start >= 1000 && i % 10 == 0) {
		// 	std::cout << "big loop i: " << i << std::endl;
		// }
		for (uint64_t j = i + 1; j < end; j++) {
			uint64_t j_idx = group[j] * num_samples;   // index for row j
			uint64_t ManhattanDist = 0;
			for (uint64_t k = 0; k < num_samples; k++) {
				if (iv_buffer[j_idx + k] > iv_buffer[i_idx + k]) {
					ManhattanDist += iv_buffer[j_idx + k] - iv_buffer[i_idx + k];
				} else {
					ManhattanDist += iv_buffer[i_idx + k] -  iv_buffer[j_idx + k];
				}
			}
			std::pair<uint64_t, uint64_t> p = {group[i], group[j]};  // always i < j
			if (group[i] > group[j]) {
				p = {group[j], group[i]};
			}
			row_vec_dists[p] = ManhattanDist;
		}
	}

	//std::cout << "Finishes cal manhattan distance for each two row vectors in group" << std::endl;

	for (auto it = row_vec_dists.begin(); it != row_vec_dists.end(); it++) {
		if (it->second <= ROW_DIST_THRESHOLD) {
			uint64_t eq_class_id1 = it->first.first;
			uint64_t eq_class_id2 = it->first.second;
			bool in_a_set_flag = true;
			// row1 and row2 do not belong to any set neither
			// MAX_LL_INTEGER means not assign set id, every set has a id
			if (combined_set_id_map[eq_class_id1] == MAX_LL_INTEGER &&
						combined_set_id_map[eq_class_id2 == MAX_LL_INTEGER]) {
				//std::cout << "row1 and row2 do not belong to any set neither" << std::endl;
				combined_set_id_map[eq_class_id1] = combined_set_id;
                combined_set_id_map[eq_class_id2] = combined_set_id;
				combined_set[combined_set_id].push_back({eq_class_id1});
                combined_set[combined_set_id].push_back({eq_class_id2});
				combined_set_id++;
			}
			// row1 belong to one set, row2 do not belong to any set, check if row2 is close to
			// other vectors in the set row1 belong
			else if (combined_set_id_map[eq_class_id1] != MAX_LL_INTEGER &&
						combined_set_id_map[eq_class_id2] == MAX_LL_INTEGER) {
				in_a_set_flag = true;
				//std::cout << "row1 belong to one set, row2 do not belong to any set" << std::endl;
				uint64_t id = combined_set_id_map[eq_class_id1];
				if (combined_set.find(id) == combined_set.end()) {
					std::cout << "cannot find id:" << id << " in combined_set" << std::endl;
				}
				for (auto &node : combined_set[id]) {
					std::pair<uint64_t, uint64_t> p = {node, eq_class_id2};
					if (node > eq_class_id2) {
						p = {eq_class_id2, node};
					}
					if (row_vec_dists.find(p) == row_vec_dists.end()) {
						std::cout << "cannot find pair in row_vec_dists" << std::endl;
					}
					assert(row_vec_dists.find(p) != row_vec_dists.end());
					if (row_vec_dists[p] > ROW_DIST_THRESHOLD) {
						in_a_set_flag = false;
						break;
					}
				}
				assert(combined_set.find(id) != combined_set.end());
				if (in_a_set_flag) {
					combined_set[id].push_back(eq_class_id2);
					combined_set_id_map[eq_class_id2] = id;
				}
			}
			// row2 belong to one set, row1 do not belong to any set yet, check if row1 is close to the set row2 belong
			else if (combined_set_id_map[eq_class_id2] != MAX_LL_INTEGER
						&& combined_set_id_map[eq_class_id1] == MAX_LL_INTEGER) {
				in_a_set_flag = true;
				//std::cout << "row2 belong to one set, row1 do not belong to any set" << std::endl;
				uint64_t id = combined_set_id_map[eq_class_id2];
				if (combined_set.find(id) == combined_set.end()) {
					std::cout << "cannot find id:" << id << " in combined_set" << std::endl;
				}
				for (auto &node : combined_set[id]) {
					std::pair<uint64_t, uint64_t> p = {node, eq_class_id1};
					if (node > eq_class_id1) {
						p = {eq_class_id1, node};
					}
					if (row_vec_dists.find(p) == row_vec_dists.end()) {
						std::cout << "cannot find pair in row_vec_dists" << std::endl;
					}
					assert(row_vec_dists.find(p) != row_vec_dists.end());
					if (row_vec_dists[p] > ROW_DIST_THRESHOLD) {
						in_a_set_flag = false;
						break;
					}
				}
				assert(combined_set.find(id) != combined_set.end());
				if (in_a_set_flag) {
					combined_set[id].push_back(eq_class_id1);
					combined_set_id_map[eq_class_id1] = id;
				}
			}
			// When eq_class_id1 and eq_class_id2 belongs to different set, try to merge them.
			else {
				in_a_set_flag = true;
				//std::cout << " When eq_class_id1 and eq_class_id2 belongs to different set, try to merge them." << std::endl;
				uint64_t id1 = combined_set_id_map[eq_class_id1];
				uint64_t id2 = combined_set_id_map[eq_class_id2];
				if (combined_set.find(id1) == combined_set.end()) {
					std::cout << "cannot find id1:" << id1 << " in combined_set" << std::endl;
				}
				if (combined_set.find(id2) == combined_set.end()) {
					std::cout << "cannot find id2:" << id2 << " in combined_set" << std::endl;
				}
				assert(combined_set.find(id1) != combined_set.end());
				assert(combined_set.find(id2) != combined_set.end());
				if (combined_set[id1].size() >= 5000 || combined_set[id1].size() >= 5000) {
					std::cout << "combined_set[id1].size() > 10000 && combined_set[id1].size() > 1000" << std::endl;
					continue;
				}
				for (auto &node1 : combined_set[id1]) {
					for (auto &node2 : combined_set[id2]) {
						std::pair<uint64_t, uint64_t> p = {node1, node2};
						if (node1 > node2) {
							p = {node2, node1};
						}
						if (row_vec_dists.find(p) == row_vec_dists.end() || row_vec_dists[p] > ROW_DIST_THRESHOLD) {
							in_a_set_flag = false;
							break;
						}
					}
					if (in_a_set_flag == false) {
						break;
					}
				}
				if (in_a_set_flag) {
					for (auto &node2 : combined_set[id2]) {
						combined_set[id1].push_back(node2);
						combined_set_id_map[node2] = id1;
					}
					combined_set[id2] = {};
				}
			}
		}
	}
	std::cout << "Finishes combine row vectors in group" << std::endl;
}


template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::update_dbg() {
	struct Iterator2 {
		QFi qfi;
		typename key_obj::kmer_t kmer;
		uint64_t count;
		uint64_t value;
		bool do_madvice{false};
		Iterator2(const QF* cqf, bool flag): do_madvice(flag) {
			if (qf_iterator_from_position(cqf, &qfi, 0) != QFI_INVALID) {
				get_key();
				if (do_madvice) {
					qfi_initial_madvise(&qfi);
				}
			}
		}
		const typename key_obj::kmer_t& key() const { return kmer; }
		uint64_t cnt() const {return count;}
		uint64_t val() const {return value;}
		bool next() {
			if (do_madvice) {
				if (qfi_next_madvise(&qfi) == QFI_INVALID) return false;
			} else {
				if (qfi_next(&qfi) == QFI_INVALID) return false;
			}
			get_key();
			return true;
		}
		private:
		void get_key() {
			qfi_get_hash(&qfi, &kmer, &value, &count);
		}
		bool end() const {
			return qfi_end(&qfi);
		}
	};
	//std::cout << "log2(dbg.numslots();)" << std::endl;
	uint64_t qbits = log2(dbg.numslots());
	//std::cout << "dbg.keybits();" << std::endl;
	uint64_t keybits = dbg.keybits();
	//std::cout << " dbg.hash_mode();" << std::endl;
	enum qf_hashmode hashmode = dbg.hash_mode();
	//std::cout << " dbg.seed();" << std::endl;
	uint64_t seed = dbg.seed();
	//std::cout << "start to create newCqf " << std::endl;
	CQF<key_obj> newCqf(qbits, keybits, hashmode, seed, prefix + std::to_string(num_serializations) + "_" + mantis::CQF_FILE);
	//std::cout << "start to create iterator " << std::endl;
	Iterator2 dbg_iter = Iterator2(dbg.get_cqf(), true);
	std::cout << "start to update dbg" << std::endl;
	uint64_t tmp_cnt = 0;
	do {
		//std::cout << "tmp_cnt:" << tmp_cnt << " get key"<< std::endl;
		typename key_obj::kmer_t key = dbg_iter.key();
		//std::cout << "get cnt"<< std::endl;
		uint64_t cnt = dbg_iter.cnt();
		//std::cout << "get val"<< std::endl;
		uint64_t val = 0;
		//std::cout << "insert"<< std::endl;
		// check: the k-mer should not already be present.
		uint64_t count = newCqf.query(KeyObject(key,0,cnt), QF_NO_LOCK |
													   QF_KEY_IS_HASH);
		if (count > 0) {
			console->error("in update dbg K-mer was already present. kmer: {} eqid: {}", key, count);
			exit(1);
		}
		auto it = prev_eqclassid_to_new_eqclassid.find(cnt);
		if (it != prev_eqclassid_to_new_eqclassid.end()) {
			//std::cout << "old eqid" << cnt << " new id:" << it->second << std::endl;
			cnt = it->second;
		}

		int ret = newCqf.insert(KeyObject(key, 0, cnt), QF_NO_LOCK | QF_KEY_IS_HASH);
		if (ret == QF_NO_SPACE) {
			std::cout << "tmp_cnt:" << tmp_cnt << std::endl;
			// This means that auto_resize failed. 
			console->error("The CQF is full and auto resize failed. Please rerun build with a bigger size.");
			exit(1);
		}
		tmp_cnt++;
	} while (dbg_iter.next());
	std::cout << "tmp_cnt:" << tmp_cnt << std::endl;
	std::cout << "finishes updating dbg" << std::endl;
	dbg.delete_file();
	dbg = newCqf;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::compress_iv_buffer() {
	uint64_t rows = iv_buffer.size() / num_samples;    // 20,000,000  
	// When iv buffer is not full
	if (get_num_eqclasses() % mantis::NUM_IV_BUFFER != 0) {
		rows = get_num_eqclasses() - prevEqId; //16,325,803
	}

	std::cout << "Start to compress iv buffer" << std::endl;
	std::cout << "rows:" << rows << std::endl;
	std::unordered_map<uint64_t, std::vector<uint64_t>> normalized_row_vecs_hash;
	for (uint64_t i = 0; i < rows; i++) {
		std::vector<uint64_t> normalized_counter(num_samples + 1, 0);
		uint64_t tmp_sum = 0;
		for (uint64_t j = 0; j < num_samples; j++) {
			normalized_counter[j] = iv_buffer[i * num_samples + j] / (ROW_DIST_THRESHOLD + 1);
			tmp_sum += iv_buffer[i * num_samples + j];
		}
		normalized_counter[num_samples] = tmp_sum / (ROW_DIST_THRESHOLD + 1);
		uint64_t hash_val = MurmurHash64A((void *)normalized_counter.data(),
											normalized_counter.size() * 8, 2038074743);
		normalized_row_vecs_hash[hash_val].push_back(i + 1);
	}

	std::cout << "Finishes hashing row vectors into different group" << std::endl;
	// Form sets for compression. The distance of the any two rows vectors in the set is <= threshold
	//<key, value> : <eq class id : combined set id>
	// MAX_LL_INTEGER means not assign set id, every set has a id
	std::vector<uint64_t> combined_set_id_map(rows, MAX_LL_INTEGER);
	// <key, value> : <set id : row set>
	std::unordered_map<uint64_t, std::vector<uint64_t>> combined_set;
	uint64_t combined_set_id = 0;
	int choice = 1;
	if (choice == 0) {
		uint64_t tmp_cnt = 0;
		for (auto it = normalized_row_vecs_hash.begin(); it != normalized_row_vecs_hash.end(); it++) {
			if (it->second.size() > 1) {
				uint64_t tmp_start = 0, tmp_end = 0, sz = it->second.size();
				while (tmp_end < sz) {
					tmp_start = tmp_end;
					if (tmp_end + 1000 >= sz) {
						tmp_end = sz;
					} else {
						tmp_end += 1000;
					}
					process_row_vec_group(it->second, tmp_start, tmp_end, combined_set_id_map, combined_set, combined_set_id);
				}
			}
			tmp_cnt += it->second.size();
			// std::cout << "completed percent:" << 1.0 * tmp_cnt / rows << std::endl;
		}
	} else if (choice == 1) { // approximate
		uint64_t tmp_cnt = 0;
		for (auto it = normalized_row_vecs_hash.begin(); it != normalized_row_vecs_hash.end(); it++) {
			combined_set[combined_set_id] = it->second;
			for (auto &data : it->second) {
				combined_set_id_map[data - 1] = combined_set_id;
			}
			tmp_cnt += it->second.size();
			combined_set_id++;
			// std::cout << "completed percent:" << 1.0 * tmp_cnt / rows << std::endl;
		}
	}

	std::cout << "Start generating compressed iv_buffer" << std::endl;
	// Generate compressed iv_buffer
	IntVector new_iv_buffer(rows * num_samples, 0);

	new_eqclass_id = 1 + num_serializations * mantis::NUM_IV_BUFFER;

	// SECOND_MAX_LL_INTEGER means the row is compressed
	// MAX_LL_INTEGER means the row form a set with one element
	for (uint64_t i = 0; i < combined_set_id_map.size(); i++) {
		uint64_t id = combined_set_id_map[i];
		if (id != MAX_LL_INTEGER && id != SECOND_MAX_LL_INTEGER  && combined_set.find(id) == combined_set.end()) {
			std::cout << "combined_set.find(id) == combined_set.end()" << std::endl;
			std::cout << "i=" << i << std::endl;
			while (1) {}
		}
		if (id != MAX_LL_INTEGER && id != SECOND_MAX_LL_INTEGER && combined_set[id].size() > 0) {
			std::vector<uint64_t> row_set = combined_set[id];
			std::vector<uint64_t> new_row_vec(num_samples, 0);
			if (new_eqclass_id == 216) {
				std::cout << "row_set.size()" << row_set.size() << std::endl;
			}
			for (auto &row : row_set) {
				for (uint64_t j = 0; j < num_samples; j++) {
					if (iv_buffer[(row - 1) * num_samples + j] > new_row_vec[j]) {
						new_row_vec[j] = iv_buffer[(row - 1) * num_samples + j];
					}
				}
				combined_set_id_map[row - 1] = SECOND_MAX_LL_INTEGER;
				prev_eqclassid_to_new_eqclassid[row + prevEqId] = new_eqclass_id;
			}
			// std::cout << " new row vector for new eqclass id" << new_eqclass_id << std::endl;
			for (uint64_t j = 0; j < num_samples; j++) {
			// 	std::cout << new_row_vec[j] << " ";
				new_iv_buffer[((new_eqclass_id - 1) % mantis::NUM_IV_BUFFER) * num_samples + j] = new_row_vec[j];
			}
			// std::cout << std::endl;
			new_eqclass_id++;
		}
		// set only with one element, the row vector is same
		else if (combined_set_id_map[i] == MAX_LL_INTEGER) {
			for (uint64_t j = 0; j < num_samples; j++) {
				new_iv_buffer[((new_eqclass_id - 1) % mantis::NUM_IV_BUFFER) * num_samples + j] = iv_buffer[i * num_samples + j];
			}
			combined_set_id_map[i] = SECOND_MAX_LL_INTEGER;
			prev_eqclassid_to_new_eqclassid[i + 1 + prevEqId] = new_eqclass_id;
			new_eqclass_id++;
		}
	}
	
	// for (auto it = prev_eqclassid_to_new_eqclassid.begin(); it != prev_eqclassid_to_new_eqclassid.end(); it++) {
	// 	std::cout << "old eqid" << it->first << " new eqid" << it->second << std::endl;
	// }
	iv_buffer = new_iv_buffer;
	std::cout << "Finishes generate compressed iv_buffer 2" << std::endl;
	std::cout << "new_eqclass_id:" << new_eqclass_id << std::endl;
	update_dbg();
	prevEqId += rows;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::iv_buffer_serialize() {
	std::string iv_file(prefix + std::to_string(num_serializations) + "_" +
											mantis::EQCLASS_FILE);
    std::ofstream outFile(iv_file, std::ios::binary);
    if (outFile) {
		CompressedIntVector vlc(iv_buffer);
        size_t size = iv_buffer.size();
		console->info("iv_buffer_serialize, output file name {}, iv_buffer size: {}", iv_file, size);
        sdsl::store_to_file(vlc, iv_file);
        outFile.close();
    } else {
		console->error("Error opening file when serializing cv buffer!");
    }
	if (INT_VECTOR_BIT_NUM == 0) {
		iv_buffer = IntVector(iv_buffer.bit_size(), 0);
	} else {
		iv_buffer = IntVector(iv_buffer.bit_size() / INT_VECTOR_BIT_NUM, 0);
	}
	num_serializations++;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::serialize() {
	// serialize the bv buffer last time if needed
	if (get_num_eqclasses() % mantis::NUM_IV_BUFFER > 0) {
		compress_iv_buffer();
		iv_buffer_serialize();
	}

	// serialize the CQF
	if (dbg_alloc_flag == MANTIS_DBG_IN_MEMORY)
		dbg.serialize(prefix + mantis::CQF_FILE);
	else {
		dbg.close();
		// std::vector<std::string> dbg_files = mantis::fs::GetFilesExt(prefix.c_str(), mantis::CQF_FILE);
		// std::cout << "rename dbg file" << std::endl;
		// assert(dbg_files.size() == 1);
		// std::string newFileName = (prefix + mantis::CQF_FILE);
		// std::rename(dbg_files[0].c_str(), newFileName.c_str());
	}

	//serialize the eq class id map
	std::ofstream opfile(prefix + mantis::SAMPLEID_FILE);
	for (auto sample : sampleid_map) {
		opfile << sample.first << " " << sample.second << std::endl;
	}
	opfile.close();

	if (flush_eqclass_dis) {
		// dump eq class abundance dist for further analysis.
		std::ofstream tmpfile(prefix + "eqclass_dist.lst");
		for (auto sample : eqclass_map)
			tmpfile << sample.second.first << " " << sample.second.second <<
				std::endl;
		tmpfile.close();
	}
}

template <class qf_obj, class key_obj>
std::vector<uint64_t>
ColoredDbg<qf_obj,key_obj>::find_samples(const mantis::QuerySet& kmers) {
	// Find a list of eq classes and the number of kmers that belong those eq
	// classes.
	// find the eq id of kmers.
	std::unordered_map<uint64_t, uint64_t> query_eqclass_map;
	for (auto k : kmers) {
		key_obj key(k, 0, 0);
		uint64_t eqclass = dbg.query(key, 0);
		if (eqclass) {
			query_eqclass_map[eqclass] += 1;
		}
	}

	std::vector<uint64_t> sample_map(num_samples, 0);
	std::vector<std::unordered_map<uint64_t, uint64_t>> sample_kmers_eqid_count(num_samples);

	for (auto it = query_eqclass_map.begin(); it != query_eqclass_map.end();
			 ++it) {
		auto eqclass_id = it->first;
		auto count = it->second;
		// counter starts from 1.
		uint64_t start_idx = (eqclass_id - 1);
		uint64_t bucket_idx = start_idx / mantis::NUM_IV_BUFFER;
		uint64_t bucket_offset = (start_idx % mantis::NUM_IV_BUFFER) * num_samples;
		//std::cout << "num samples" << num_samples << std::endl;
		//std::cout << "For k-mer start" << std::endl;
		//std::cout << "1: bucket_idx:" << bucket_idx << " bucket_offset:" << bucket_offset << std::endl;
		for (uint32_t i = 0; i < num_samples; i++) {
			//std::cout << "i:" << i << std::endl;
			if (eqclasses[bucket_idx][bucket_offset + i] > 0) {
				sample_map[i] += count;
				// std::cout << "eqclass_id:" << eqclass_id << "sample:" << i << std::endl;
				// std::cout << "Occurrence:" << eqclasses[bucket_idx][bucket_offset + i] << std::endl;
				// store the occurence of each eqclass_id of k-mers for each sample file.
				sample_kmers_eqid_count[i][eqclass_id] += eqclasses[bucket_idx][bucket_offset + i];
			}
		}
		//std::cout << "For k-mer end" << std::endl;
	}
	
	std::ofstream opfile("query_sample_file_kmer_eq_id_occurrence_count_res.txt");
	for (int i = 0; i < sample_map.size(); i++) {
		if (sample_map[i] > 0) {
			opfile << "Occurence of eqid in file:" <<get_sample(i) << '\n';
			for (auto it = sample_kmers_eqid_count[i].begin(); it != sample_kmers_eqid_count[i].end(); it++) {
				opfile << "eq_id:" << it->first << " count:" << it->second << '\n';
			}
		}
	}
	return sample_map;
}

template <class qf_obj, class key_obj>
uint64_t
ColoredDbg<qf_obj,key_obj>::getEqclassid(uint64_t kmer) {
	key_obj key(kmer, 0, 0);
	uint64_t eqclass = dbg.query(key, 0);
	return eqclass;
}

template <class qf_obj, class key_obj>
std::unordered_map<mantis::KmerHash, std::vector<uint64_t>>
ColoredDbg<qf_obj,key_obj>::find_samples3(const mantis::QuerySet& kmers) {
	// Find a list of eq classes and the number of kmers that belong those eq
	// classes.
	// find the eq id of kmers.
	
	std::unordered_map<uint64_t, uint64_t> query_eqclass_map;
	std::unordered_map<mantis::KmerHash, uint64_t> kmer_eqid_map;
	for (auto k : kmers) {
		key_obj key(k, 0, 0);
		uint64_t eqclass = dbg.query(key, 0);
		if (eqclass) {
			query_eqclass_map[eqclass] += 1;
			kmer_eqid_map[k] = eqclass;
		}
	}

	std::unordered_map<uint64_t, std::vector<uint64_t>> eqid_res_map;
	for (auto it = query_eqclass_map.begin(); it != query_eqclass_map.end(); ++it) {
		std::vector<uint64_t> sample_kmers_count (num_samples, 0);
		auto eqclass_id = it->first;
		auto count = it->second;
		// counter starts from 1.
		uint64_t start_idx = (eqclass_id - 1);
		uint64_t bucket_idx = start_idx / mantis::NUM_IV_BUFFER;
		uint64_t bucket_offset = (start_idx % mantis::NUM_IV_BUFFER) * num_samples;
		for (uint32_t i = 0; i < num_samples; i++) {
			if (eqclasses[bucket_idx][bucket_offset + i] > 0) {
				//std::cout << "count:" << count << std::endl;
				// store the occurence of each eqclass_id of k-mers for each sample file.
				sample_kmers_count[i] += eqclasses[bucket_idx][bucket_offset + i];
			}
		}
		eqid_res_map[eqclass_id] = sample_kmers_count;
		// if (eqclass_id == 1) {
		// 	std::cout << "sample_kmers_count for eq id 0" << std::endl;
		// 	for (int i = 0; i < num_samples; i++) {
		// 		std::cout << sample_kmers_count[i] << std::endl;
		// 	}
		// }
	}

	std::unordered_map<mantis::KmerHash, std::vector<uint64_t>> res;
	for (auto k : kmers) {
		if (kmer_eqid_map.find(k) != kmer_eqid_map.end() && 
				eqid_res_map.find(kmer_eqid_map[k]) != eqid_res_map.end()) {
			res[k] = eqid_res_map[kmer_eqid_map[k]];
		}
	}
	return res;
}

template <class qf_obj, class key_obj>
std::unordered_map<uint64_t, std::vector<uint64_t>>
ColoredDbg<qf_obj,key_obj>::find_samples(const std::unordered_map<mantis::KmerHash, uint64_t> &uniqueKmers) {
	// Find a list of eq classes and the number of kmers that belong those eq
	// classes.
	std::unordered_map<uint64_t, std::vector<uint64_t>> query_eqclass_map;
	for (auto kv : uniqueKmers) {
		key_obj key(kv.first, 0, 0);
		uint64_t eqclass = dbg.query(key, 0);
		if (eqclass) {
		    kv.second = eqclass;
            query_eqclass_map[eqclass] = std::vector<uint64_t>();
        }
	}

	std::vector<uint64_t> sample_map(num_samples, 0);
	for (auto it = query_eqclass_map.begin(); it != query_eqclass_map.end();
		 ++it) {
		auto eqclass_id = it->first;
		auto &vec = it->second;
		// counter starts from 1.
		uint64_t start_idx = (eqclass_id - 1);
		uint64_t bucket_idx = start_idx / mantis::NUM_IV_BUFFER;
		uint64_t bucket_offset = (start_idx % mantis::NUM_IV_BUFFER) * num_samples;
		//std::cout << "num samples" << num_samples << std::endl;
		//std::cout << "For k-mer start" << std::endl;
		//std::cout << "1 bucket_idx:" << bucket_idx << " bucket_offset:" << bucket_offset << std::endl;
		for (uint32_t i = 0; i < num_samples; i++) {
			if (eqclasses[bucket_idx][bucket_offset+i] > 0) {
				//std::cout << " i:" << i << std::endl;
				vec.push_back(i);
				//std::cout << "eqclass_id:" << eqclass_id << "sample:" << i << std::endl;
			}
		}
		//std::cout << "For k-mer end" << std::endl;
	}
	return query_eqclass_map;
}

template <class qf_obj, class key_obj>
cdbg_bv_map_t<__uint128_t, std::pair<uint64_t, uint64_t>>& ColoredDbg<qf_obj,
	key_obj>::construct(qf_obj *incqfs, uint64_t num_kmers)
{
	uint64_t counter = 0;
	bool is_sampling = (num_kmers < std::numeric_limits<uint64_t>::max());

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
        if (do_madvice)
          qfi_initial_madvise(&qfi);
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

    typename CQF<key_obj>::Iterator walk_behind_iterator;
  
	struct Minheap_PQ {
		void push(const Iterator& obj) {
			c.emplace_back(obj);
			std::push_heap(c.begin(), c.end(), std::greater<Iterator>());
		}
		void pop() {
			std::pop_heap(c.begin(), c.end(), std::greater<Iterator>());
			c.pop_back();
		}
		void replace_top(const Iterator& obj) {
			c.emplace_back(obj);
			pop();
		}
		Iterator& top() { return c.front(); }
		bool empty() const { return c.empty(); }
		private:
		std::vector<Iterator> c;
	};
	Minheap_PQ minheap;

	for (uint32_t i = 0; i < num_samples; i++) {
		Iterator qfi(i, incqfs[i].obj->get_cqf(), true);
		if (qfi.end()) continue;
		minheap.push(qfi);
	}

	while (!minheap.empty()) {
		// BitVector eq_class(num_samples);
		// fix bug: make the size of eq_class2 to be even
		//std::vector<int> eq_class2((num_samples + 1) / 2 * 2, 0);
		std::vector<uint64_t> eq_class2(num_samples, 0);

		KeyObject::kmer_t last_key;
		do {
			Iterator& cur = minheap.top();
			last_key = cur.key();
			// console->info("cur.id {}", cur.id); is 0 and 1
			eq_class2[cur.id] += cur.get_count();
			if (cur.next()) {
				minheap.replace_top(cur);
			}
			else {
				minheap.pop();
			}
		} while(!minheap.empty() && last_key == minheap.top().key());
		
		bool added_eq_class = add_kmer3(last_key, eq_class2);
		++counter;

    	if (counter == 4096) {
      		walk_behind_iterator = dbg.begin();
    	} else if (counter > 4096) {
      		++walk_behind_iterator;
    	}
    
		// Progress tracker
		static uint64_t last_size = 0;
		if (dbg.dist_elts() % 10000000 == 0 &&
				dbg.dist_elts() != last_size) {
			last_size = dbg.dist_elts();
			console->info("Kmers merged: {}  Num eq classes: {}  Total time: {}",
										dbg.dist_elts(), get_num_eqclasses(), time(nullptr) -
										start_time_);
		}

		// Check if the bit vector buffer is full and needs to be serialized.
		if (added_eq_class and (get_num_eqclasses() % mantis::NUM_IV_BUFFER == 0))
		{
			// Check if the process is in the sampling phase.
			if (is_sampling) {
				break;
			} else {
				// The bit vector buffer is full.
				console->info("Serializing int vector with {} eq classes.",
											get_num_eqclasses());
				console->info("new_eqclass_id: {}", new_eqclass_id);
				compress_iv_buffer();
				iv_buffer_serialize();
			}
		} else if (counter > num_kmers) {
			// Check if the sampling phase is finished based on the number of k-mers.
			break;
		}

		//while(!minheap.empty() && minheap.top().end()) minheap.pop();
	}
	return eqclass_map;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::build_sampleid_map(qf_obj *incqfs) {
	for (uint32_t i = 0; i < num_samples; i++) {
		std::pair<uint32_t, std::string> pair(incqfs[i].id, incqfs[i].sample_id);
		sampleid_map.insert(pair);
	}
}

template <class qf_obj, class key_obj>
ColoredDbg<qf_obj, key_obj>::ColoredDbg(uint64_t qbits, uint64_t key_bits,
																				enum qf_hashmode hashmode,
																				uint32_t seed, std::string& prefix,
																				uint64_t nqf, int flag) :
				prefix(prefix), num_samples(nqf),
				num_serializations(0), start_time_(std::time(nullptr)),
				iv_buffer(mantis::NUM_BV_BUFFER * nqf, 0) {

	if (flag == MANTIS_DBG_IN_MEMORY) {
		CQF<key_obj> cqf(qbits, key_bits, hashmode, seed);
		dbg = cqf;
		dbg_alloc_flag = MANTIS_DBG_IN_MEMORY;
	}
	else if (flag == MANTIS_DBG_ON_DISK) {
		CQF<key_obj> cqf(qbits, key_bits, hashmode, seed, prefix + mantis::CQF_FILE);
		dbg = cqf;
		dbg_alloc_flag = MANTIS_DBG_ON_DISK;
	} else {
		ERROR("Wrong Mantis alloc mode.");
		exit(EXIT_FAILURE);
	}
	dbg.set_auto_resize();
}

template <class qf_obj, class key_obj>
ColoredDbg<qf_obj, key_obj>::ColoredDbg(std::string& cqf_file,
										std::vector<std::string>& eqclass_files,
										std::string& sample_file, int flag) : 
										iv_buffer(),
										start_time_(std::time(nullptr)) {
	// Cannot use console->info() to print things here but can use PRINT(), why?
	PRINT("Constructing ColoredDbg.\n");
	num_samples = 0;
	num_serializations = 0;

	if (flag == MANTIS_DBG_IN_MEMORY) {
		CQF<key_obj>cqf(cqf_file, CQF_FREAD);
		dbg = cqf;
		dbg_alloc_flag = MANTIS_DBG_IN_MEMORY;
	} else if (flag == MANTIS_DBG_ON_DISK) {
		CQF<key_obj>cqf(cqf_file, CQF_MMAP);
		dbg = cqf;
		dbg_alloc_flag = MANTIS_DBG_ON_DISK;
	} else {
		ERROR("Wrong Mantis alloc mode.");
		exit(EXIT_FAILURE);
	}

	std::map<int, std::string> sorted_files;
	for (std::string file : eqclass_files) {
		int id = std::stoi(first_part(last_part(file, '/'), '_'));
		sorted_files[id] = file;
	}

	eqclasses.reserve(sorted_files.size());
	CompressedIntVector civ;

	for (auto file : sorted_files) {
		sdsl::load_from_file(civ, file.second); 
		eqclasses.push_back(civ);
		num_serializations++;
	}

	std::ifstream sampleid(sample_file.c_str());
	std::string sample;
	uint32_t id;
	while (sampleid >> id >> sample) {
		std::pair<uint32_t, std::string> pair(id, sample);
		sampleid_map.insert(pair);
		num_samples++;
	}
	sampleid.close();
}

#endif
