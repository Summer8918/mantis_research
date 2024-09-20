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
#define INT_VECTOR_BIT_NUM (0x20)

typedef std::vector<int> CountVector;

typedef sdsl::int_vector<INT_VECTOR_BIT_NUM> IntVector;
typedef sdsl::vlc_vector<sdsl::coder::fibonacci> CompressedIntVector;

typedef sdsl::bit_vector BitVector;
typedef sdsl::rrr_vector<63> BitVectorRRR;

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

		void serialize();
		
		void reinit(default_cdbg_bv_map_t& map);
		void set_flush_eqclass_dist(void) { flush_eqclass_dis = true; }

	private:
		// returns true if adding this k-mer increased the number of equivalence
		// classes
		// and false otherwise.
		bool add_kmer(const typename key_obj::kmer_t& hash, const BitVector&
									vector);
		bool add_kmer2(const typename key_obj::kmer_t& hash, const CountVector&
									vector);
		bool add_kmer3(const typename key_obj::kmer_t& hash, const std::vector<int>&
									vector);
		void add_countvector(const CountVector& vector, uint64_t eq_id);
		void add_intvector(const std::vector<int> & vector, uint64_t eq_id);
		void add_bitvector(const BitVector& vector, uint64_t eq_id);
		void add_eq_class(BitVector vector, uint64_t id);
		uint64_t get_next_available_id(void);
		void bv_buffer_serialize();
		void cv_buffer_serialize();
		void iv_buffer_serialize();
		//bool deserialize_eqclass(const std::string& filename, CountVector &vec);
		//bool deserialize_eqclass2(const std::string& filename, IntVector &vec);
		void reshuffle_bit_vectors(cdbg_bv_map_t<__uint128_t, std::pair<uint64_t,
															 uint64_t>>& map);
		void reshuffle_count_vectors(cdbg_bv_map_t<__uint128_t, std::pair<uint64_t,
															 uint64_t>>& map);
		void reshuffle_int_vectors(cdbg_bv_map_t<__uint128_t, std::pair<uint64_t,
															 uint64_t>>& map);
		std::unordered_map<uint64_t, std::string> sampleid_map;
		// bit_vector --> <eq_class_id, abundance>
		cdbg_bv_map_t<__uint128_t, std::pair<uint64_t, uint64_t>> eqclass_map;
		CQF<key_obj> dbg;
		BitVector bv_buffer;

		CountVector cv_buffer;
		IntVector iv_buffer;

		//std::vector<BitVectorRRR> eqclasses;
		//std::vector<CountVector> eqclasses;
		std::vector<CompressedIntVector> eqclasses;

		std::string prefix;
		uint64_t num_samples;
		uint64_t num_serializations;
		int dbg_alloc_flag;
		bool flush_eqclass_dis{false};
		std::time_t start_time_;
		spdlog::logger* console;
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
void ColoredDbg<qf_obj,key_obj>::reshuffle_count_vectors(cdbg_bv_map_t<__uint128_t, 
			std::pair<uint64_t, uint64_t>>& map) {
	CountVector new_cv_buffer(mantis::NUM_CV_BUFFER * num_samples);
	for (auto &it_input : map) {
		auto it_local = eqclass_map.find(it_input.first);

		if (it_local == eqclass_map.end()) {
			console->error("Can't find the vector hash during shuffling");
			exit(1);
		}
		assert(it_local->second.first <= mantis::NUM_CV_BUFFER &&
							it_input.second.first <= mantis::NUM_CV_BUFFER);
		uint64_t src_idx = ((it_local->second.first - 1) * num_samples);
		uint64_t dest_idx = ((it_input.second.first - 1) * num_samples);
		for (uint32_t i = 0; i < num_samples; i++, src_idx++, dest_idx++) {
			if (cv_buffer[src_idx]) {
				new_cv_buffer[dest_idx] = 1;
			}
		}
	}
	cv_buffer = new_cv_buffer;
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
			// if (cv_buffer[src_idx]) {
			// 	new_cv_buffer[dest_idx] = 1;
			// }
			if (iv_buffer[src_idx]) {
				new_iv_buffer[dest_idx] = 1;
			}
		}
	}
	iv_buffer = new_iv_buffer;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj,
		 key_obj>::reshuffle_bit_vectors(cdbg_bv_map_t<__uint128_t,
										std::pair<uint64_t, uint64_t>>& map) {
			 BitVector new_bv_buffer(mantis::NUM_BV_BUFFER * num_samples);
			 for (auto& it_input : map) {
				 auto it_local = eqclass_map.find(it_input.first);
				 if (it_local == eqclass_map.end()) {
					 console->error("Can't find the vector hash during shuffling");
					 exit(1);
				 } else {
					 assert(it_local->second.first <= mantis::NUM_BV_BUFFER &&
									it_input.second.first <= mantis::NUM_BV_BUFFER);
					 uint64_t src_idx = ((it_local->second.first - 1) * num_samples);
					 uint64_t dest_idx = ((it_input.second.first - 1) * num_samples);
					 for (uint32_t i = 0; i < num_samples; i++, src_idx++, dest_idx++)
						 if (bv_buffer[src_idx])
							 new_bv_buffer[dest_idx] = 1;
				 }
			 }
			 bv_buffer = new_bv_buffer;
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

	//reshuffle_count_vectors(map);
	reshuffle_int_vectors(map);
	//reshuffle_bit_vectors(map);
	// Check if the current bit vector buffer is full and needs to be serialized.
	// This happens when the sampling phase fills up the bv buffer.
	// if (get_num_eqclasses() % mantis::NUM_CV_BUFFER == 0) {
	if (get_num_eqclasses() % mantis::NUM_IV_BUFFER == 0) {
		// The bit vector buffer is full.
		console->info("Serializing int vector with {} eq classes.",
									get_num_eqclasses());
		//bv_buffer_serialize();
		// cv_buffer_serialize();
		iv_buffer_serialize();
	}
	eqclass_map = map;
}

template <class qf_obj, class key_obj>
bool ColoredDbg<qf_obj, key_obj>::add_kmer2(const typename key_obj::kmer_t& key, const CountVector&
									vector) {
    // A kmer (hash) is seen only once during the merge process.
	// So we insert every kmer in the dbg
	uint64_t eq_id;
	//console->info("vector.size() {}", vector.size());
	// (((vector.size() * 32 + 63) >> 6) << 6) / 8 is to make it align with 64 bits
	__uint128_t vec_hash = MurmurHash128A((void*)vector.data(),
											//(((vector.size() * 32 + 63) >> 6) << 6) / 8, 2038074743,
											(((vector.size() * 32 + 63) >> 6) << 6) / 8, 2038074743,
											2038074751);
	auto it = eqclass_map.find(vec_hash);

	bool added_eq_class{false};
	// Find if the eqclass of the kmer is already there.
	// If it is there then increment the abundance.
	// Else create a new eq class.
	if (it == eqclass_map.end()) {
		// Extracting the lower and higher 64 bits
    	//uint64_t lower_part = (uint64_t)vec_hash;  // Lower 64 bits
    	//uint64_t higher_part = (uint64_t)(vec_hash >> 64);  // Higher 64 bits

    	// Printing the two parts
    	//std::cout << "vec_hash (lower 64 bits): " << lower_part << std::endl;
    	//std::cout << "vec_hash (higher 64 bits): " << higher_part << std::endl;
		eq_id = get_next_available_id();
		eqclass_map.emplace(std::piecewise_construct,
												std::forward_as_tuple(vec_hash),
												std::forward_as_tuple(eq_id, 1));
		add_countvector(vector, eq_id - 1);
		added_eq_class = true;
	} else { // eq class is seen before so increment the abundance.
		eq_id = it->second.first;
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
bool ColoredDbg<qf_obj, key_obj>::add_kmer3(const typename key_obj::kmer_t& key, const std::vector<int>&
									vector) {
    // A kmer (hash) is seen only once during the merge process.
	// So we insert every kmer in the dbg
	uint64_t eq_id;
	//console->info("vector.size() {}", vector.size());
	// (((vector.size() * 32 + 63) >> 6) << 6) / 8 is to make it align with 64 bits
	__uint128_t vec_hash = MurmurHash128A((void*)vector.data(),
											//(((vector.size() * 32 + 63) >> 6) << 6) / 8, 2038074743,
											(((vector.size() * 32 + 63) >> 6) << 6) / 8, 2038074743,
											2038074751);
	auto it = eqclass_map.find(vec_hash);

	bool added_eq_class{false};
	// Find if the eqclass of the kmer is already there.
	// If it is there then increment the abundance.
	// Else create a new eq class.
	if (it == eqclass_map.end()) {
		// Extracting the lower and higher 64 bits
    	//uint64_t lower_part = (uint64_t)vec_hash;  // Lower 64 bits
    	//uint64_t higher_part = (uint64_t)(vec_hash >> 64);  // Higher 64 bits

    	// Printing the two parts
    	//std::cout << "vec_hash (lower 64 bits): " << lower_part << std::endl;
    	//std::cout << "vec_hash (higher 64 bits): " << higher_part << std::endl;
		eq_id = get_next_available_id();
		eqclass_map.emplace(std::piecewise_construct,
												std::forward_as_tuple(vec_hash),
												std::forward_as_tuple(eq_id, 1));
		add_intvector(vector, eq_id - 1);
		added_eq_class = true;
	} else { // eq class is seen before so increment the abundance.
		eq_id = it->second.first;
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
bool ColoredDbg<qf_obj, key_obj>::add_kmer(const typename key_obj::kmer_t&
																					 key, const BitVector& vector) {
	// A kmer (hash) is seen only once during the merge process.
	// So we insert every kmer in the dbg
	uint64_t eq_id;
	__uint128_t vec_hash = MurmurHash128A((void*)vector.data(),
																				vector.capacity()/8, 2038074743,
																				2038074751);

	auto it = eqclass_map.find(vec_hash);
	bool added_eq_class{false};
	// Find if the eqclass of the kmer is already there.
	// If it is there then increment the abundance.
	// Else create a new eq class.
	if (it == eqclass_map.end()) {
		// eq class is seen for the first time.
		eq_id = get_next_available_id();
		eqclass_map.emplace(std::piecewise_construct,
												std::forward_as_tuple(vec_hash),
												std::forward_as_tuple(eq_id, 1));
		add_bitvector(vector, eq_id - 1);
		added_eq_class = true;
	} else { // eq class is seen before so increment the abundance.
		eq_id = it->second.first;
		// with standard map
		it->second.second += 1; // update the abundance.
	}

	// check: the k-mer should not already be present.
	uint64_t count = dbg.query(KeyObject(key,0,eq_id), QF_NO_LOCK |
													   QF_KEY_IS_HASH);
	if (count > 0) {
		console->error("K-mer was already present. kmer: {} eqid: {}", key, count);
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
void ColoredDbg<qf_obj, key_obj>::add_countvector(const CountVector& vector, uint64_t eq_id) {
	uint64_t start_idx = (eq_id % mantis::NUM_CV_BUFFER) * num_samples;
	console->info("vector size:{}, num_samples:{}", vector.size(), num_samples);
	console->info("eq_id:{}", eq_id);
	for (uint32_t i = 0; i < num_samples; i++) {
		cv_buffer[start_idx + i] = vector[i];
		console->info("vector i:{}, val: {} ", i, vector[i]);
	}
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::add_intvector(const std::vector<int> & vector, uint64_t eq_id) {
	uint64_t start_idx = (eq_id % mantis::NUM_IV_BUFFER) * num_samples;
	console->info("vector size:{}, num_samples:{}", vector.size(), num_samples);
	console->info("eq_id:{}", eq_id);
	for (uint32_t i = 0; i < num_samples; i++) {
		iv_buffer[start_idx + i] = vector[i];
		console->info("vector i:{}, val: {} ", i, vector[i]);
	}
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::add_bitvector(const BitVector& vector,
																								uint64_t eq_id) {
	uint64_t start_idx = (eq_id  % mantis::NUM_BV_BUFFER) * num_samples;
	for (uint32_t i = 0; i < num_samples/64*64; i+=64)
		bv_buffer.set_int(start_idx+i, vector.get_int(i, 64), 64);
	if (num_samples%64)
		bv_buffer.set_int(start_idx+num_samples/64*64,
											vector.get_int(num_samples/64*64, num_samples%64),
											num_samples%64);
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::bv_buffer_serialize() {
	BitVector bv_temp(bv_buffer);
	if (get_num_eqclasses() % mantis::NUM_BV_BUFFER > 0) {
		bv_temp.resize((get_num_eqclasses() % mantis::NUM_BV_BUFFER) *
									 num_samples);
	}

	BitVectorRRR final_com_bv(bv_temp);
	std::string bv_file(prefix + std::to_string(num_serializations) + "_" +
											mantis::EQCLASS_FILE);
	sdsl::store_to_file(final_com_bv, bv_file);
	bv_buffer = BitVector(bv_buffer.bit_size());
	num_serializations++;
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
	iv_buffer = IntVector(iv_buffer.bit_size() / INT_VECTOR_BIT_NUM, 0);
	num_serializations++;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::cv_buffer_serialize() {
	std::string cv_file(prefix + std::to_string(num_serializations) + "_" +
											mantis::EQCLASS_FILE);
    std::ofstream outFile(cv_file, std::ios::binary);
    if (outFile) {
        size_t size = cv_buffer.size();
		console->info("cv_buffer_serialize, output file name {}, cv_buffer size: {}", cv_file, cv_buffer.size());
        outFile.write(reinterpret_cast<const char*>(&size), sizeof(size));
        outFile.write(reinterpret_cast<const char*>(cv_buffer.data()), size * sizeof(int));
        outFile.close();
    } else {
		console->error("Error opening file when serializing cv buffer!");
    }
	num_serializations++;
}

template <class qf_obj, class key_obj>
void ColoredDbg<qf_obj, key_obj>::serialize() {
	// serialize the CQF
	if (dbg_alloc_flag == MANTIS_DBG_IN_MEMORY)
		dbg.serialize(prefix + mantis::CQF_FILE);
	else
		dbg.close();

	// serialize the bv buffer last time if needed
	if (get_num_eqclasses() % mantis::NUM_IV_BUFFER > 0)
		iv_buffer_serialize();

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

// template <class qf_obj, class key_obj> 
// bool ColoredDbg<qf_obj,key_obj>::deserialize_eqclass(const std::string& filename, CountVector &vec) {
//     std::ifstream inFile(filename, std::ios::binary);
//     if (inFile) {
//         // size_t size;
//         // inFile.read(reinterpret_cast<char*>(&size), sizeof(size));
//         // vec.resize(size);
//         // inFile.read(reinterpret_cast<char*>(vec.data()), size * sizeof(int));
//         // inFile.close();

//     } else {
//         std::cerr << "Error opening file for reading!" << std::endl;
//         return false;
//     }
//     return true;
// }

template <class qf_obj, class key_obj>
std::vector<uint64_t>
ColoredDbg<qf_obj,key_obj>::find_samples(const mantis::QuerySet& kmers) {
	// Find a list of eq classes and the number of kmers that belong those eq
	// classes.
	std::unordered_map<uint64_t, uint64_t> query_eqclass_map;
	for (auto k : kmers) {
		key_obj key(k, 0, 0);
		uint64_t eqclass = dbg.query(key, 0);
		if (eqclass) {
			query_eqclass_map[eqclass] += 1;
		}
	}

	std::vector<uint64_t> sample_map(num_samples, 0);

	for (auto it = query_eqclass_map.begin(); it != query_eqclass_map.end();
			 ++it) {
		auto eqclass_id = it->first;
		auto count = it->second;
		// counter starts from 1.
		uint64_t start_idx = (eqclass_id - 1);
		// uint64_t bucket_idx = start_idx / mantis::NUM_BV_BUFFER;
		uint64_t bucket_idx = start_idx / mantis::NUM_IV_BUFFER;
		// uint64_t bucket_offset = (start_idx % mantis::NUM_BV_BUFFER) * num_samples;
		uint64_t bucket_offset = (start_idx % mantis::NUM_IV_BUFFER) * num_samples;
		std::cout << "num samples" << num_samples << std::endl;
		std::cout << "For k-mer start" << std::endl;
		// for (uint32_t w = 0; w <= num_samples / 64; w++) {
		// 	uint64_t len = std::min((uint64_t)64, num_samples - w * 64);
		// 	std::cout << "len:" << len << std::endl;
		// 	//CountVector::const_iterator first = eqclasses[bucket_idx].begin() + bucket_offset;
		// 	//CountVector::const_iterator last = eqclasses[bucket_idx].begin() + bucket_offset + len;
		// 	uint64_t wrd = eqclasses[bucket_idx].get_int(bucket_offset, len);
		// 	CountVector wrd(first, last);
		// 	for (uint32_t i = 0, sCntr = w * 64; i < len; i++, sCntr++) {
		// 		//if ((wrd >> i) & 0x01)
		// 		if (wrd[i] & 0x1) {
		// 			sample_map[sCntr] += count;
		// 			std::cout << "eqclass_id:" << eqclass_id << "sample:" << sCntr << std::endl;
		// 		}
		// 	}
		// 	bucket_offset += len;
		// }
		std::cout << "1: bucket_idx:" << bucket_idx << " bucket_offset:" << bucket_offset << std::endl;
		for (uint32_t i = 0; i < num_samples; i++) {
			std::cout << "i:" << i << std::endl;
			if (eqclasses[bucket_idx][bucket_offset + i] > 0) {
				sample_map[i] += count;
				std::cout << "eqclass_id:" << eqclass_id << "sample:" << i << std::endl;
			}
		}
		std::cout << "For k-mer end" << std::endl;
	}
	return sample_map;
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
		//uint64_t bucket_idx = start_idx / mantis::NUM_BV_BUFFER;
		uint64_t bucket_idx = start_idx / mantis::NUM_IV_BUFFER;
		//uint64_t bucket_offset = (start_idx % mantis::NUM_BV_BUFFER) * num_samples;
		uint64_t bucket_offset = (start_idx % mantis::NUM_IV_BUFFER) * num_samples;
		std::cout << "num samples" << num_samples << std::endl;
		std::cout << "For k-mer start" << std::endl;
		// for (uint32_t w = 0; w <= num_samples / 64; w++) {
		// 	uint64_t len = std::min((uint64_t)64, num_samples - w * 64);
		// 	std::cout << "len:" << len << std::endl;
		// 	uint64_t wrd = eqclasses[bucket_idx].get_int(bucket_offset, len);
		// 	//CountVector::const_iterator first = eqclasses[bucket_idx].begin() + bucket_offset;
		// 	//CountVector::const_iterator last = eqclasses[bucket_idx].begin() + bucket_offset + len;
		// 	//CountVector wrd(first, last);
		// 	for (uint32_t i = 0, sCntr = w * 64; i < len; i++, sCntr++) {
		// 		//if ((wrd >> i) & 0x01)
		// 		if (wrd[i] & 0x1) {
		// 			vec.push_back(sCntr);
		// 			std::cout << "eqclass_id:" << eqclass_id << "sample:" << sCntr << std::endl;
		// 		}
				
		// 	}
		// 	bucket_offset += len;
		// }
		std::cout << "1 bucket_idx:" << bucket_idx << " bucket_offset:" << bucket_offset << std::endl;
		for (uint32_t i = 0; i < num_samples; i++) {
			if (eqclasses[bucket_idx][bucket_offset+i] > 0) {
				std::cout << " i:" << i << std::endl;
				vec.push_back(i);
				std::cout << "eqclass_id:" << eqclass_id << "sample:" << i << std::endl;
			}
		}
		std::cout << "For k-mer end" << std::endl;
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
		private:
		void get_key() {
			uint64_t value, count;
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
		CountVector eq_class2((num_samples + 1) / 2 * 2, 0);

		KeyObject::kmer_t last_key;
		do {
			Iterator& cur = minheap.top();
			last_key = cur.key();
			// eq_class[cur.id] = 1;
			// console->info("cur.id {}", cur.id); is 0 and 1
			eq_class2[cur.id] = 1;
			if (cur.next()) {
				minheap.replace_top(cur);
			}
			else {
				minheap.pop();
			}
		} while(!minheap.empty() && last_key == minheap.top().key());
		//bool added_eq_class = add_kmer(last_key, eq_class);
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
				bv_buffer(mantis::NUM_BV_BUFFER * nqf), prefix(prefix), num_samples(nqf),
				num_serializations(0), start_time_(std::time(nullptr)), cv_buffer(mantis::NUM_BV_BUFFER * nqf)
				, iv_buffer(mantis::NUM_BV_BUFFER * nqf, 0) {

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
										std::string& sample_file, int flag) : cv_buffer(),
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
	//BitVectorRRR bv;
	CompressedIntVector civ;
	//CountVector cv;

	for (auto file : sorted_files) {
		// deserialize_eqclass(file.second, cv);
		//sdsl::load_from_file(bv, file.second);
		//eqclasses.push_back(bv);
		sdsl::load_from_file(civ, file.second); 
		//eqclasses.push_back(cv);
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
