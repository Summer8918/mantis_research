/*
 * =====================================================================================
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
 */

#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <unordered_set>
#include <bitset>
#include <cassert>
#include <fstream>

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <openssl/rand.h>

#include "MantisFS.h"
#include "ProgOpts.h"
#include "kmer.h"
#include "coloreddbg.h"
#include "mantisconfig.hpp"
#include "squeakrconfig.h"
#include "validateMantisUtils.h"


#include <stdlib.h>


bool get_cdbg_query_res(std::unordered_map<mantis::KmerHash, std::vector<uint64_t>> &query_res, 
				std::string mantis_index_dir, const mantis::QuerySet & kmers) {
    if (mantis_index_dir.back() != '/') {
		mantis_index_dir += '/';
	}

	// make the output directory if it doesn't exist
	if (!mantis::fs::DirExists(mantis_index_dir.c_str())) {
		std::cout << "Error: mantis directory " <<  mantis_index_dir << " not exist." << std::endl;
		return false;
	}

	// Read the colored dBG
	std::cout << "Reading colored dbg from disk." << std::endl;
	std::vector<std::string> dbg_files = mantis::fs::GetFilesExt(mantis_index_dir.c_str(),
																	mantis::CQF_FILE);
	assert(dbg_files.size() == 1);
	std::string dbg_file = dbg_files[0];
	std::string sample_file(mantis_index_dir + mantis::SAMPLEID_FILE);
	std::vector<std::string> eqclass_files = mantis::fs::GetFilesExt(mantis_index_dir.c_str(),
																	mantis::EQCLASS_FILE);

	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg(dbg_file,
															eqclass_files,
															sample_file,
															MANTIS_DBG_IN_MEMORY);

	std::cout << "Read colored dbg with" << cdbg.get_cqf()->dist_elts() << " k-mers and" 
				<< cdbg.get_num_bitvectors() << "color classes" << std::endl;
	
	query_res = cdbg.find_samples3(kmers);
	return true;
}

bool test_get_cdbg_query_res(const mantis::QuerySet & kmers, std::unordered_map<mantis::KmerHash, std::vector<uint64_t>> res1) {
	std::string mantis_index_dir = "/home/jie/code_base/mantis_research/raw_approx";
	std::unordered_map<mantis::KmerHash, std::vector<uint64_t>> test_q;
	get_cdbg_query_res(test_q, mantis_index_dir, kmers);

	// The size of two result should be same
	assert(test_q.size() == res1.size());
	for (auto it1 = res1.begin(); it1 != res1.end(); it1++) {
		auto it2 = test_q.find(it1->first);
		assert(it2 != test_q.end() && it2->second.size() == it1->second.size());
		for (int i = 0; i < it1->second.size(); i++) {
			assert(it1->second[i] == it2->second[i]);
		}
	}
	std::cout << "Compare res1 and test_q is same" << std::endl;
	return true;
}

void compare_exact_and_approx_query_res(const mantis::QuerySet & kmers) {
	std::string approx_mantis = "/home/jie/code_base/mantis_research/raw_approx";
	std::unordered_map<mantis::KmerHash, std::vector<uint64_t>> approx_res;
	get_cdbg_query_res(approx_res, approx_mantis, kmers);

	std::string exact_mantis = "/home/jie/code_base/mantis_research/raw_exact";
	std::unordered_map<mantis::KmerHash, std::vector<uint64_t>> exact_res;
	get_cdbg_query_res(exact_res, exact_mantis, kmers);

	// The absolute differences of median, mean and std dev between approx and exact mantis
	//  for different kmer row (the occurrence of kmer in different samples)
	
	std::vector<double> medians, means, stdDevs;
	int count = 0;
	for (auto kmer : kmers) {
		auto it1 = approx_res.find(kmer);
		auto it2 = exact_res.find(kmer);
		//std::cout << "kmer " << count << std::endl;
		count++;
		if (it1 == approx_res.end() && it2 == exact_res.end()) {
			//std::cout << "Cann't find kmer" << std::endl;
		} else if (it1 == approx_res.end() || it2 == exact_res.end()) {
			std::cout << "Find kmer in one, not exist in another" << std::endl;
			if (it1 != approx_res.end()) {
				std::cout << "find in appro res" << std::endl;
				for (int i = 0; i < it1->second.size(); i++) {
					if (it1->second[i] != 0) {
						std::cout << "sample index:" << i << " Occurrence:" << it1->second[i] << std::endl;
					}
				}
			}
			if (it2 != exact_res.end()) {
				std::cout << "find in exact res" << std::endl;
			}
		} else if(it2->second.size() == it1->second.size()) {
			std::cout << "Find kmer in exact and approx mantis" << std::endl;
			double approxMean = calculateMean(it1->second);
			double exactMean = calculateMean(it2->second);
			means.push_back(abs(approxMean - exactMean));

			double approxMedian = calculateMedian(it1->second);
			double exactMedian = calculateMedian(it2->second);
			medians.push_back((abs(approxMedian - exactMedian)));

			double approxStdDev = calculateStandardDeviation(it1->second, approxMean);
			double exactStdDev = calculateStandardDeviation(it2->second, exactMean);
			stdDevs.push_back((abs(approxStdDev - exactStdDev)));

			if (approxMean != exactMean) {
				std::cout << "approxMean" << approxMean << " exactMean" << exactMean << std::endl;
			}

			if (approxMedian != exactMedian) {
				std::cout << "approxMedian" << approxMedian << " exactMedian" << exactMedian << std::endl;
			}

			if (approxStdDev != exactStdDev) {
				std::cout << "approxStdDev" << approxStdDev << " exactStdDev" << exactStdDev << std::endl;
			}

			std::vector<int> approxRanks = getRanks(it1->second);
			std::vector<int> exactRanks = getRanks(it2->second);
			assert(approxRanks.size() == exactRanks.size());
			for (int i = 0; i < approxRanks.size(); ++i) {
				if (approxRanks[i] != exactRanks[i]) {
					std::cout << "i:" << i << ": approxRank:" << approxRanks[i] << "exactRanks" 
							<< exactRanks[i] << std::endl;
				}
			}
		} else {
			std::cout << "it2->second.size() == it1->second.size()" << std::endl;
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
	int
validate_main ( ValidateOpts& opt )
{

	spdlog::logger* console = opt.console.get();
	// Read experiment CQFs
	std::ifstream infile(opt.inlist);
  uint64_t num_samples{0};
  if (infile.is_open()) {
    std::string line;
    while (std::getline(infile, line)) { ++num_samples; }
    infile.clear();
    infile.seekg(0, std::ios::beg);
  } else {
    console->error("Input filter list {} does not exist or could not be opened.", opt.inlist);
    std::exit(1);
  }

	std::vector<SampleObject<CQF<KeyObject>*>> inobjects;
  std::vector<CQF<KeyObject>> cqfs;

	// reserve QF structs for input CQFs
  inobjects.reserve(num_samples);
  cqfs.reserve(num_samples);

	// mmap all the input cqfs
	std::string squeakr_file;
	uint32_t nqf = 0;
	uint64_t kmer_size;
	while (infile >> squeakr_file) {
		if (!mantis::fs::FileExists(squeakr_file.c_str())) {
			console->error("Squeakr file {} does not exist.", squeakr_file);
			exit(1);
		}
		squeakr::squeakrconfig config;
		int ret = squeakr::read_config(squeakr_file, &config);
		if (ret == squeakr::SQUEAKR_INVALID_VERSION) {
			console->error("Squeakr index version is invalid. Expected: {} Available: {}",
										 squeakr::INDEX_VERSION, config.version);
			exit(1);
		}
		if (cqfs.size() == 0)
			kmer_size = config.kmer_size;
		else {
			if (kmer_size != config.kmer_size) {
				console->error("Squeakr file {} has a different k-mer size. Expected: {} Available: {}",
											 squeakr_file, kmer_size, config.kmer_size);
				exit(1);
			}
		}
		if (config.cutoff == 1) {
			console->warn("Squeakr file {} is not filtered.", squeakr_file);
		}

    cqfs.emplace_back(squeakr_file, CQF_FREAD);
		std::string sample_id = first_part(first_part(last_part(squeakr_file, '/'),
																									'.'), '_');
		console->info("Reading CQF {} Seed {}", nqf, cqfs[nqf].seed());
		console->info("Sample id {}", sample_id);
		cqfs[nqf].dump_metadata();
    inobjects.emplace_back(&cqfs[nqf], sample_id, nqf);
		if (!cqfs.front().check_similarity(&cqfs.back())) {
			console->error("Passed Squeakr files are not similar.", squeakr_file);
			exit(1);
		}
		nqf++;
	}

	std::string prefix = opt.prefix;
	if (prefix.back() != '/') {
		prefix += '/';
	}
	// make the output directory if it doesn't exist
	if (!mantis::fs::DirExists(prefix.c_str())) {
		mantis::fs::MakeDir(prefix.c_str());
	}

	// Read the colored dBG
	console->info("Reading colored dbg from disk.");
	std::vector<std::string> dbg_files = mantis::fs::GetFilesExt(prefix.c_str(),
																	mantis::CQF_FILE);
	assert(dbg_files.size() == 1);
	std::string dbg_file = dbg_files[0];
	std::string sample_file(prefix + mantis::SAMPLEID_FILE);
	std::vector<std::string> eqclass_files = mantis::fs::GetFilesExt(prefix.c_str(),
																	mantis::EQCLASS_FILE);

	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg(dbg_file,
															eqclass_files,
															sample_file,
															MANTIS_DBG_IN_MEMORY);

	console->info("Read colored dbg with {} k-mers and {} color classes",
								cdbg.get_cqf()->dist_elts(), cdbg.get_num_bitvectors());

	std::string query_file = opt.query_file;
	console->info("Reading query kmers from disk.");
	uint64_t total_kmers = 0;
	std::unordered_map<mantis::KmerHash, uint64_t> _dummy_uniqueKmers;
	mantis::QuerySets multi_kmers = Kmer::parse_kmers(query_file.c_str(),
														kmer_size,
														total_kmers,
														false,
														_dummy_uniqueKmers);
	console->info("Total k-mers to query: {}", total_kmers);

	// Query kmers in each experiment CQF
	// Maintain the fraction of kmers present in each experiment CQF.
	std::vector<std::unordered_map<uint64_t, float>> ground_truth;
	std::vector<std::vector<uint64_t>> cdbg_output;
	bool fail{false};
	for (auto kmers : multi_kmers) {
		std::unordered_map<uint64_t, float> fraction_present;
		for (uint64_t i = 0; i < nqf; i++) {
			for (auto kmer : kmers) {
				KeyObject k(kmer, 0, 0);
				uint64_t count = cqfs[i].query(k, 0);
				if (count > 0)
					fraction_present[inobjects[i].id] += 1;
			}
		}
		// Query kmers in the cdbg
		std::vector<uint64_t> result = cdbg.find_samples(kmers);

		// Validate the cdbg output
		for (uint64_t i = 0; i < nqf; i++) {
			if (fraction_present[i] != result[i]) {
				console->info("Failed for sample: {} original CQF {} cdbg {}",
											inobjects[i].sample_id, fraction_present[i], result[i]);
				fail = true;
				//abort();
			}
		}
		ground_truth.push_back(fraction_present);
		cdbg_output.push_back(result);
	}
	if (fail) {
		console->info("Mantis validation 1 failed!");
	}
	else {
		console->info("Mantis validation 1 passed!");
	}


	fail = false;
	int total_cnt = 0, correct_cnt = 0;
	
	for (auto kmers : multi_kmers) {
		std::unordered_map<mantis::KmerHash, std::vector<uint64_t>> cdbg_output = cdbg.find_samples3(kmers);
		for (auto kmer : kmers) {
			mantis::QuerySet tmpKmerSet;
			tmpKmerSet.insert(kmer);
			KeyObject k(kmer, 0, 0);
			for (uint64_t i = 0; i < nqf; i++) {
				uint64_t count = cqfs[i].query(k, 0);
				uint64_t cdbg_count = 0;
				if (cdbg_output.find(kmer) == cdbg_output.end()) {
					cdbg_count = 0;
				} else if (cdbg_output[kmer].size() < nqf) {
					console->info("cdbg_output[kmer].size() < nqf");
				} else {
					cdbg_count = cdbg_output[kmer][i];
				}

				if (cdbg_count != count) {

					fail = true;
					console->info("Failed for kmer {} in sample: {} original CQF {} cdbg {}",
											kmer, inobjects[i].sample_id, count, cdbg_count);
					uint64_t eq_id = cdbg.getEqclassid(kmer);
					console->info("The eqid of the kmer is {}", eq_id);
				} else {
					correct_cnt += 1;
				}
				total_cnt += 1;
				// else if (count > 0){
				// 	console->info("Passed kmer {} for sample: {} original CQF {} cdbg {}",
				// 							kmer, inobjects[i].sample_id, count, cdbg_output[i]);
				// }
			}

		}
		//test_get_cdbg_query_res(kmers, cdbg_output);
		//compare_exact_and_approx_query_res(kmers);
		console->info("correct ratio: {}", 1.0 * correct_cnt / total_cnt);
		if (fail) {
		    console->info("Mantis validation 3 failed!");
		}
		else {
		    console->info("Mantis validation 3 passed!");
		}
	}

#if 0
	// This is x-axis
	// For one query set
	// Divide aggregate by the total kmers in the query to get the fraction.
	uint64_t cnt = 0;
	std::unordered_map<uint64_t, float> fraction_present = ground_truth[0];
	for (uint64_t i = 0; i < nqf; i++)
		fraction_present[i] = (fraction_present[i] / (float)multi_kmers[cnt++].size()) * 100;

	std::vector<std::vector<std::string>> buckets;
	for (auto it : fraction_present)
		buckets[it.second / 5].push_back(cdbg.get_sample(it.first));
#endif

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

