#ifndef __MANTIS_COMMON_TYPES__
#define __MANTIS_COMMON_TYPES__

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <cstdint>


namespace mantis {
  using KmerHash = uint64_t;
  using ExperimentID = uint64_t;
  using QuerySet = std::unordered_set<KmerHash>;
  using QuerySets = std::vector<QuerySet>;
  using QueryMap = std::unordered_map<KmerHash, uint64_t>;
  using EqMap = std::unordered_map<KmerHash, std::vector<ExperimentID>>;
  struct BulkQuery {
    QuerySet qs;
     EqMap qmap;
  };


  using QueryResult = std::vector<uint64_t>;//std::unordered_map<uint64_t, uint64_t>;
  using QueryResults = std::vector<QueryResult>;
  using QueryResult2 = std::unordered_map<uint64_t, std::vector<uint64_t>>;
}

#endif //__MANTIS_COMMON_TYPES__
