/*
 * SingleGenomeUnique.hpp
 *
 *  Created on: Apr 3, 2023
 *      Author: cjustin
 */

#ifndef SRC_SINGLEGENOMEUNIQUE_HPP_
#define SRC_SINGLEGENOMEUNIQUE_HPP_
#include <vector>
#include <string>
#include <omp.h>
#include <stdio.h>
#include <zlib.h>
#include <iostream>

#include "Options.h"

#include "vendor/tsl/robin_map.h"
#include "vendor/tsl/robin_set.h"
#include "vendor/KseqHashIterator.hpp"
//#include "vendor/concurrentqueue.h"
//#include "vendor/kseq_util.h"
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include "vendor/kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;

class SingleGenomeUnique {
public:
	SingleGenomeUnique(const vector<string> &filenames) :
			m_kmers(), m_filenames(filenames) {
#pragma omp parallel for
		for (unsigned i = 0; i < filenames.size(); ++i) {
			gzFile fp;
			fp = gzopen(filenames[i].c_str(), "r");
			if (fp == Z_NULL) {
#pragma omp critical(cerr)
				std::cerr << "file " << filenames[i] << " cannot be opened"
						<< std::endl;
				exit(1);
			} else if (opt::verbose) {
#pragma omp critical(cerr)
				std::cerr << "Opening " << filenames[i] << std::endl;
			}
			//create temp file to count uniqueMGKmerSet
			SGKmerSet tmp;

			//read in seq
			kseq_t *seq = kseq_init(fp);
			int l = kseq_read(seq);
			while (l >= 0) {
				insertRead(*seq, tmp);
				l = kseq_read(seq);
			}
			kseq_destroy(seq);
			gzclose(fp);
			for (SGKmerSet::const_iterator itr = tmp.begin(); itr != tmp.end();
					++itr) {
				if (itr->second == 0) {
					//TODO use locks + buckets to optimize?
#pragma omp critical(m_kmers)
					{
						if (m_kmers.find(itr->first) == m_kmers.end()) {
							m_kmers[itr->first] = unique_ptr<vector<bool>>(
									new vector<bool>(filenames.size(), 0));
						}
						(*m_kmers[itr->first])[i] = true;
					}
				}
			}
		}
	}

	/*
	 * prints summary stats for unique and shared unique kmers
	 */
	void printStats() {
		vector<uint64_t> sgCount(m_filenames.size(), 0);
		vector<uint64_t> mgCount(m_filenames.size(), 0);

		//TODO use locks + buckets to optimize?
		for (MGKmerSet::const_iterator itr = m_kmers.begin();
				itr != m_kmers.end(); ++itr) {
			unsigned sum = 0;
			for (unsigned i = 0; i < m_filenames.size(); ++i) {

				sgCount[i] += unsigned(itr->second->at(i));
				sum += itr->second->at(i);
			}
			mgCount[sum - 1] += 1;
		}

		//genomeID, k, count unique
		ofstream sFH = ofstream(opt::prefix + "_sstats.tsv");

		for (unsigned i = 0; i < m_filenames.size(); ++i){
			sFH << m_filenames[i] << "\t" << opt::k << "\t" << sgCount[i] << endl;
		}
		sFH.close();

		//number of genomes, k, count shared unique
		ofstream mFH = ofstream(opt::prefix + "_mstats.tsv");

		for (unsigned i = 0; i < mgCount.size(); ++i) {
			mFH << (i+1) << "\t" << opt::k << "\t" << mgCount[i]
					<< endl;
		}
		mFH.close();
	}

	~SingleGenomeUnique() {
		// TODO Auto-generated destructor stub
	}

private:
	enum KmerType {
		unique, multiple
	};

	typedef uint16_t GenomeID;
	typedef tsl::robin_map<uint64_t, KmerType> SGKmerSet;
	//TODO use faster bitset?
	typedef tsl::robin_map<uint64_t, unique_ptr<vector<bool>>> MGKmerSet;

	MGKmerSet m_kmers;
	const vector<string> &m_filenames;
//	uint64_t m_totalKmers;


	/*
	 * Inserts k-mers into the set and counts up to 2 occurrences
	 * This is not thread safe, one k-mer set per thread please
	 */
	void insertRead(const kseq_t seq, SGKmerSet &kmerSet) {
		for (KseqHashIterator itr(seq.seq.s, seq.seq.l, opt::k);
				itr != itr.end(); ++itr) {
			if (kmerSet.find(*itr) == kmerSet.end()) {
				kmerSet[*itr] = unique;
			} else {
				kmerSet[*itr] = multiple;
			}
//#pragma omp atomic
//			++m_totalKmers;
		}
	}
}
;

#endif /* SRC_SINGLEGENOMEUNIQUE_HPP_ */
