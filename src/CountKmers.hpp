/*
 * CountKmers.hpp
 *
 *  Created on: Feb. 8, 2024
 *      Author: cjustin
 */

#ifndef SRC_COUNTKMERS_HPP_
#define SRC_COUNTKMERS_HPP_

#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <zlib.h>
#include <omp.h>

#include "vendor/tsl/robin_set.h"
#include "vendor/concurrentqueue.h"
#include "Options.h"

#include "vendor/kseq_util.h"
#include "khash/KseqHashIterator.hpp"

#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include "vendor/kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;

class CountKmers {
public:
	CountKmers(const string &kmerFile) :
			m_totalCounts(0), m_matchCounts(0) {
		// TODO Auto-generated constructor stub
		initCountsHash(kmerFile);
		if (m_kmers.size() == 0) {
			cerr << "No k-mers of k" << opt::k << " found in " << kmerFile
					<< endl;
			exit(1);
		}
	}

	void computeCounts(const vector<string> &filenames) {
		if (filenames.size() < opt::threads) {
			computeCountsProducerConsumer(filenames);
		} else {
#pragma omp parallel for
			for (unsigned i = 0; i < filenames.size(); ++i) {
				gzFile fp;
				fp = gzopen(filenames[i].c_str(), "r");
				if (fp == Z_NULL) {
#pragma omp critical (stderr)
					{
						std::cerr << "file " << filenames[i]
								<< " cannot be opened" << std::endl;
					}
					exit(1);
				} else if (opt::verbose) {
#pragma omp critical (stderr)
					{
						std::cerr << "Opening " << filenames[i] << std::endl;
					}
				}
				//read in seq
				kseq_t *seq = kseq_init(fp);
				int l = kseq_read(seq);
				while (l >= 0) {
					incrementCount(seq);
					l = kseq_read(seq);
				}
				kseq_destroy(seq);
				gzclose(fp);
			}
		}
	}

	void printInfo() {
		cout << "Total k-mer Counts: " << m_totalCounts << endl;
		cout << "K-mers Matching Counts: " << m_matchCounts << endl;
		cout << "K-mers in initial set: " << m_kmers.size() << endl;
		cout << "Error Rate:" << computeErrorRate(m_totalCounts, m_kmers.size(), m_matchCounts) << endl;
		cout << "Coverage:" << (double(m_matchCounts)/double(m_kmers.size())) << endl;
	}

	~CountKmers() {
		// TODO Auto-generated destructor stub
	}

private:
	tsl::robin_set<uint64_t> m_kmers;
	uint64_t m_totalCounts;
	uint64_t m_matchCounts;

	void initCountsHash(const string &kmerFile){
		gzFile fp = gzopen(kmerFile.c_str(), "r");
		tsl::robin_set<uint64_t> dupes;
		if (fp == Z_NULL) {
#pragma omp critical (stderr)
			{
				std::cerr << "file " << kmerFile.c_str() << " cannot be opened"
						<< std::endl;
			}
			exit(1);
		} else if (opt::verbose) {
#pragma omp critical (stderr)
			{
				std::cerr << "Opening " << kmerFile.c_str() << std::endl;
			}
		}
		kseq_t *seq = kseq_init(fp);
		int l = kseq_read(seq);
		size_t entryNum = 0;
		while (l >= 0) {
			//k-merize and insert
			for (KseqHashIterator itr(seq->seq.s, seq->seq.l, opt::k);
					itr != itr.end(); ++itr) {
				uint64_t hv = *itr;
				m_kmers.insert(hv);
			}
			l = kseq_read(seq);
			entryNum++;
		}
		kseq_destroy(seq);
		gzclose(fp);
	}

	//use only if threads > number of files
	void computeCountsProducerConsumer(const vector<string> &filenames) {
		const static size_t s_bulkSize = 1024;
		if (opt::threads <= filenames.size()) {
			//not enough threads to saturate
			computeCounts(filenames);
		} else {
			uint64_t numReads = 0, processedCount = 0;

			moodycamel::ConcurrentQueue<kseq_t> workQueue(
					opt::threads * s_bulkSize);
			moodycamel::ConcurrentQueue<kseq_t> recycleQueue(
					opt::threads * s_bulkSize * 2);
			bool good = true;
			typedef std::vector<kseq_t>::iterator iter_t;

			//fill recycleQueue with empty objects
			{
				std::vector<kseq_t> buffer(opt::threads * s_bulkSize * 2,
						kseq_t());
				recycleQueue.enqueue_bulk(
						std::move_iterator<iter_t>(buffer.begin()),
						buffer.size());
			}

#pragma omp parallel
			{
				std::vector<kseq_t> readBuffer(s_bulkSize);
				string outBuffer;
				if (unsigned(omp_get_thread_num()) < filenames.size()) {
					//file reading init
					gzFile fp;
					fp = gzopen(filenames.at(omp_get_thread_num()).c_str(),
							"r");
					kseq_t *seq = kseq_init(fp);

					//per thread token
					moodycamel::ProducerToken ptok(workQueue);

					//tokens for recycle queue
					moodycamel::ConsumerToken rctok(recycleQueue);
					moodycamel::ProducerToken rptok(recycleQueue);

					unsigned dequeueSize = recycleQueue.try_dequeue_bulk(rctok,
							std::move_iterator<iter_t>(readBuffer.begin()),
							s_bulkSize);
					while (dequeueSize == 0) {
						dequeueSize = recycleQueue.try_dequeue_bulk(rctok,
								std::move_iterator<iter_t>(readBuffer.begin()),
								s_bulkSize);
					}

					unsigned size = 0;
					while (kseq_read(seq) >= 0) {
						cpy_kseq(&readBuffer[size++], seq);
						if (dequeueSize == size) {
							//try to insert, if cannot queue is full
							while (!workQueue.try_enqueue_bulk(ptok,
									std::move_iterator<iter_t>(
											readBuffer.begin()), size)) {
								//try to work
								if (kseq_read(seq) >= 0) {
									//------------------------WORK CODE START---------------------------------------
									incrementCount(seq);
									//------------------------WORK CODE END-----------------------------------------
								} else {
									goto fileEmpty;
								}
							}
							//reset buffer
							dequeueSize = recycleQueue.try_dequeue_bulk(rctok,
									std::move_iterator<iter_t>(
											readBuffer.begin()), s_bulkSize);
							while (dequeueSize == 0) {
								//try to work
								if (kseq_read(seq) >= 0) {
									//------------------------WORK CODE START---------------------------------------
									incrementCount(seq);
									//------------------------WORK CODE END-----------------------------------------
								} else {
									goto fileEmpty;
								}
								dequeueSize = recycleQueue.try_dequeue_bulk(
										rctok,
										std::move_iterator<iter_t>(
												readBuffer.begin()),
										s_bulkSize);
							}
							size = 0;
						}
					}
					fileEmpty:
					//finish off remaining work
					for (unsigned i = 0; i < size; ++i) {
						//------------------------WORK CODE START---------------------------------------
						incrementCount(seq);
						//------------------------WORK CODE END-----------------------------------------
					}
					assert(
							recycleQueue.enqueue_bulk(rptok,
									std::move_iterator<iter_t>(
											readBuffer.begin()), size));
					if (processedCount < numReads) {
						moodycamel::ConsumerToken ctok(workQueue);
						//join in if others are still not finished
						while (processedCount < numReads) {
							size_t num = workQueue.try_dequeue_bulk(ctok,
									std::move_iterator<iter_t>(
											readBuffer.begin()), s_bulkSize);
							if (num) {
								for (unsigned i = 0; i < num; ++i) {
									//------------------------WORK CODE START---------------------------------------
									incrementCount(seq);
									//------------------------WORK CODE END-----------------------------------------
								}
								assert(
										recycleQueue.enqueue_bulk(rptok,
												std::move_iterator<iter_t>(
														readBuffer.begin()),
												num));
							}
						}
					}
#pragma omp atomic update
					good &= false;
					kseq_destroy(seq);
					gzclose(fp);
				} else {
					moodycamel::ConsumerToken ctok(workQueue);
					moodycamel::ProducerToken rptok(recycleQueue);
					while (good) {
						if (workQueue.size_approx() >= s_bulkSize) {
							size_t num = workQueue.try_dequeue_bulk(ctok,
									std::move_iterator<iter_t>(
											readBuffer.begin()), s_bulkSize);
							if (num) {
								for (unsigned i = 0; i < num; ++i) {
									//------------------------WORK CODE START---------------------------------------
									incrementCount(&readBuffer[i]);
									//------------------------WORK CODE END-----------------------------------------
								}
								assert(
										recycleQueue.enqueue_bulk(rptok,
												std::move_iterator<iter_t>(
														readBuffer.begin()),
												num));
							}
						}
					}
				}
			}
		}
	}

	void incrementCount(kseq_t *seq) {
		for (KseqHashIterator itr(seq->seq.s, seq->seq.l, opt::k); itr != itr.end();
				++itr) {
			if (m_kmers.find(*itr) != m_kmers.end()) {
#pragma omp atomic update
				++m_matchCounts;
			}
#pragma omp atomic update
			++m_totalCounts;
		}
	}

	double computeErrorRate(uint64_t rawCount, uint64_t distinctCount,
			uint64_t matchingCount) const {
		if (rawCount > 0) {
			double expected = double(rawCount) * double(distinctCount)
					/ double(opt::genomeSize);
			return (1.0
					- pow(double(matchingCount) / expected,
							1.0 / double(opt::k)));
		} else {
			return -1.0;
		}
	}
};

#endif /* SRC_COUNTKMERS_HPP_ */

