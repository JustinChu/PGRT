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

#include <bitset>

#include "KmerInfo.hpp"
#include "Options.h"

#include "vendor/tsl/robin_map.h"
#include "vendor/tsl/robin_set.h"
#include "khash/KseqHashIterator.hpp"
#include "khash/khashutil.hpp"
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
	SingleGenomeUnique(const vector<string> &filenames) : m_mask(KHashUtil::getMask(opt::k)),
			m_kmers(), m_filenames(filenames), m_unitigs(), m_shift(((opt::k) - 1) * 2) {
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
							m_kmers[itr->first] = KmerInfo();
						}
					}
					m_kmers[itr->first].incrementCount();
				}
			}
		}
	}

	/*
	 * Uses only k-mers in the reference provided
	 */
	SingleGenomeUnique(const vector<string> &filenames, const string &ref) :
			m_mask(KHashUtil::getMask(opt::k)), m_kmers(), m_filenames(
					filenames), m_unitigs(), m_shift(((opt::k) - 1) * 2) {
		//k-merize
		gzFile refFH;
		refFH = gzopen(ref.c_str(), "r");
		if (refFH == Z_NULL) {
			std::cerr << "file " << ref << " cannot be opened" << endl;
			exit(1);
		} else if (opt::verbose) {
			std::cerr << "Opening " << ref << endl;
		}
		kseq_t *seq = kseq_init(refFH);

		if (opt::bed.empty()) {
			//create temp file to count uniqueMGKmerSet
			SGKmerSet tmp;
			//read in seq
			int l = kseq_read(seq);
			while (l >= 0) {
				insertRead(*seq, tmp);
				l = kseq_read(seq);
			}
			for (SGKmerSet::const_iterator itr = tmp.begin(); itr != tmp.end();
					++itr) {
				if (itr->second == 0) {
					if (m_kmers.find(itr->first) == m_kmers.end()) {
						m_kmers[itr->first] = KmerInfo();
					}
//					m_kmers[itr->first].incrementCount();
				}
			}
		} else {
			tsl::robin_map<string, unsigned> chrs;
			vector<unique_ptr<vector<pair<uint64_t, uint64_t>>>> intervals;

			//load in bed file
			ifstream fh(opt::bed.c_str());
			{
				if (!fh.good()) {
					std::cerr << "file " << opt::bed.c_str() << " cannot be opened" << endl;
					exit(1);
				} else if (opt::verbose) {
					std::cerr << "Opening " << opt::bed.c_str() << endl;
				}
				string line;
				while (getline(fh, line)) {
					std::stringstream ss(line);
					string chr;
					ss >> chr;
					uint64_t start, end;
					ss >> start;
					ss >> end;
					if (chrs.find(chr) == chrs.end()) {
						chrs[chr] = intervals.size();
						intervals.emplace_back(
								unique_ptr<vector<pair<uint64_t, uint64_t>>>(
										new vector<pair<uint64_t, uint64_t>>));
					}
					intervals[chrs.at(chr)]->push_back(
							std::make_pair(start, end));
				}
			}

			//create temp file to count uniqueMGKmerSet
			SGKmerSet tmp;
			//read in seq
			int l = kseq_read(seq);
			while (l >= 0) {
				const string &chrName = string(seq->name.s, seq->name.l);
				if (chrs.find(chrName) != chrs.end()) {
					if (opt::verbose) {
						std::cerr << "Processing: " << chrName << endl;
					}
					const vector<pair<uint64_t, uint64_t>> &startEnds =
							*(intervals[chrs.at(chrName)]);
					for (size_t i = 0; i < startEnds.size(); ++i) {
						insertRead(*seq, tmp, startEnds.at(i).first,
								startEnds.at(i).second);
					}
				}
				l = kseq_read(seq);
			}
			for (SGKmerSet::const_iterator itr = tmp.begin(); itr != tmp.end();
					++itr) {
				if (itr->second == 0) {
					if (m_kmers.find(itr->first) == m_kmers.end()) {
						m_kmers[itr->first] = KmerInfo();
					}
//					m_kmers[itr->first].incrementCount();
				}
			}
		}
		kseq_destroy(seq);
		gzclose(refFH);

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
				insertReadExistOnly(*seq, tmp);
				l = kseq_read(seq);
			}
			kseq_destroy(seq);
			gzclose(fp);
			for (SGKmerSet::const_iterator itr = tmp.begin(); itr != tmp.end();
					++itr) {
				if (itr->second == 0) {
					m_kmers[itr->first].incrementCount();
				}
			}
		}
	}

//	/*
//	 * prints summary stats for unique and shared unique kmers
//	 */
//	void printStats() {
//		vector<uint64_t> sgCount(m_filenames.size(), 0);
//		vector<uint64_t> mgCount(m_filenames.size(), 0);
//
//		//TODO use locks + buckets to optimize?
//		for (KMGmerSet::const_iterator itr = m_kmers.begin();
//				itr != m_kmers.end(); ++itr) {
//			unsigned sum = 0;
//			for (unsigned i = 0; i < m_filenames.size(); ++i) {
//
//				sgCount[i] += unsigned(itr->second->at(i));
//				sum += itr->second->at(i);
//			}
//			mgCount[sum - 1] += 1;
//		}
//
//		//genomeID, k, count unique
//		ofstream sFH = ofstream(opt::prefix + "_sstats.tsv");
//
//		for (unsigned i = 0; i < m_filenames.size(); ++i) {
//			sFH << m_filenames[i] << "\t" << opt::k << "\t" << sgCount[i]
//					<< endl;
//		}
//		sFH.close();
//
//		//number of genomes, k, count shared unique
//		ofstream mFH = ofstream(opt::prefix + "_mstats.tsv");
//
//		for (unsigned i = 0; i < mgCount.size(); ++i) {
//			mFH << (i + 1) << "\t" << opt::k << "\t" << mgCount[i] << endl;
//		}
//		mFH.close();
//	}

	void createConnections() {
		if (opt::verbose) {
			std::cerr << "Finding edges" << std::endl;
		}
#pragma omp parallel for
		for (unsigned i = 0; i < m_filenames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_filenames[i].c_str(), "r");
			if (fp == Z_NULL) {
#pragma omp critical(cerr)
				std::cerr << "file " << m_filenames[i] << " cannot be opened"
						<< std::endl;
				exit(1);
			} else if (opt::verbose) {
#pragma omp critical(cerr)
				std::cerr << "Opening " << m_filenames[i] << std::endl;
			}
			//create temp file to count uniqueMGKmerSet
			SGKmerSet tmp;

			//read in seq
			kseq_t *seq = kseq_init(fp);
			int l = kseq_read(seq);
			while (l >= 0) {
				KseqHashIterator itr(seq->seq.s, seq->seq.l, opt::k);
				//get starting pos
				Kmer prevKmer(itr.getPos(), *itr, itr.getDirection());
				for (; itr != itr.end(); ++itr) {
					if (m_kmers[*itr].getCount() == m_filenames.size()) {
						if (prevKmer.m_pos + 1 == itr.getPos()) {
							char addedBase = seq->seq.s[itr.getPos()];
							char removedBase = seq->seq.s[prevKmer.m_pos
									- opt::k];
							m_kmers[prevKmer.m_hashVal].setOutBase(addedBase,
									prevKmer.m_direction);
							m_kmers[*itr].setInBase(removedBase,
									itr.getDirection());
						}
						else {
#pragma omp critical(m_lEdges)
							if (m_lEdges.find(prevKmer.m_hashVal)
									== m_lEdges.end()) {
								m_lEdges[prevKmer.m_hashVal] =
										std::unique_ptr<
												tsl::robin_map<HashedKmer,
														uint8_t>>(
												new tsl::robin_map<HashedKmer,
														uint8_t>());
							}
#pragma omp critical(m_lEdges)
							if (m_lEdges.find(*itr) == m_lEdges.end()) {
								m_lEdges[*itr] =
										std::unique_ptr<
												tsl::robin_map<HashedKmer,
														uint8_t>>(
												new tsl::robin_map<HashedKmer,
														uint8_t>());
							}

							if (prevKmer.m_direction == itr.getDirection()) {
								if (m_lEdges[prevKmer.m_hashVal]->find(*itr)
										!= m_lEdges[prevKmer.m_hashVal]->end()) {
									if ((m_lEdges[prevKmer.m_hashVal]->at(*itr)
											& 3) == Util::s_edgeRV) {
#pragma omp atomic
										(*m_lEdges[prevKmer.m_hashVal])[*itr] |=
												Util::s_edgeBoth;
#pragma omp atomic
										(*m_lEdges[*itr])[prevKmer.m_hashVal] |=
												Util::s_edgeBoth;
									}
								} else {
									//TODO use locks?
#pragma omp critical(m_lEdges)
									{
										(*m_lEdges[prevKmer.m_hashVal])[*itr] |=
												Util::s_edgeFW;
										(*m_lEdges[*itr])[prevKmer.m_hashVal] |=
												Util::s_edgeFW;
									}
								}
							} else {
								if (m_lEdges[prevKmer.m_hashVal]->find(*itr)
										!= m_lEdges[prevKmer.m_hashVal]->end()) {
									if ((m_lEdges[prevKmer.m_hashVal]->at(*itr)
											& 3) == Util::s_edgeFW) {
#pragma omp atomic
										(*m_lEdges[prevKmer.m_hashVal])[*itr] |=
												Util::s_edgeBoth;
#pragma omp atomic
										(*m_lEdges[*itr])[prevKmer.m_hashVal] |=
												Util::s_edgeBoth;
									}
								} else {
									//TODO use locks?
#pragma omp critical(m_lEdges)
									{
										(*m_lEdges[prevKmer.m_hashVal])[*itr] |=
												Util::s_edgeFW;
										(*m_lEdges[*itr])[prevKmer.m_hashVal] |=
												Util::s_edgeFW;
									}
								}
							}
						}
						prevKmer = Kmer(itr.getPos(), *itr, itr.getDirection());
					}
				}
				l = kseq_read(seq);
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
	}

	void printEdgesGV(){
		cout << "digraph D {" << endl;

		//print vertexes
		for(KmerSet::const_iterator itr = m_kmers.begin(); itr != m_kmers.end(); ++itr){
			cout << kmerToStr(KHashUtil::hash64i(itr->first, m_mask)) << endl;
		}

		//print long edges
		for(LongEdgeSet::const_iterator itr = m_lEdges.begin(); itr != m_lEdges.end(); ++itr){
			for(tsl::robin_map<HashedKmer,uint8_t>::const_iterator j = itr->second->begin();j != itr->second->end(); ++j){
				cout << kmerToStr(KHashUtil::hash64i(itr->first, m_mask)) << " -> " << kmerToStr(KHashUtil::hash64i(j->first, m_mask)) << endl;
			}
		}

//		//print short edges
//		for(KmerSet::const_iterator itr = m_kmers.begin(); itr != m_kmers.end(); ++itr){
//
//		}

		cout << "}" << endl;
	}

	void printEdgeStats(){
		vector<unsigned> counts(m_filenames.size(), 0);
		vector<unsigned> lcounts(m_filenames.size(), 0);
		for (LongEdgeSet::const_iterator itr = m_lEdges.begin();
				itr != m_lEdges.end(); ++itr) {
			cout << itr->second->size() << endl;
			++lcounts[itr->second->size()];
		}
//		for (KmerSet::const_iterator itr = m_kmers.begin();
//				itr != m_kmers.end(); ++itr) {
//			cout << getEdgeMultiplicty(itr->second) << endl;
//			++counts[getEdgeMultiplicty(itr->second)];
//		}
		for (unsigned i = 0; i > lcounts.size(); ++i) {
//			if (lcounts.at(i) > 0) {
				cout << i << " " << lcounts.at(i) << endl;
//			}
		}
//		for (unsigned i = 0; i > counts.size(); ++i) {
//			if (counts.at(i) > 0) {
//				cout << i << " " << counts.at(i) << endl;
//			}
//		}
		cout << m_lEdges.size() << endl;
	}

//	void printBasicGraph(){
//		for(auto itr = m_lEdges.begin(); itr != m_lEdges.end(); ++itr){
//			itr->first
//		}
//	}

	/*
	 * Only complete unique k-mers
	 */
	void printSubKmers(){
		if(opt::verbose){
			cerr << "Printing Sub-kmers" << endl;
		}
		string tmp = "";

		size_t kmerID = 0;
		//iterate through k-mers
		for (KmerSet::const_iterator itr = m_kmers.begin();
				itr != m_kmers.end(); ++itr) {
			if (itr->second.getCount() == m_filenames.size()) {
				tmp.clear();
				tmp += '>';
				tmp += to_string(kmerID++);
				tmp += "\n";
				tmp += kmerToStr(itr->first);
				tmp += "\n";
				cout << tmp;
			}
		}
	}

	void genUnitigs(){
		if(opt::verbose){
			cerr << "Generating Unitigs" << endl;
		}
		unsigned unitigID = 0;
		//iterate through k-mers
		for (KmerSet::const_iterator itr = m_kmers.begin();
				itr != m_kmers.end(); ++itr) {
			//find blunted edge (not including long edges)
			if (!itr->second.isTraversed() && itr->second.isBlunt()) {
				uint64_t endKmer = extend(itr->first);
				//TODO output unitig
				m_unitigs[itr->first] = unitigID;
				m_unitigs[endKmer] = unitigID;
				unitigID++;
			}
		}
	}

//	void printUnitigStats(){
//		for(auto itr = m_lEdges.begin(); itr != m_lEdges.end(); ++itr){
//			++lcounts[itr->second->size()];
//		}
//		for (auto itr = m_kmers.begin(); itr != m_kmers.end(); ++itr) {
//			++counts[getEdgeMultiplicty(itr->second)];
//		}
//		for (unsigned i = 0; i > lcounts.size(); ++i) {
//			if (lcounts[i] > 0) {
//				cout << i << " " << lcounts[i] << endl;
//			}
//		}
//		for (unsigned i = 0; i > counts.size(); ++i) {
//			if (counts[i] > 0) {
//				cout << i << " " << counts[i] << endl;
//			}
//		}
//		cout << m_lEdges.size() << endl;
//	}

	~SingleGenomeUnique() {
		// TODO Auto-generated destructor stub
	}

private:
	enum KmerType {
		e_unique, e_multiple
	};

	typedef uint16_t GenomeID;
	typedef uint64_t HashedKmer;
	typedef unsigned UnitigID;

//	struct LongEdgeInfo {
//		LongEdgeInfo(bool dir, bool inOut) :
//				m_direction(dir), m_inOut(inOut) {
//		}
//		bool m_direction; //0 = rv, 1 = fw
//		bool m_inOut; //0 = out, 1 = in
//	};

	struct Kmer {
		Kmer(uint64_t pos, uint64_t hashVal, bool direction) :
				m_pos(pos), m_hashVal(hashVal), m_direction(direction) {
		}
		uint64_t m_pos, m_hashVal;
		bool m_direction;
	};

	struct Unitig{
		HashedKmer start;
		HashedKmer end;
		string seq;
	};

	typedef tsl::robin_map<HashedKmer, KmerType> SGKmerSet;
	//TODO use faster bitset?
	//typedef tsl::robin_map<uint64_t, unique_ptr<vector<bool>>> MGKmerSet;
	typedef tsl::robin_map<HashedKmer, KmerInfo> KmerSet; //iA,iC,iG,iT,oA,oC,oG,oT
	typedef tsl::robin_map<HashedKmer,
			std::unique_ptr<tsl::robin_map<HashedKmer, uint8_t>>> LongEdgeSet; //if distance > 1, first 2 bits for status, rest unused
	typedef tsl::robin_map<HashedKmer, UnitigID> UnitigSet; //only start and end k-mers

	const uint64_t m_mask;

//	MGKmerSet m_kmers;
	KmerSet m_kmers; // counts
	const vector<string> &m_filenames;
	LongEdgeSet m_lEdges;
//	uint64_t m_totalKmers;
	vector<Unitig> m_unitigInfo;
	UnitigSet m_unitigs;
	const uint64_t m_shift;

	unsigned m_multiOutEdgeCount;

	/*
	 * Extends the k-mers given a incorporation base
	 * k-mers are unhashed, base is in numerical form (0123 <-> ACGT)
	 * returns hashed cannonical k-mer
	 */
	uint64_t extendKmer(uint64_t &fw, uint64_t &rv, uint8_t base) {
		fw = (fw << 2 | base) & m_mask;
		rv = rv >> 2 | (uint64_t) (3 - base) << m_shift;
		return KHashUtil::hash64(fw < rv ? fw : rv, m_mask);
	}

	/*
	 * TODO: actually create unitig
	 */
	uint64_t extend(HashedKmer kmer) {
		bool hasNext = true;
		uint64_t fw = KHashUtil::hash64i(kmer, m_mask); //relative to original k-mer
		uint64_t rv = KHashUtil::revComp(fw, opt::k, m_mask);
		bool prevDirection = m_kmers[kmer].anchorDirection(); //1 = FW, 2 = RV
		cout << kmerToStr(fw) << endl;
		cout << kmerToStr(rv) << endl;
		cout << std::bitset<64>(fw) << endl;
		cout << std::bitset<64>(rv) << endl;
		while (hasNext) {
			m_kmers[kmer].setTraversed();
			uint8_t base = m_kmers[kmer].getBaseExtension(prevDirection);
			if (base > 3) {
				cout << m_kmers[kmer].toBitset() << endl;
				cout << kmerToStr(kmer) << endl;
				cout << m_kmers[kmer].getNumBases(prevDirection) << " paths"
						<< endl;
				exit(1);
			} else {
				cout << std::bitset<8>(base) << endl;
				cout << std::bitset<64>(fw) << endl;
				cout << std::bitset<64>(rv) << endl;
				kmer = extendKmer(fw, rv, base);
			}
		}
//		cerr << kmerToStr(kmer) << " " << kmerToStr(prevKmer) << endl;
		return(kmer);
	}

	string kmerToStr(HashedKmer hkmer) {
		return KHashUtil::kmerToStr(KHashUtil::hash64i(hkmer, m_mask), opt::k);
	}

//	/*
//	 * Returns the base extension as a hashed k-mer
//	 */
//	uint64_t getExtension(HashedKmer startKmer, bool &direction) {
//		KmerInfo &info = m_kmers[startKmer];
//		uint64_t mask = KHashUtil::getMask(opt::k);
//		//converted back to kmer
//		uint64_t kmer = KHashUtil::hash64i(startKmer, mask);
//		//shift off end
//		uint64_t newKmer = direction ? kmer >> 2 : kmer << 2;
//		//add new base
//		uint64_t base = info.getBaseExt(direction);
//		//generate revComp version
//
//		//rehash
//		//check if cannonical
//	}

	unsigned getEdgeMultiplicty(KmerInfo status){
		unsigned count = 0;
		for(unsigned i = 0 ; i < 8; ++i ){
			count += status.getEdgeStatus(i) != 1;
		}
		return count;
	}

	/*
	 * Inserts k-mers into the set and counts up to 2 occurrences
	 * This is not thread safe, one k-mer set per thread please
	 */
	void insertRead(const kseq_t seq, SGKmerSet &kmerSet) {
		for (KseqHashIterator itr(seq.seq.s, seq.seq.l, opt::k);
				itr != itr.end(); ++itr) {
			if (kmerSet.find(*itr) == kmerSet.end()) {
				kmerSet[*itr] = e_unique;
			} else {
				kmerSet[*itr] = e_multiple;
			}
//#pragma omp atomic
//			++m_totalKmers;
		}
	}

	/*
	 * Inserts k-mers into the set and counts up to 2 occurrences
	 * This is not thread safe, one k-mer set per thread please
	 */
	void insertRead(const kseq_t seq, SGKmerSet &kmerSet, size_t start, size_t end) {
		for (KseqHashIterator itr(seq.seq.s, seq.seq.l, opt::k, start);
				itr != itr.end(); ++itr) {
			if (itr.getPos() > end) {
				return;
			}
			if (kmerSet.find(*itr) == kmerSet.end()) {
				kmerSet[*itr] = e_unique;
			} else {
				kmerSet[*itr] = e_multiple;
			}
		}
	}


	/*
	 * Inserts k-mers into the set and counts up to 2 occurrences
	 * Only insert if it already exists in k-mers set (increments count)
	 */
	void insertReadExistOnly(const kseq_t seq, SGKmerSet &kmerSet) {
		for (KseqHashIterator itr(seq.seq.s, seq.seq.l, opt::k);
				itr != itr.end(); ++itr) {
			if (m_kmers.find(*itr) != m_kmers.end()) {
				if (kmerSet.find(*itr) == kmerSet.end()) {
					kmerSet[*itr] = e_unique;
				} else {
					kmerSet[*itr] = e_multiple;
				}
			}
		}
	}

	void traverse(){

	}
};

#endif /* SRC_SINGLEGENOMEUNIQUE_HPP_ */
