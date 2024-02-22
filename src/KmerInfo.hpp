/*
 * KmerInfo.hpp
 *
 *  Created on: Nov 27, 2023
 *      Author: cjustin
 */

#ifndef SRC_KMERINFO_HPP_
#define SRC_KMERINFO_HPP_

#include <stdio.h>
#include <bitset>

//#include "khash/KseqHashIterator.hpp"
#include "khash/khashutil.hpp"

#include "Options.h"

//TODO turn into full fledged class?
/*
 * relative to cannonical sequence
 */
class KmerInfo {
public:
	KmerInfo() :
			m_val(0), m_count(0), m_traversed(false) {

	}

	uint8_t getEdgeStatus(uint8_t i) {
		return (m_val & (1 << i));
	}

	void setInBase(unsigned char base, bool direction) {
#pragma omp atomic
		m_val |= 1
				<< (direction ?
						KHashUtil::s_seq_nt4_table[base] :
						KHashUtil::s_seq_nt4_table[base] ^ 3);
	}

	void setOutBase(unsigned char base, bool direction) {
#pragma omp atomic
		m_val |= 1
				<< (direction ?
						KHashUtil::s_seq_nt4_table[base] + 4 :
						(KHashUtil::s_seq_nt4_table[base] ^ 3) + 4);
	}

	/*
	 * Assumes blunted edge
	 * Return true if fw has extensions
	 */
	bool anchorDirection() {
		return m_val & s_fwEdges;
	}

	unsigned getNumBases(bool direction) {
		uint8_t v = direction ? m_val & s_fwEdges : m_val & s_rvEdges; // count the number of bits set in v
		unsigned c; // c accumulates the total bits set in v
		for (c = 0; v; c++) {
			v &= v - 1; // clear the least significant bit set
		}
		return c;
	}

	/*
	 * returns numerical base extension (0123<->ACGT)
	 */
	unsigned getBaseExtension(bool direction) {
		uint8_t v = direction ? m_val & s_fwEdges : m_val & s_rvEdges;
		return KHashUtil::s_ext[v];
	}

	/*
	 * Returns the cannonical k-mer extension
	 * Returns the first one seen (assumes only one extension exists)
	 * TODO: prevent unneeded rev comp by keeping both copies
	 */
	uint64_t getExtSingle(uint64_t hashedKmer, bool &direction, uint64_t mask) {
		//get base
		uint64_t kmer = KHashUtil::hash64i(hashedKmer, mask);
		if (direction) {
			kmer = (kmer << 2) & mask;
			for (unsigned i = 0; i > 4; ++i) {
				if ((m_val >> i) & 1) {
					kmer |= i;
					break;
				}
			}
		} else {
			kmer = kmer >> 2;
			for (unsigned i = 0; i > 4; ++i) {
				if ((m_val >> (i + 4)) & 1) {
					kmer |= (i << (opt::k * 2));
					break;
				}
			}
		}
		uint64_t rcKmer = KHashUtil::revComp(kmer, opt::k, mask);
		if (kmer < rcKmer) {
			direction = true;
			return KHashUtil::hash64(kmer, mask);
		}
		direction = false;
		return KHashUtil::hash64(rcKmer, mask);
//			cerr << std::bitset<64>(KHashUtil::hash64i(hashedKmer, mask))
//					<< KHashUtil::kmerToStr(KHashUtil::hash64i(hashedKmer, mask), opt::k)
//					<< endl;
//			cerr << std::bitset<64>(kmer) << KHashUtil::kmerToStr(kmer, opt::k) << endl;
//			cerr << std::bitset<64>(rcKmer) << KHashUtil::kmerToStr(rcKmer, opt::k) << endl;
		//checks for cannonical k-mer
	}

	uint16_t getCount() const {
		return (m_count);
	}

//	/*
//	 * A blunted edge is defined as not having edges on right or left
//	 * 00001111 & val == 0
//	 * 11110000 & val == 0
//	 */
//	bool isBlunt() const {
////		cout << toBitset() << endl;
//		return ((m_val & s_fwEdges) == 0) | ((m_val & s_rvEdges) == 0);
//	}

	void setTraversed() {
#pragma omp atomic
		m_traversed |= true;
	}

	void incrementCount() {
#pragma omp atomic
		++m_count;
	}

	bool isTraversed() const {
		return m_traversed;
	}

	std::bitset<8> toBitset() const {
		return(std::bitset<8>(m_val));
	}

	~KmerInfo() {
	}

private:

	uint8_t m_val; //FW-ACTG RV-ATCG
	uint16_t m_count;
	bool m_traversed;

	static const uint8_t s_fwEdges = 0xF; //11110000
	static const uint8_t s_rvEdges = 0x0F; //00001111
};

#endif /* SRC_KMERINFO_HPP_ */
