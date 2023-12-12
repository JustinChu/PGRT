/*
 * khashutil.hpp
 *
 *  Created on: Apr 28, 2023
 *      Author: cjustin
 *
 *  TODO rename to khash.h?
 */

#ifndef SRC_KHASHUTIL_HPP_
#define SRC_KHASHUTIL_HPP_

#include <stdint.h>
#include <string>

using namespace std;

namespace KHashUtil {

const static uint8_t s_ext[256] = { // translates to base extension value
		0, 4, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, //0
		0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //1
		1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //2
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //3
		2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //4
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //5
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //6
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //7
		3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //8
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //9
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //10
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //11
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //12
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //13
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //14
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //15
};

char s_ntToSeqTab[4] = {'A','C','G','T'}; // translate 0123 to ACG

const unsigned char s_seq_nt4_table[256] = { // translate ACGT to 0123
		0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4,
				2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

inline uint64_t hash64(uint64_t key, uint64_t mask) // invertible integer hash function
		{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}


/*
 * Taken from https://github.com/bcgsc/abyss/blob/master/Common/Kmer.cpp
 */
static const uint8_t swapBases[256] = {
	0x00, 0x40, 0x80, 0xc0, 0x10, 0x50, 0x90, 0xd0,
	0x20, 0x60, 0xa0, 0xe0, 0x30, 0x70, 0xb0, 0xf0,
	0x04, 0x44, 0x84, 0xc4, 0x14, 0x54, 0x94, 0xd4,
	0x24, 0x64, 0xa4, 0xe4, 0x34, 0x74, 0xb4, 0xf4,
	0x08, 0x48, 0x88, 0xc8, 0x18, 0x58, 0x98, 0xd8,
	0x28, 0x68, 0xa8, 0xe8, 0x38, 0x78, 0xb8, 0xf8,
	0x0c, 0x4c, 0x8c, 0xcc, 0x1c, 0x5c, 0x9c, 0xdc,
	0x2c, 0x6c, 0xac, 0xec, 0x3c, 0x7c, 0xbc, 0xfc,
	0x01, 0x41, 0x81, 0xc1, 0x11, 0x51, 0x91, 0xd1,
	0x21, 0x61, 0xa1, 0xe1, 0x31, 0x71, 0xb1, 0xf1,
	0x05, 0x45, 0x85, 0xc5, 0x15, 0x55, 0x95, 0xd5,
	0x25, 0x65, 0xa5, 0xe5, 0x35, 0x75, 0xb5, 0xf5,
	0x09, 0x49, 0x89, 0xc9, 0x19, 0x59, 0x99, 0xd9,
	0x29, 0x69, 0xa9, 0xe9, 0x39, 0x79, 0xb9, 0xf9,
	0x0d, 0x4d, 0x8d, 0xcd, 0x1d, 0x5d, 0x9d, 0xdd,
	0x2d, 0x6d, 0xad, 0xed, 0x3d, 0x7d, 0xbd, 0xfd,
	0x02, 0x42, 0x82, 0xc2, 0x12, 0x52, 0x92, 0xd2,
	0x22, 0x62, 0xa2, 0xe2, 0x32, 0x72, 0xb2, 0xf2,
	0x06, 0x46, 0x86, 0xc6, 0x16, 0x56, 0x96, 0xd6,
	0x26, 0x66, 0xa6, 0xe6, 0x36, 0x76, 0xb6, 0xf6,
	0x0a, 0x4a, 0x8a, 0xca, 0x1a, 0x5a, 0x9a, 0xda,
	0x2a, 0x6a, 0xaa, 0xea, 0x3a, 0x7a, 0xba, 0xfa,
	0x0e, 0x4e, 0x8e, 0xce, 0x1e, 0x5e, 0x9e, 0xde,
	0x2e, 0x6e, 0xae, 0xee, 0x3e, 0x7e, 0xbe, 0xfe,
	0x03, 0x43, 0x83, 0xc3, 0x13, 0x53, 0x93, 0xd3,
	0x23, 0x63, 0xa3, 0xe3, 0x33, 0x73, 0xb3, 0xf3,
	0x07, 0x47, 0x87, 0xc7, 0x17, 0x57, 0x97, 0xd7,
	0x27, 0x67, 0xa7, 0xe7, 0x37, 0x77, 0xb7, 0xf7,
	0x0b, 0x4b, 0x8b, 0xcb, 0x1b, 0x5b, 0x9b, 0xdb,
	0x2b, 0x6b, 0xab, 0xeb, 0x3b, 0x7b, 0xbb, 0xfb,
	0x0f, 0x4f, 0x8f, 0xcf, 0x1f, 0x5f, 0x9f, 0xdf,
	0x2f, 0x6f, 0xaf, 0xef, 0x3f, 0x7f, 0xbf, 0xff
};

//TODO possible to alter lookup table to inverse to remove 1 operation
inline uint64_t revComp(uint64_t kmer, unsigned k, uint64_t mask){
	uint64_t rc;
	unsigned char * p = (unsigned char *) &kmer;
	unsigned char * q = (unsigned char *) &rc;
	q[7] = swapBases[p[0]];
	q[6] = swapBases[p[1]];
	q[5] = swapBases[p[2]];
	q[4] = swapBases[p[3]];
	q[3] = swapBases[p[4]];
	q[2] = swapBases[p[5]];
	q[1] = swapBases[p[6]];
	q[0] = swapBases[p[7]];
	return (~rc >> (64-(k*2)) & mask) ;
}

uint8_t s_seqMask = 3; // translate 0123 to ACGT

inline uint64_t getMask(unsigned k){
	return (1ULL << k * 2) - 1;
}

inline uint64_t strToKmer(string kmer, uint64_t mask) {
	uint64_t ntFW, ntRV = 0;
	for (unsigned i = 0; i < kmer.size(); ++i) {
		int c = s_seq_nt4_table[(uint8_t) kmer[i]];
		if (c < 4) { // not an "N" base
			ntFW = (ntFW << 2 | c) & mask;       // forward strand
			ntRV = ntRV >> 2 | (uint64_t) (s_seqMask - c) << (kmer.size() - 1) * 2; // reverse strand
		}
	}
	return(ntFW < ntRV ? ntFW : ntRV);
}

inline string kmerToStr(uint64_t kmer, unsigned k) {
	string kmerStr(k,' ');
	for (unsigned i = 0; i < k; ++i) {
		kmerStr[k-i-1] = s_ntToSeqTab[kmer & s_seqMask];
		kmer = kmer >> 2;
	}
	return(kmerStr);
}

// The inversion of hash_64(). Modified from <https://naml.us/blog/tag/invertible>
inline uint64_t hash64i(uint64_t key, uint64_t mask) {
	uint64_t tmp;

	// Invert key = key + (key << 31)
	tmp = (key - (key << 31));
	key = (key - (tmp << 31)) & mask;

	// Invert key = key ^ (key >> 28)
	tmp = key ^ key >> 28;
	key = key ^ tmp >> 28;

	// Invert key *= 21
	key = (key * 14933078535860113213ull) & mask;

	// Invert key = key ^ (key >> 14)
	tmp = key ^ key >> 14;
	tmp = key ^ tmp >> 14;
	tmp = key ^ tmp >> 14;
	key = key ^ tmp >> 14;

	// Invert key *= 265
	key = (key * 15244667743933553977ull) & mask;

	// Invert key = key ^ (key >> 24)
	tmp = key ^ key >> 24;
	key = key ^ tmp >> 24;

	// Invert key = (~key) + (key << 21)
	tmp = ~key;
	tmp = ~(key - (tmp << 21));
	tmp = ~(key - (tmp << 21));
	key = ~(key - (tmp << 21)) & mask;

	return key;
}
}
#endif /* SRC_KHASHUTIL_HPP_ */
