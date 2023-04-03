/*
 * GenPangenome.hpp
 *
 *  Created on: Oct. 16, 2020
 *      Author: cjustin
 */

#ifndef SRC_GENPANGENOME_HPP_
#define SRC_GENPANGENOME_HPP_
#include <bifrost/ColoredCDBG.hpp>
#include "Options.h"
//#include "kseq_util.h"

class GenPangenome {
public:
	/*
	 * Generate a dbg pangenome graph
	 * Colours indicate sample origin information
	 */
	GenPangenome(const vector<string> &filenames) :
			m_filenames(filenames) {
		ColoredCDBG<> colCDBG(opt::k);
		colCDBG.add()

	};



	//Read in files
	//Gather file IDs -> associate as ID
	//print out ID file
	//Break file k-mers
	//load into bifrost
	//Read file again for path information
	//print out GFA file? with path information according to k-mer overlaps

	//load in reads file

	/*
	 * Destructor
	 */
	virtual ~GenPangenome(){

	};

private:
	const vector<string> &m_filenames;
};

#endif /* SRC_GENPANGENOME_HPP_ */

