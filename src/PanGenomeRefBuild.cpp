/*
 * PanGenomeRefBuild.cpp
 *
 *  Created on: Oct. 13, 2020
 *      Author: cjustin
 */

#include <sstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <iostream>
#include <stdlib.h>
#include "config.h"
#include "src/Options.h"
#include "src/Util.h"
#if _OPENMP
# include <omp.h>
#endif

#include "SingleGenomeUnique.hpp"

using namespace std;

#define PROGRAM "pgrt"

void printVersion()
{
	const char VERSION_MESSAGE[] = PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
	"Written by Justin Chu\n"
	"\n"
	"Copyright 2023 Dana-Farber Cancer Institute\n";
	cerr << VERSION_MESSAGE << endl;
	exit(EXIT_SUCCESS);
}

void printHelpDialog()
{
	const string dialog =
	"Usage: pgrt-build [OPTION]... [FILES]...\n"
	"The input is expected to be a set of FASTA files\n\n"
	"  -t, --threads=INT      threads to use. [1]\n"
	"  -p, --prefix=STR       Prefix of output files.\n"
	"  -k, --kmer=INT         Size of k-mer used.[" + to_string(opt::k) + "]\n"
	"  -h, --help             Display this dialog.\n"
	"  -v, --verbose          Display verbose output.\n"
	"      --version          Print version information.\n"
	"\nReport bugs to <cjustin@ds.dfci.harvard.edu>.";

	cerr << dialog << endl;
	exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
	//switch statement variable
	int c;

	//control variables
	bool die = false;
	int OPT_VERSION = 0;

	//long form arguments
	static struct option long_options[] = { {
		"threads", required_argument, NULL, 't' }, {
		"kmer", required_argument, NULL, 'k' }, {
		"help", no_argument, NULL, 'h' }, {
		"version", no_argument, &OPT_VERSION, 1 }, {
		"verbose", no_argument, NULL, 'v' }, {
		NULL, 0, NULL, 0 } };

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "vhk:p:t:", long_options,
			&option_index)) != -1)
	{
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case 'h': {
			printHelpDialog();
			break;
		}
		case 't': {
			stringstream convert(optarg);
			if (!(convert >> opt::threads)) {
				cerr << "Error - Invalid parameter t: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'p': {
			stringstream convert(optarg);
			if (!(convert >> opt::prefix)) {
				cerr << "Error - Invalid parameter p: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'k': {
			stringstream convert(optarg);
			if (!(convert >> opt::k)) {
				cerr << "Error - Invalid parameter k: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'v': {
			opt::verbose++;
			break;
		}
		case '?': {
			die = true;
			break;
		}
		}
	}

#if defined(_OPENMP)
	if (opt::threads > 0)
	omp_set_num_threads(opt::threads);
#endif

	if (OPT_VERSION) {
		printVersion();
	}

	vector<string> inputFiles;
	while (optind < argc) {
		inputFiles.emplace_back(argv[optind]);
		Util::fexists(inputFiles.back());
		optind++;
		//check if file exists
	}
	//Check needed options
	if (inputFiles.size() == 0) {
		cerr << "Error: Need Input File" << endl;
		die = true;
	}

	if (die) {
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}
	double time = omp_get_wtime();
	SingleGenomeUnique gu(inputFiles);
	gu.printStats();
	cerr << "Time: " << omp_get_wtime() - time << " s Memory: " << Util::getRSS()
			<< " kbytes" << endl;
	return 0;
}

