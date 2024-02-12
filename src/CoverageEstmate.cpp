#include <sstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <iostream>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include "config.h"
#include "src/Options.h"
#include "src/Util.h"
#include "CountKmers.hpp"

#if _OPENMP
# include <omp.h>
#endif

using namespace std;

#define PROGRAM "pgrt-count"

void printVersion() {
	const char VERSION_MESSAGE[] =
	PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
	"Written by Justin Chu <cjustin@ds.dfci.harvard.edu>\n"
	"\n"
	"Copyright 2023 Dana-Farber Cancer Institute\n";
	cerr << VERSION_MESSAGE << endl;
	exit(EXIT_SUCCESS);
}

void printHelpDialog() {
	const string dialog =
			"Usage: " PROGRAM " -s [FASTA] [OPTION]... [FILES...]\n"
			"  -t, --threads = INT    Number of threads to run.[1]\n"
			"  -s, --sites = STR      Fasta of sites to k-merize.\n"
			"                         [required]\n"
			"  -g, --gsize = INT      Genome size for error rate\n"
			"                         estimation. ["+ to_string(opt::genomeSize)+"]\n"
			"  -k, --kmer = INT       k-mer size used. [19]\n"
			"  -h, --help             Display this dialog.\n"
			"  -v, --verbose          Display verbose output.\n"
			"      --version          Print version information.\n";

	cerr << dialog << endl;
	exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[]) {
	//switch statement variable
	int c;

	//control variables
	bool die = false;
	int OPT_VERSION = 0;

	//long form arguments
	static struct option long_options[] = {
			{ "sites", required_argument, NULL, 's' },
			{ "kmer", required_argument, NULL, 'k' },
			{ "help", no_argument, NULL, 'h' },
			{ "version", no_argument,&OPT_VERSION, 1 },
			{ "verbose", no_argument, NULL, 'v' },
			{NULL, 0, NULL, 0 } };

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "s:t:vhk:", long_options,
			&option_index)) != -1) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case 'h': {
			printHelpDialog();
			break;
		}
		case 's': {
			stringstream convert(optarg);
			if (!(convert >> opt::sites)) {
				cerr << "Error - Invalid parameter s: " << optarg << endl;
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
		case 't': {
			stringstream convert(optarg);
			if (!(convert >> opt::threads)) {
				cerr << "Error - Invalid parameter t: " << optarg << endl;
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

	if (opt::k > 32) {
		die = true;
		cerr << "k cannot be greater than 32" << endl;
	}

	vector<string> inputFiles;
	while (optind < argc) {
		inputFiles.emplace_back(argv[optind]);
		assert(Util::fexists(inputFiles.back()));
		optind++;
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

	//read in fasta
	//count files
	//output counts, error rate and coverage, genome size

	CountKmers counter(opt::sites);
	counter.computeCountsProducerConsumer(inputFiles);
	counter.printInfo();
	cerr << "Time: " << omp_get_wtime() - time << " s Memory: " << Util::getRSS()
			<< " kbytes" << endl;
	return 0;
}

