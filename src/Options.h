/*
 * Options.h
 *
 *  Created on: Oct 14, 2020
 *      Author: cjustin
 */

#ifndef OPTIONS_H
#define OPTIONS_H 1

#include <stdint.h>
#include <string>

using namespace std;

/**
 * Global variables that are mostly constant for the duration of the
 * execution of the program.
 */
namespace opt {
int verbose = 0;
int k = 31;
string prefix = "";
unsigned threads = 1;
string bed = "";
string ref = "";
bool graph = false;
size_t genomeSize;
string sites = "";
}
#endif
