/* 
 * File:   Param.cpp
 * Author: david
 * 
 * Created on October 16, 2013, 3:39 PM
 */
/*
 * 
 *  Copyright (c) 2008, 2011, 2012, 2013
 *  David Hernandez, Patrice Francois, Jacques Schrenzel
 * 
 *  This file is part of EDENA.
 *
 *  EDENA is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EDENA is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EDENA.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Param.h"
#include "globalFunc.h"
#include "logWriter.h"
#include <sstream>
#include <cstdlib>
using namespace std;

extern logWriter LOG;
extern bool DEV_INFO;

Param::Param() {
    init();
}

Param::Param(const Param& orig) {
}

Param::~Param() {
}

void Param::init() {
    DEV_INFO = false;
    minOverlap = 0; //computed during overlapping step
    overlapCutoff = 0; //considered during assembly
    truncateTo = 0;
    maxDeadEndLength = 0;
    minNodeCoverage = 0;
    minContigCoverage = 0;
    nThreads = 2;
    minContigSize = 0;
    minNPair = 5;
    minRatio = 0.95;
    maxContigRedundancy = 0;
    shortPeHorizon = 1000;
    longPeHorizon = 15000;
    trimContigEnds = 4; //coverage cutoff for contig ends;
    targetSize = 0;

    prefix = "out";
    singleEndFiles.clear();
    drPairedEndFiles.clear();
    rdPairedEndFiles.clear();
    ovlFile = "";

    contextualCleaning = true;
    cleanGraph = true;
    interactiveShell = false;
    discardNonUsable = true;
    discardSuspicousNodes = true;
    cleanBubbles = true;
    trimRed = true;

    //non documented dev flags
    writeOSG = false; //write the cleanded and condensed version of the OSG
    checkGraph = false; //chech graph for consistency
}

void Param::parseCommandLine(int argc, char**argv) {

    string flag;
    istringstream iss;

    if (argc == 1) {
        edenaUsage(EXIT_SUCCESS);
    }

    for (int argn = 1; argn < argc; argn++) {
        if (argv[argn][0] != '-') {
            cerr << "[err] check command line\n";
            exit(EXIT_FAILURE);
        }

        flag = argv[argn] + 1;

        if (argn > argc-2 && flag != "shell"
                && flag != "dev"
                && flag != "h"
                && flag != "help") {
                cerr << "[err] check command line\n";
                exit(EXIT_FAILURE);
        }

        if (flag == "h" || flag == "help") {
            edenaUsage(EXIT_SUCCESS);
        } else if (flag == "r" || flag == "singleEnd") //one or more read file
        {
            argn++;

            while (argn < argc && argv[argn][0] != '-') {
                this->singleEndFiles.push_back(argv[argn]);
                argn++;
            }
            argn--;
        } else if (flag == "paired" || flag == "DRpairs") {
            //direct-reverse paired read files
            argn++;

            while (argn < argc && argv[argn][0] != '-') {
                this->drPairedEndFiles.push_back(argv[argn]);
                argn++;
            }
            argn--;

        } else if (flag == "matePairs" || flag == "RDpairs") {
            //reverse-direct paired read files

            argn++;

            while (argn < argc && argv[argn][0] != '-') {
                rdPairedEndFiles.push_back(argv[argn]);
                argn++;
            }
            argn--;

        } else if (flag == "e" || flag == "edenaFile") {
            argn++;
            ovlFile = argv[argn];
        } else if (flag == "m" || flag == "overlapCutoff") //for assembling mode
        {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> overlapCutoff;
        } else if (flag == "M" || flag == "minOverlap") //for overlapping mode
        {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> minOverlap;
        } else if (flag == "minNodeCov") { //undocumented
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> minNodeCoverage;
        } else if (flag == "minCoverage") //minimum coverage required for the contigs
        {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> minContigCoverage;
        } else if (flag == "c" || flag == "minContigSize") //min contig size
        {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> minContigSize;
        } else if (flag == "p" || flag == "prefix") // prefix
        {
            argn++;
            prefix = argv[argn];
        } else if (flag == "trim") //low covered contig ends trimming
        {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> trimContigEnds;
            if (trimContigEnds<0)
                trimContigEnds=0;
        } else if (flag == "d" || flag == "deadEnds") //for dead-ends removal
        {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> maxDeadEndLength;
        } else if (flag == "targetSize") {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> targetSize;
        } else if (flag == "discardNonUsable") {
            argn++;
            flag = argv[argn];
            if (flag == "no")
                discardNonUsable = false;
            else if (flag == "yes")
                discardNonUsable = true;
            else {
                cerr << "\n-discardNonUsable: possible values are \"yes\" and \"no\"\n";
                exit(EXIT_FAILURE);
            }
        } else if (flag == "discardSuspicious") { //currently undocumented
            argn++;
            flag = argv[argn];
            if (flag == "no")
                discardSuspicousNodes = false;
            else if (flag == "yes")
                discardSuspicousNodes = true;
            else {
                cerr << "\n-discardSuspicious: possible values are \"yes\" and \"no\"\n";
                exit(EXIT_FAILURE);
            }
        } else if (flag == "cc" || flag == "contextualCleaning") {
            argn++;
            flag = argv[argn];
            if (flag == "no")
                contextualCleaning = false;
            else if (flag == "yes")
                contextualCleaning = true;
            else {
                cerr << "\n-contextualCleaning: possible values are \"yes\" and \"no\"\n";
                exit(EXIT_FAILURE);
            }
        } else if (flag == "clearBubbles") {
            argn++;
            flag = argv[argn];
            if (flag == "no")
                cleanBubbles = false;
            else if (flag == "yes")
                cleanBubbles = true;
            else {
                cerr << "\n-cleanBubbles: possible values are \"yes\" and \"no\"\n";
                exit(EXIT_FAILURE);
            }
        } else if (flag == "peHorizon") {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> shortPeHorizon;
            longPeHorizon = shortPeHorizon;
        } else if (flag == "shortPeHorizon" || flag == "sph") {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> shortPeHorizon;
        } else if (flag == "longPeHorizon" || flag == "lph") {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> longPeHorizon;
        } else if (flag == "t" || flag == "truncate") {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> truncateTo;
        } else if (flag == "nThreads") {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> nThreads;
            if (nThreads < 1)
                nThreads = 1;
        } else if (flag == "minNPairs") {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> minNPair;
        } else if (flag == "minRatio") {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> minRatio;
        } else if (flag == "maxRed") {
            argn++;
            iss.clear();
            iss.str(argv[argn]);
            iss >> maxContigRedundancy;
        } else if (flag == "trimRed") {
            argn++;
            flag = argv[argn];
            if (flag == "no")
                trimRed = false;
            else if (flag == "yes")
                trimRed = true;
            else {
                cerr << "\n-trimRed: possible values are \"yes\" and \"no\"\n";
                exit(EXIT_FAILURE);
            }
        } else if (flag == "devInfo") //undocumented
        {
            DEV_INFO = true; //global flag
        } else if (flag == "dev" || flag == "shell") //undocumented
        {
            interactiveShell = true;
        } else if (flag == "v") {
            edenaVersion(cout);
            exit(EXIT_SUCCESS);
        } else if (flag == "wOSG") //write overlaps string graph
        {
            writeOSG = true;
        } else if (flag == "checkGraph") {
            checkGraph = true;
        } else if (flag == "cleanGraph") {
            argn++;
            flag = argv[argn];
            if (flag == "no")
                cleanGraph = false;
            else if (flag == "yes")
                cleanGraph = true;
            else {
                cerr << "\n--cleanGraph: possible values are \"yes\" and \"no\"\n";
                exit(EXIT_FAILURE);
            }
        } else {
            cerr << "-" << flag << " : unknown argument\n";
            exit(EXIT_FAILURE);
        }
    }

    checkParamConsistency();
}

void Param::checkParamConsistency() {

    readsProvided = (singleEndFiles.size() > 0 ||
            drPairedEndFiles.size() != 0 ||
            rdPairedEndFiles.size() != 0);

    if (readsProvided && ovlFile != "") {
        cerr << "[err] You must either provide reads file(s) or an ovl file\n";
        exit(EXIT_FAILURE);
    }

    if (!readsProvided && ovlFile == "") {
        cerr << "[err] You must either provide reads file(s) or an ovl file\n";
        exit(EXIT_FAILURE);
    }

    if (readsProvided && interactiveShell) {
        cerr << "[err] You most provide an ovl file to use the interactive shell\n";
        exit(EXIT_FAILURE);
    }

    if ((drPairedEndFiles.size() + rdPairedEndFiles.size()) % 2 != 0) {
        cerr << "[err] Number of paired reads files must be odd\n";
        exit(EXIT_FAILURE);
    }

    //to prevent copy/paste errors in the command line
    for (size_t i = 0; i < drPairedEndFiles.size(); i += 2) {
        if (drPairedEndFiles[i] == drPairedEndFiles[i + 1]) {
            cerr << "[err] Paired reads must be specified using two different files\n";
            cerr << "[err] (you specified the same file twice)\n";
            exit(EXIT_FAILURE);
        }
    }
    for (size_t i = 0; i < rdPairedEndFiles.size(); i += 2) {
        if (rdPairedEndFiles[i] == rdPairedEndFiles[i + 1]) {
            cerr << "[err] Mate pairs must be specified using two different files\n";
            cerr << "[err] you specified the same file twice\n";
            exit(EXIT_FAILURE);
        }
    }

    if (shortPeHorizon > longPeHorizon) {
        cerr << "[err] check PEhorizon setting\n";
        exit(EXIT_FAILURE);
    }

    if (truncateTo > 0 && truncateTo <= 30) {
        LOG.oss << "[WARNING] check -truncate parameter (" << truncateTo << ")\n";
        LOG.oss << "[WARNING] This parameter specifies the desired read length";
        LOG.flushStream(TOSTDOUT);
    }
}
