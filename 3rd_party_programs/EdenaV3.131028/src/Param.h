/* 
 * File:   Param.h
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

#ifndef PARAM_H
#define	PARAM_H

#include <string>
#include <vector>
using namespace std;

class Param {
public:
    Param();
    void init(); //setDefaultValues;
    void parseCommandLine(int argc, char**argv);
    
    Param(const Param& orig);
    virtual ~Param();
    
    inline bool overlappingMode(){return readsProvided;}
    inline bool assemblyMode(){return !readsProvided;}

    bool readsProvided;
    
    //user parameters
    int minOverlap; //computed during overlapping step
    int overlapCutoff; //considered during assembly
    int truncateTo;
    unsigned int maxDeadEndLength;
    double minNodeCoverage;
    double minContigCoverage;
    double targetSize;
    int nThreads;
    int minContigSize;
    int minNPair;
    double minRatio;
    int maxContigRedundancy;
    unsigned int shortPeHorizon;
    unsigned int longPeHorizon;
    int trimContigEnds; //coverage cutoff for contig ends;
    
    
    string prefix;
    vector<string> singleEndFiles;
    vector<string> drPairedEndFiles;
    vector<string> rdPairedEndFiles;
    string ovlFile;
    
    bool contextualCleaning;
    bool cleanGraph;
    bool interactiveShell;
    bool discardNonUsable;
    bool discardSuspicousNodes;
    bool cleanBubbles;
    bool trimRed;
    
    //non documented dev flags
    bool writeOSG;
    bool checkGraph;
    
private:
    void checkParamConsistency();
    
};

#endif	/* PARAM_H */

