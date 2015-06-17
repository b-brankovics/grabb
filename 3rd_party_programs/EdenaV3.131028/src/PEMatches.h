/* 
 * File:   PEMatches.h
 * Author: david
 *
 * Created on April 8, 2012, 11:34 PM
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

#ifndef PEMATCHES_H
#define	PEMATCHES_H

#include "Pairing.h"
#include <vector>
#include <cmath>

using namespace std;

class PEMatches;

class PEMatchesStorage {
    
public:
    
    PEMatchesStorage(Pairing *p);
    ~PEMatchesStorage();
       
    void allocate();
    void doubleAllocatedSize();
    void freeMemory();
    
    void initStorage();
    
    //return an initialized element
    PEMatches* getPPEMatches(unsigned int nChoice);
    
private:
    
    PEMatches  **V;
    Pairing *P;
    size_t nElement;
    size_t maxSize;
};

class PEMatches {
public:
    PEMatches();
    PEMatches(const PEMatches& orig);
    virtual ~PEMatches();
    
    void init(unsigned int n);
    void allocate();
    void doubleAllocationSize();
    void freeMemory();
    inline unsigned int getCount(unsigned int lib, unsigned i){return counts[lib-1][i];}
    inline unsigned int getNChoice(){return nBranch;}
    inline double getMean(unsigned int lib, unsigned int i){return sumDist[lib-1][i]/counts[lib-1][i];}
    double getVar(unsigned int lib,unsigned int i);
    inline double getSD(unsigned int lib, unsigned int i){return sqrt( getVar(lib,i));}
    inline void addMatch(unsigned int lib, unsigned int i, double distance){counts[lib-1][i]++;sumDist[lib-1][i]+=distance;sumDistSquare[lib-1][i]+=(distance*distance);}
   
    void getMaxCount(unsigned int lib, unsigned int &sumCount, unsigned int &max, size_t&maxIndex);
    inline int getSumCount(unsigned int lib){unsigned int sum=0; for(unsigned int i=0; i<nBranch; i++)sum+=counts[lib-1][i];return sum;}

    static Pairing *P;
    static unsigned int nLibrary;
    
private:
    
    
    unsigned int **counts;
    double **sumDist;
    double **sumDistSquare;
    unsigned int maxBranch;
    unsigned int nBranch;
};

#endif	/* PEMATCHES_H */

