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

#ifndef PAIRING_H
#define PAIRING_H

#include "readsStorage.h"
#include "ActualDistribution.h"
#include <vector>
#include <iostream>

using namespace std;
class ReadsStorage;

class PELibrary
{
    friend class Pairing;
public:
    PELibrary();
    virtual ~PELibrary();
    
    void save(ostream&);
    bool load(istream&);
    
    inline void setMateOrientation(int mo){mateOrientation=mo;}
    inline void setStartIndex(size_t i){startIndex=i;}
    inline void setEndIndex(size_t i){endIndex=i;}
    inline void setMeanCloneSize(double v){meanCloneSize=v;}
    inline void setSDCloneSize(double v){sdCloneSize=v;}
    inline void setNEchant(unsigned int n){nEchant=n;}
    inline void setDistanceMin(unsigned int d){distanceMin=d;}
    inline void setDistanceMax(unsigned int d){distanceMax=d;}
    inline void setNUsableMates(unsigned int n){nUsableMates=n;}
    inline void setExpectedMateCoverage(double v){expectedMateCoverage=v;}
   
    inline int getMateOrientation(){return mateOrientation;}
    inline size_t getStartIndex(){return startIndex;}
    inline size_t getEndIndex(){return endIndex;}
    inline double getMeanCloneSize(){return meanCloneSize;}
    inline double getSdCloneSize(){return sdCloneSize;}
    inline unsigned int getNEchant(){return nEchant;}
    inline unsigned int getDistanceMin(){return distanceMin;}
    inline unsigned int getDistanceMax(){return distanceMax;}
    inline unsigned int getNUsableMates(){return nUsableMates;}
    inline double getExpectedMateCoverage(){return expectedMateCoverage;}
    
    inline unsigned int getSize(){return endIndex-startIndex+1;}
   
   
    inline int operator< (const PELibrary &rhs) const
	{
		if( meanCloneSize < rhs.meanCloneSize)
			return 1;
		else
			return 0;
        }
    static int readLength;
  
     ActualDistribution lengthDistr;
private:
    int mateOrientation; //1 = >< (direct-reverse)
                         //2 = <> (reverse-direct)
                         //0 = unpaired
    size_t startIndex;
    size_t endIndex;
    double meanCloneSize; //full sample
    double sdCloneSize;   //full sample
    unsigned int nEchant;  //fullSample
    
    unsigned int distanceMin; //sample inf bound
    unsigned int distanceMax; //sample sup bound
    
    double sampleSize; //bounded sample, not used yet
    unsigned int nUsableMates; // an estimation
    double expectedMateCoverage;
};

class Pairing
{
public:
    Pairing();
    virtual ~Pairing();

    void init();
    unsigned int getNPairing() const;
    unsigned int getR1(size_t index) const;
    unsigned int getR2(size_t index) const;
    unsigned int getPairing(unsigned int readId, bool &pair2) const;
    void getDistanceRange(unsigned int readId, unsigned int &min, unsigned int &max, int &mateOrientation);
    inline unsigned int getPairing(unsigned int readId) const {return fastPairing[readId];}
    
    inline unsigned int getNLibrary(){return vPeLibrary.size();}
    inline size_t getPeLibraryID(unsigned int read){return libraryId[read];} // 1-based
    inline int getMateOrientation(unsigned int read){return libraryId[read]<=nDRLib?1:2;}
    inline vector<PELibrary>::iterator getLibraryIt(size_t i){return vPeLibrary.begin()+i;}
    void updatePERange(double nsd);
    
    inline unsigned int getGlobalMaxAllowedDistance(){return globalMaxAllowedDistance;}
    inline void setGlobalMaxAllowedDistance(unsigned int v){globalMaxAllowedDistance=v;}
    inline unsigned int getPairedEndsMaxD(){return pairedEndsMaxD;}
    inline unsigned int getMatePairsMaxD(){return matePairsMaxD;}

    void addPair(unsigned int R1, unsigned int R2);
    void startNewLibrary(int mateOrientation);
    void endLibrary();

    void buildIndex(unsigned int);
    void updatePeLibraryIndex();
    void cleanMemory();

    void save(ostream&);
    bool load(istream&);
    
    //remove the pairing information (for test purpose only!)
    void unPair(size_t lib);
    
   unsigned int countNValidPair(size_t lib); //number of pairs for which both are non-discarded
   inline void setRpointer(ReadsStorage *RS){R=RS;}
   
private:

    vector<unsigned int> R1;//paired1
    vector<unsigned int> R2;//paired2
    size_t headp;
    ReadsStorage *R;
    unsigned int *pairingIndex;
    unsigned int *fastPairing;
    unsigned char *libraryId;
    unsigned int numberOfReads;
    vector<PELibrary> vPeLibrary;
    unsigned int nDRLib;
    unsigned int globalMaxAllowedDistance;
    unsigned int pairedEndsMaxD;
    unsigned int matePairsMaxD;
};


#endif // PAIRING_H
