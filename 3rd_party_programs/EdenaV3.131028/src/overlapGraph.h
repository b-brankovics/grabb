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

#ifndef OVERLAPS_GRAPH_H
#define OVERLAPS_GRAPH_H
#include "node.h"
#include "readsStorage.h"
#include "BackwardsWalker.h"
#include "globalFunc.h"
#include <string>
#include <vector>
#include <iostream>
#include <ctime>
#include <cstring>
#include <iomanip>

using namespace std;


class PathTree;
class ReadsLayout;
class Pairing;
class ReadsStorage;

extern ofstream outGlob;
class AssemblyProgress;

class OverlapsGraph
{
    public:
        OverlapsGraph();
        ~OverlapsGraph();

    ReadsStorage *R;  //reads
    Pairing *P; //pairing information
    ReadsLayout *L; //Layouts

    unsigned int maxNodes; //maximum number of nodes that can be stored
    unsigned int nNodes; //actual number of nodes
    unsigned int nEdges;
    int minOverlap;
    double targetSizeEstimation;
    string ovlVersion;
    static AssemblyProgress AP;
    static deque<unsigned int> nodeQueue;
    static unsigned int g_count;

    void init(ReadsStorage*, Pairing*, ReadsLayout*);
    void allocateNodeTab(unsigned int maxN);

    bool loadData(string);
    bool saveData(string);

    
    void computeOverlaps(unsigned int nThreads);

    //DEV
    void computeIrreductibeEdges();

    //START quick multithreading implementation
    struct Thread_param
    {
        unsigned int start;
        unsigned int end;
        OverlapsGraph* _this;
        void (OverlapsGraph::*funcPTR)(void*);
    };
   
    void overlapNodeRange(void* dat_s);
    static void* thread_maker(void* dat_s);
    //END  quick multithreading implementation
    
    
    void condense(bool verbose, bool sortEdges);
    void renumber(bool verbose, bool sortEdges);

    void assemble(string fileName, unsigned int minSize , double minCoverage, int trim, int minNPair, double minRatio, unsigned int maxRedundancy );
   
    int PEelongate(vector<unsigned int>& path,
        unsigned int node,
        bool direction,
        int minNPair,
        double minRatio,
        unsigned int maxRedundancy);
    
    int getPathBetweenTwoReads(unsigned int id1, bool d1, unsigned int id2, bool d2, int minLength, int maxLength, double maxLogNLeaves);

    unsigned int getNOverlap(unsigned int id);
    unsigned int getNRightOverlap(unsigned int id);
    unsigned int getNLeftOverlap(unsigned int id);

    unsigned int getReadLength();

    inline unsigned int getNNodes() const {return nNodes;}
    inline unsigned int getNEdges() const {return nEdges;}
    inline  int getMinOverlap() const {return minOverlap;}

    void printSummary(unsigned int id, bool direct, bool right);
    void nicePrint(unsigned int id, bool direct);

    void overlapSizeCutoff(unsigned int cutoff);
    
    //Remove short overlaps according to the context
    void computeEdgesProb(double cutOffProb);//min prob to be considered as reliable
    unsigned int removeEdgesByValue(float suspectCutoff, float GIncoherentCutoff);

    unsigned int coverageCutoff(unsigned int cutoff);
    void removeTransitiveEdges();
    unsigned int countEdges();
    bool checkConsistency();

    unsigned int discardShortOrphanNodes(unsigned int &nDiscardedReads);
    unsigned int discardSuspiciousNode(double minCov);
    unsigned int identifyDeadEnd(unsigned int &dLimit, unsigned int &nrreads);
    bool identifyDeadEnd(unsigned int nodeId,
                                 bool dir,
                                 unsigned int distance,
                                 unsigned int &nDeadEnd,
                                 unsigned int dLimit);
    
    unsigned int bubbles(double minCov);
    
    void estimatePairedDistance(unsigned int PEhorizon, unsigned int MPhorizon, double nsd, string prefix);
    void testChimera(unsigned int maxLength);
    
    
    
    // an "ePath" is a path defined by a list of edges. Edges are given by their absolute index
    // this is the unambiguous, machine readable, format.
    
    // a "nodePath" is a list of oriented nodes, provided by a list of node and a list of orientation.
    // this format is used to represent a path in a way easily readable for a human being.
    
    string ePathToSeq(const vector<unsigned int> &path);
    
    unsigned int getNReadsInEPath(const vector<unsigned int> &path);
    
    void ePathToNodePath(const vector<unsigned int> &path, vector<unsigned int> &nodePath, vector<bool> &vDir);
    void nodePathToEPath(const vector<unsigned int> &nodePath, const vector<bool> &vDir, vector<unsigned int> &path);
    void reverseEPath(vector<unsigned int> & path);
    void ePathToCov(const vector<unsigned int> &path, vector<unsigned int> &cov);
    void ePathTabularInfo(const vector<unsigned int> &path, ostream &out);
    
    
    void nodeListToEPath(string args, vector<unsigned int> &ePath);//for DEV mode
    void estimateCoverage(double &minCoverage, double &targetSize);

    static Node * nodesTab;//start at 1
    
    struct nodeRank
    {
	unsigned int nodeIndex;
	inline int operator< (const nodeRank &rhs) const
	{
		if( nodesTab[this->nodeIndex].getSequenceLength() >
                        nodesTab[rhs.nodeIndex].getSequenceLength())
			return 1;
		else
			return 0;
        }
    };
    vector<nodeRank> rank;
    void sortNodesByLength();

    struct orderCov
    {
        double cov;
        unsigned int weight;

        inline int operator<(const orderCov & rhs) const
        {
            if (cov < rhs.cov)
                return 1;
            else
                return 0;
        }
    };
    
    private:

        void freeMemory();
        unsigned int depthLimit;
};


class AssemblyProgress
{
public:
    AssemblyProgress(){ init();}
    ~AssemblyProgress(){};
    
    inline void addKb(unsigned int kb) {totKb+=kb;}
    inline void setKb(unsigned int kb) {totKb=kb;}
    inline void incrementNContigs(){nContigs++;}
    inline void setNContigs(unsigned int nc) {nContigs=nc;
    }

    void printProgress(ostream &out) {printProgress(out,false);}
    
    void printProgress(ostream &out, bool update) {
        
        currentClock = clock();
        if (currentClock - lastClock > CLOCKS_PER_SEC / 10) {
            pulseIndex++;
            lastClock = currentClock;
            update=true;
        }
       
        if (update)
        {
            out << "\r"
                    << setfill(' ') << setw(8) << nContigs << " (nContigs)"
                    << setfill(' ') << setw(11) << smartDNALength(totKb) << " (totSize)"
                    //<< "   " << pulse[pulseIndex % 4] << " (pulse)"
                    << flush;
         }
    }
    void init()
    {
        nContigs=0;
        totKb=0;
        lastClock=0;
        pulseIndex=0;
        strcpy(pulse,".oOo");
    }
    
private:
    
    unsigned int nContigs;
    unsigned int totKb;
    char pulse[5];
    int pulseIndex;
    clock_t lastClock;
    clock_t currentClock;
};

#endif // OVERLPAPS_GRAPH_H
