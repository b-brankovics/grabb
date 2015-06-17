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

#ifndef NODE_H
#define NODE_H
#include <cstddef>
#include <vector>
#include <deque>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "readsLayout.h"
//#include "VPaths.h"


//a parameter that limit the depth during the graph path searches
//#define MAX_LOG_N_LEAVES 16.0
#define MAX_LOG_N_LEAVES 25.0
using namespace std;

class ReadsStorage;
class Routing;
//class Layouts;
class NodeIt;
class OverlapsGraph;
class Pairing;

struct NodeMem {
    unsigned int id;
    unsigned int incomIndex;
    int overHang;
    bool dir;
    bool isDelegated;
};


//quick addition to be used by the new "sortById" method
//the whole graph structure should be recoded one day...
struct Edge {
    unsigned int ovId;
    short ovSize;
    float edgeValue;
    unsigned char edgeBitFlag;

    inline int operator<(const Edge &e) const {
        if (ovId < e.ovId)
            return 1;
        else
            return 0;
    }
};

class Node {
    friend class OverlapsGraph;
    friend class NodeIt;

private:

    float *edgeValue;
    unsigned char *edgeBitsFlag;
    //signed char *ovSize;
    short int *ovSize;
    unsigned int *ovId;

    unsigned int nOv; //nEdges
    unsigned int nOvRight;

    unsigned short int value;
    unsigned int layout;
    char flagsBit;

    static Node *N;
    static ReadsStorage *R;
    static ReadsLayout *L;
    static OverlapsGraph *G;
    static Pairing *P;
    static unsigned int nNodes;
    static unsigned int S_counter;
    static bool edgeSorted;

    static vector<unsigned int> path; //used to store a path: PATH[0]=starting node, then indexes of the edges forming the path;

    static unsigned int nPathFound;
    static double pathLengthMean, pathLengthSD; //used to store the mean and SD of the length distribution
    static vector<unsigned int> nodeList;
    static vector<unsigned int> dotNodeList;

    static deque<NodeMem> stack;

public:

    Node();
    ~Node();

    static unsigned int maxLog;

    void allocate(int n); //allocate memory for n edges
    void allocateEdgeValue(); //only used in assembly mode
    void allocateEdgeFlag();
    void reallocMemory();
    void freeMemory();

    inline bool isAllocated() {
        return ovId != 0x0;
    }
    static void setEdgeSortedFlag(bool v);

    void computeOverlaps();

    int searchPathToPaired(
            bool currentDir,
            unsigned int targetNode,
            bool targetDir,
            int minDistance,
            int maxDistance,
            int shift,
            double maxLogNLeaves);

    int searchPathToPaired(
            bool currentDir,
            unsigned int targetNode,
            bool targetDir,
            int minDistance,
            int maxDistance,
            int shift,
            int currentDistance,
            double &cancel,
            double maxLogNLeaves);

    //for graph condensing
    unsigned int getLongestNonAmbiguousPath(string &s,
            vector<unsigned int> &coverage,
            unsigned int &leftId,
            bool &leftDir,
            unsigned int &rightId,
            bool &rightDir);

    void sortByIds(); //assume that edges are already sorted by size
    void sortByIds(unsigned int, unsigned int);

    //DOT
    static void initDot(ostream &out);

    inline void dotAround(ostream &out, int depthLimit) {
        dotAround(out, depthLimit, 0);
    }
    void dotAround(ostream &out, int depthLimit, int currentDepth);
    static void closeDot(ostream &out);
    void dotLocalGraph(int depth, string fileName);
    void writeDotGraph(ostream &out, int &depthLimit, int currentDepth);

    //accessors

    inline unsigned int getThisId() const {
        return (this-N);
    }

    inline unsigned int getNRightOverlap() const {
        return nOvRight;
    };

    inline unsigned int getNLeftOverlap() const {
        return nOv - nOvRight;
    };

    inline unsigned int getNOverlap() const {
        return nOv;
    };

    inline unsigned int getNEdges(bool dir) const {
        if (dir)
            return nOvRight;
        else
            return nOv - nOvRight;
    };

    unsigned int getSequenceLength();

    double getCoverage();
    double getUnbiaisedCoverage(double& nodeLength);
    void getVCoverage(vector<unsigned int> &cov, bool direction);
    unsigned int getNReads();

    string getDirectSequence();
    string getReverseSequence();

    inline string getSequence(bool dir) {
        if (dir)return getDirectSequence();
        else return getReverseSequence();
    }

    inline unsigned int getLayout() {
        return layout;
    }

    inline void setLayout(size_t l) {
        layout = l;
    }
    unsigned int markReads(bool dir);

    // "dir" corresponds to the neighbor orientation, and not to the overlap orientation

    inline void getNeighbor(bool direction, size_t index, unsigned int &id, unsigned int &size, bool &dir) {
        if (!direction)
            index += nOvRight;

        id = ovId[index];
        int s = (int) ovSize[index];
        if (s > 0) {
            size = s;
            dir = direction;
        } else {
            size = -s;
            dir = !direction;
        }
    }
    void getNeighbor(bool direction, size_t index, NodeIt&);

    unsigned int getRightOverlap(unsigned int index, unsigned int &id, unsigned int&size, bool &direction);
    unsigned int getLeftOverlap(unsigned int index, unsigned int &id, unsigned int&size, bool &direction);
    unsigned int getOverlap(unsigned int index, unsigned int &id, unsigned int&size, bool &direction);
    unsigned int getRightOverlapSize(unsigned int index);
    unsigned int getLeftOverlapSize(unsigned int index);
    unsigned int getOverlapSize(unsigned int index);
    unsigned int getRightOverlapId(unsigned int index);
    unsigned int getLeftOverlapId(unsigned int index);
    unsigned int getOverlapId(unsigned int index);
    bool getOverlapDirection(unsigned int index);

    unsigned int getNext(bool dir, unsigned int index, unsigned int& id, unsigned int& size, bool&direction);

    //return edge index if it exists or nOv if not;
    //uses the size to speed up the search in the edge list
    unsigned int getEdgeIndex(bool right, unsigned int id, unsigned int size, bool direction);

    void overlapSizeCutoff(unsigned int cutoff); //do not preserve consistance
    void removeShortBranchingOverlaps(unsigned int cutoff); //preserve consistance

    void markTransitiveEdges(); //mark transitive edges as discarded
    void anotherMarkTransitiveEdges();
    bool testReciprocal(); //test
    bool testReciprocal2(); //test

    //must be applied to ALL node
    int removeMarkedEdge();
    int removeEdgesByValue(float suspectCutoff, float GIncoherentCutoff);

    void computeEdgesProb(unsigned int maxD, unsigned int max, double reliableCutoff);

    unsigned int bubble(bool dir, double minCov);

    //the one to be called
    unsigned int sampleOverHangs(
            bool direction,
            unsigned int maxD,
            unsigned int maxN,
            unsigned int *distr1);
    //do not call this one!
    unsigned int sampleOverHangsRec(
            bool direction,
            unsigned int maxD,
            unsigned int maxN,
            unsigned int currentD,
            unsigned int currentN,
            unsigned int *distr1);

    unsigned int getReciprocal(unsigned int index);

    void removeEdge(unsigned int index); //remove edge i and its reciprocal
    void removeEdgeNR(unsigned int index); //do not care of the reciprocal
    void removeEdges(); //remove all edges of a node and reciprocal

    inline void isolate() {
        removeEdges();
    };

    void flagReads(char state);
    void removeAllEdges();

    void setEdgeFlag(unsigned int index, bool state); //mark edge i and its reciprocal
    void setEdgeFlagUnrec(unsigned int index, bool state); //do not take care of the reciprocal
    void initEdgeFlags(bool state);

    inline bool getEdgeFlag(unsigned int index) {
        return edgeBitsFlag[index];
    }


    void initializeEdgeValues(float v); //set all edge values to zero
    void setEdgeValue(unsigned int index, float value);
    float getEdgeValue(size_t i);
    void incrementEdgeValue(size_t i);
    void decrementEdgeValue(size_t i);
    void incrementRightEdgeValue(size_t i);
    void decrementRightEdgeValue(size_t i);
    void incrementLeftEdgeValue(size_t i);
    void decrementLeftEdgeValue(size_t i);

    bool isBranching(); //true if branching on the right or left hand side
    bool isBranching(bool dir);
    bool isInternal(); //true if not branching
    bool hasSingleEdge(bool dir);
    bool isEnding();
    bool isEnding(bool dir);
    bool hasNoEdge();
    bool hasNeighbors(); //at least connected to another node (not self connected)
    bool hasSingleSideOV(); //Ov on a single side
    bool hasBothSidesOV();

    //bit flags:
    //BEWARE OF POSSIBLE FLAGS COLLISIONS

    //used during graph condensing

    inline bool isExtended() {
        return (flagsBit & 1);
    }

    inline void setExtended() {
        flagsBit |= 1;
    }

    inline void unsetExtended() {
        flagsBit &= ~1;
    }

    //graph dot

    inline bool isVisited() {
        return flagsBit & 6;
    }//Either L or R

    inline void setVisited() {
        flagsBit |= 6;
    }//Either L or R

    inline void unsetVisited() {
        flagsBit &= ~6;
    }//Either L or R

    //flag nodes in history (backward mapping)

    inline bool isVisited(bool dir) {
        if (dir) return flagsBit & 2;
        else return flagsBit & 4;
    }

    inline void setVisited(bool dir) {
        dir ? flagsBit |= 2 : flagsBit |= 4;
    }

    inline void unsetVisited(bool dir) {
        dir ? flagsBit &= ~2 : flagsBit &= ~4;
    }

    //for graph dot

    inline bool isTargeted() {
        return flagsBit & 8;
    }

    inline void setTargeted() {
        flagsBit |= 8;
    }

    inline void unsetTargeted() {
        flagsBit &= ~8;
    }

    //transitive edges reduction

    inline bool isInplay(bool dir) {
        if (dir) return flagsBit & 2;
        else return flagsBit & 4;
    }

    inline void setInplay(bool dir) {
        dir ? flagsBit |= 2 : flagsBit |= 4;
    }

    inline void unsetInplay(bool dir) {
        dir ? flagsBit &= ~2 : flagsBit &= ~4;
    }

    //heuristic for redundancy detection during layout phase

    inline void setAlreadyUsed(bool dir) {
        if (dir) flagsBit |= 8;
        else flagsBit |= 16;
    }

    inline bool isAlreadyUsed(bool dir) {
        if (dir) return flagsBit & 8;
        else return flagsBit & 16;
    }

    inline void setAlreadyUsed() {
        flagsBit |= 24;
    }

    inline bool isAlreadyUsed() {
        return flagsBit & 24;
    }

    inline void unsetAlreadyUsed() {
        flagsBit &= !24;
    }

    //used during transitive edge reduction

    inline void setMultipleEdges(bool dir) {
        if (dir) flagsBit |= 8;
        else flagsBit |= 16;
    }

    inline bool hasMultipleEdges(bool dir) {
        if (dir) return flagsBit & 8;
        else return flagsBit & 16;
    }

    inline void unsetMultiplesEdges() {
        flagsBit &= !24;
    }

    inline bool isEliminated(bool dir) {
        if (dir) return flagsBit & 32;
        else return flagsBit & 64;
    }

    inline void setEliminated(bool dir) {
        if (dir) flagsBit |= 32;
        else flagsBit |= 64;
    }

    inline void unsetEliminated(bool dir) {
        if (dir) flagsBit &= ~32;
        else flagsBit &= ~64;
    }

    inline bool isEliminated() {
        return (flagsBit & 96);
    }

    inline void setEliminated() {
        flagsBit |= 96;
    }

    inline void unsetEliminated() {
        flagsBit &= ~96;
    }

    inline bool isDiscarded() {
        return (flagsBit & 128);
    }

    inline void setDiscarded() {
        flagsBit |= 128;
        flagReads(1);
    }

    inline void unsetDiscarded() {
        flagsBit &= ~128;
        flagReads(0);
    }

    inline bool ovComputed() {
        return (flagsBit & 128);
    }

    inline void setOvComputed() {
        flagsBit |= 128;
    }

    inline void unsetOvComputed() {
        flagsBit &= ~128;
    }

    inline bool isReduced() {
        return flagsBit & 1;
    }

    inline void setReduced() {
        flagsBit |= 1;
    }

    inline void unsetReduced() {
        flagsBit &= ~1;
    }

    //edge flags

    inline bool getEdgeFlagA(size_t i) {
        return (edgeBitsFlag[i] & 1);
    }

    inline void setEdgeFlagA(size_t i) {
        edgeBitsFlag[i] |= 1;
    }

    inline void unsetEdgeFlagA(size_t i) {
        edgeBitsFlag[i] &= ~1;
    }

    inline bool getEdgeFlagB(size_t i) {
        return (edgeBitsFlag[i] & 2);
    }

    inline void setEdgeFlagB(size_t i) {
        edgeBitsFlag[i] |= 2;
    }

    inline void unsetEdgeFlagB(size_t i) {
        edgeBitsFlag[i] &= ~2;
    }

    inline bool getEdgeFlagC(size_t i) {
        return (edgeBitsFlag[i] & 4);
    }

    inline void setEdgeFlagC(size_t i) {
        edgeBitsFlag[i] |= 4;
    }

    inline void unsetEdgeFlagC(size_t i) {
        edgeBitsFlag[i] &= ~4;
    }

    inline bool getEdgeFlagD(size_t i) {
        return (edgeBitsFlag[i] & 8);
    }

    inline void setEdgeFlagD(size_t i) {
        edgeBitsFlag[i] |= 8;
    }

    inline void unsetEdgeFlagD(size_t i) {
        edgeBitsFlag[i] &= ~8;
    }

    inline bool getEdgeFlagE(size_t i) {
        return (edgeBitsFlag[i] & 16);
    }

    inline void setEdgeFlagE(size_t i) {
        edgeBitsFlag[i] |= 16;
    }

    inline void unsetEdgeFlagE(size_t i) {
        edgeBitsFlag[i] &= ~16;
    }

    inline bool getEdgeFlagF(size_t i) {
        return (edgeBitsFlag[i] & 32);
    }

    inline void setEdgeFlagF(size_t i) {
        edgeBitsFlag[i] |= 32;
    }

    inline void unsetEdgeFlagF(size_t i) {
        edgeBitsFlag[i] &= ~32;
    }

    inline bool getEdgeFlagG(size_t i) {
        return (edgeBitsFlag[i] & 64);
    }

    inline void setEdgeFlagG(size_t i) {
        edgeBitsFlag[i] |= 64;
    }

    inline void unsetEdgeFlagG(size_t i) {
        edgeBitsFlag[i] &= ~64;
    }

    inline bool getEdgeFlagH(size_t i) {
        return (edgeBitsFlag[i] & 128);
    }

    inline void setEdgeFlagH(size_t i) {
        edgeBitsFlag[i] |= 128;
    }

    inline void unsetEdgeFlagH(size_t i) {
        edgeBitsFlag[i] &= ~128;
    }
};

#endif // NODE_H
