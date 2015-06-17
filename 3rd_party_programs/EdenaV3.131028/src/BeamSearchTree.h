/* 
 * File:   BeamSearchTree.h
 * Author: david
 *
 * Created on October 3, 2011, 11:17 AM
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

#ifndef BEAMSEARCHTREE_H
#define	BEAMSEARCHTREE_H

#include "NodeIt.h"
#include "PEMatches.h"
#include "Pairing.h"


#include<vector>
#include<deque>
#include<iostream>
using namespace std;

class BSNode;

class BeamSearchTree
{
     friend class PathFinder;
     friend class BackwardsWalker;
public:
    BeamSearchTree(Pairing *p);
    BeamSearchTree(const BeamSearchTree& orig);
    virtual ~BeamSearchTree();

    void init(NodeIt n);
    void init(vector<unsigned int> &path);
    
    //setting maxTreeDepth to 0 disable the depth limit
    unsigned int stepUpdateTree(unsigned int maxTreeDepth);
    
    void dot(string fileName, string title);
    void dot(string fileName, string title, string format);
    inline void setMaxBeamWidth(unsigned int n){maxBeamWidth=n;};
    
    inline unsigned int getUsableDistance(){return usableDistance;}
    inline size_t getQueueSize() const {return nodeQueue.size();}
    inline unsigned int getQueueAt(size_t i) const {return nodeQueue[i];}
    void getPath(unsigned int leaf, vector<unsigned int> &pp);

    int cleanTree(unsigned int chosenNode, unsigned int currentHead);
    void cleanTree(unsigned int currentHead);

   
    void computeExpectedNHits(unsigned int lib,
                            unsigned int bsNode,
                            vector<double>&summedProbs);
    
    void computeExpectedMeans(unsigned int lib,
                                        unsigned int bsNode,
                                        vector<double> &eMean);
    void computeExpectedMeans(unsigned int lib,
                                        vector<double> &eMean);

    inline void initPEMatchesStorage(){PEMS->initStorage();}

   
    
private:
    
    void dot(ostream &out,size_t node);
    
    void updateQueue(unsigned int node);
    void updateQueue_(unsigned int node);
            
public: //because of code completion issues in Netbeans...
    vector<BSNode> V;
    PEMatches * summedPEM;
    PEMatchesStorage *PEMS;
    
    unsigned int headBSNode; //required for the dot visualization
     //used for dot
    unsigned int chosen;
    
private:
    deque<size_t> nodeQueue, nodeQueueTmp;

    //safeguard
    unsigned int maxBeamWidth;
    
    //distance achieved during last tree update
    unsigned int usableDistance;
    
    //for information only
    static unsigned int maxAchievedBeamWidth;

    vector<size_t> nodeToMap;
    Pairing *P;   
};

class BSNode {
    friend class BeamSearchTree;
    friend class PathFinder;
    friend class BackwardsWalker;
    
public:
    BSNode();
    void init();
    void init(NodeIt &n);
    
    inline void getNodeIt(NodeIt &n){n.initNodeIt(nodeId,nodeDir);}
   // inline unsigned int getEdgeIndex() const {return edgeIndex;}
   // inline unsigned int getArrivalEdgesSize() const {return 
    inline unsigned int getNodeId() const {return nodeId;}
    inline bool getNodeDir() const {return nodeDir;}
    inline unsigned int getNodeLength() const {return nodeLength;} 
    inline unsigned int getDistance() const {return distance;}
    inline bool getFlag() const {return flag;}
    inline size_t getPLayout() const {return pLayout;}
    inline void setPLayout(size_t v) {pLayout=v;}
    
    inline unsigned int getSon(){return son;}
    inline unsigned int getBrother(){return brother;}
    inline unsigned int getFather(){return father;}
    inline unsigned int getOGNodeId(){return nodeId;}

    inline bool isActive(){return active;}
    inline void setActive(){active=true;}
    inline void setInactive(){active=false;}
    inline void setBranchingNodeIndex(size_t i){branchingNodeIndex=i;}
    inline size_t getBranchingNodeIndex(){return branchingNodeIndex;}
    
    void getNodeIt(NodeIt&) const;
    
    inline void initPEMatches(unsigned int n){ppeMatches = PEMS->getPPEMatches(n);}
    inline void setPEMatchesNull(){ppeMatches=0x0;}
    inline bool PEMatchesAllocated(){return ppeMatches!=0x0;}
    void addDistanceSample(unsigned int lib, unsigned int branching, unsigned int d);
    
    static vector<unsigned int> branchingNodeIds;
    static vector<unsigned int> nReadsChecked;
private:
    
    //tree connections
    unsigned int son;
    unsigned int brother;
    unsigned int father;
    
    //node values
    unsigned int edgeIndex;// edge index
    short ovSize;
    unsigned int nodeId; //nodeId in the graph
    bool nodeDir; //node direction r
    unsigned int nodeLength;
    unsigned int distance; //sequence length from the root
    
    PEMatches* ppeMatches;
    static PEMatchesStorage *PEMS;

    size_t branchingNodeIndex;
    bool active;
    
    size_t pLayout; //last read checked for backwards walk
    
    bool flag;//to see
    static vector<unsigned int> markedNode;
};

#endif	/* BEAMSEARCHTREE_H */

