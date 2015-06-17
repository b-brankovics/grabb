/* 
 * File:   PathFinder.h
 * Author: david
 *
 * Created on March 4, 2011, 2:17 PM
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


#ifndef BACKWARDWALKER_H
#define	BACKWARDWALKER_H
#include "NodeIt.h"
#include "BeamSearchTree.h"

#include <vector>
#include <set>
#include <deque>
using namespace std;


struct BsNodeIdentifier
{
    int nIdDir; //This is the index key
                //both nodeID and nodeDir are coded in the same variable
                 
    unsigned int BSNode; //point to a node in the BeamSearchTree
  //  mutable int pathDistance; 

    inline bool operator< (const BsNodeIdentifier& b) const
    {      
        return this->nIdDir < b.nIdDir;
    }

    inline unsigned int getNodeId() const {return nIdDir > 0 ? nIdDir : -nIdDir;} 
    inline bool getDir() const {return nIdDir > 0 ? true : false;}
    inline void setNodeId(unsigned int n, bool d){nIdDir= d ? (int)n : (int)-n;}
    inline unsigned int getBSNode() const {return BSNode;}

};

class BackwardsWalker
{
public:
    BackwardsWalker();
    BackwardsWalker(const BackwardsWalker& orig);
    virtual ~BackwardsWalker();

    void init(OverlapsGraph *);
    inline bool isAllocated(){return G!=0x0;}
    void clean();
    unsigned int initWalker(unsigned source,
                            bool dir,
                            unsigned int minNPair,
                            float minRatio,
                            unsigned int maxRedundancy,
                            unsigned int maxBeamWidth);
    
    unsigned int initWalker(vector <unsigned int> &ePath,
                            unsigned int minNPair,
                            float minRatio,
                            bool checkD,
                            unsigned int maxBeamWidth);
    inline bool isInitialized(){return headBSNode >=1;}
    
    void initHistoryMap();
    int acceptChoice(unsigned int chosenNode);
    int elongate(vector<unsigned int> &path);
    
    //return value = update status
    //BSNode get the value of the chosen node is the "decision zone"
    //BSNode=0 if no choice has be made
    int stepElongate(unsigned int maxD, unsigned int &BSNode);
    int updateHistoryMap(unsigned int BSNode);
  
    void backwardsMap();
    unsigned int getDecision();
    void CountPEMatches(vector<unsigned int> path, size_t p1, size_t p2);
   
    inline unsigned int getMaxAchievedBeamWidth(){return BST->maxAchievedBeamWidth;}   
    inline unsigned int getUsableDistance(){return BST->getUsableDistance();}
    inline void setMinNPair(unsigned int n) {minNPair=n;}
    inline void setMinRatio(float v){minRatio=v;}
    inline void setMaxRedundancy(unsigned int v){maxRedundancy=v;}
    inline void setCheckDistances(bool f){checkDistances = f;}
    inline bool getCheckDistances(){return checkDistances;}
   
    inline void dot(string fileName, string format){BST->dot(fileName,"BKWDTree",format);}
    inline void dot(string fileName){BST->dot(fileName,"BKWDTree");}
   
    void getPath(unsigned int treeNode, vector<unsigned int> &p);
    void getPath(vector<unsigned int> &p); //path from headBSNode
    static unsigned int st_count;
    
private:

    unsigned int headBSNode;
    unsigned int rL;
    unsigned int minNPair;
    float minRatio;
    unsigned int maxRedundancy;
    
    bool checkDistances;
    
    unsigned int nPotentialElongation;//number of sons under headBSNode
    
    multiset<BsNodeIdentifier> historyMap;
    vector<size_t> activeBsNodes; //will contain the BSNode to be examined
   
    BeamSearchTree *BST;
    OverlapsGraph *G;
    Node *N;
    Pairing *P;
    ReadsLayout *L;
};

#endif	/* BACKWARDWALKER_H */

