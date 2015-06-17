/* 
 * File:   PathFinder.cpp
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

#include "BackwardsWalker.h"
#include "overlapGraph.h"
#include "node.h"
#include "Pairing.h"
#include "readsLayout.h"
//#include "VPaths.h"
#include "globalFunc.h"

#include <algorithm>
#include <limits>
#include <sstream>


extern bool DEV_INFO;
unsigned int BackwardsWalker::st_count = 0;

BackwardsWalker::BackwardsWalker() {
    BST = 0x0;
    G = 0x0;
    P = 0x0;
    N = 0x0;
    L = 0x0;
}

BackwardsWalker::BackwardsWalker(const BackwardsWalker& orig) {
}

BackwardsWalker::~BackwardsWalker() {
    clean();
}

void BackwardsWalker::clean() {
    if (BST != 0x0)
        delete BST;

    BST = 0x0;
    G = 0x0;
    P = 0x0;
    N = 0x0;
    L = 0x0;

}

void BackwardsWalker::init(OverlapsGraph *g) {
    clean();
    BST = new BeamSearchTree(g->P);
    headBSNode = 0;
    G = g;
    P = G->P;
    N = G->nodesTab;
    L = G->L;
    rL = G->getReadLength();
}

unsigned int BackwardsWalker::initWalker(unsigned source,
        bool dir,
        unsigned int minNPair,
        float minRatio,
        unsigned int maxRedundancy,
        unsigned int maxBW) {
    NodeIt n;

    n.initNodeIt(source, dir);
    BST->init(n);
    headBSNode = 1;
    setMinNPair(minNPair);
    setMinRatio(minRatio);
    setMaxRedundancy(maxRedundancy);
    BST->setMaxBeamWidth(maxBW);

    BsNodeIdentifier bsNId;
    bsNId.setNodeId(source, dir);
    bsNId.BSNode = 1;

    initHistoryMap();
    activeBsNodes.clear();

    historyMap.insert(bsNId);
    G->nodesTab[source].setVisited(dir);
    G->nodesTab[source].setAlreadyUsed();
    checkDistances = true; //by default

    return 0;
}

unsigned int BackwardsWalker::initWalker(vector <unsigned int> &ePath,
        unsigned int minNPair,
        float minRatio,
        bool checkD,
        unsigned int maxBeamWidth) {
    BST->init(ePath);
    headBSNode = BST->headBSNode;
    setMinNPair(minNPair);
    setMinRatio(minRatio);
    BST->setMaxBeamWidth(maxBeamWidth);
    activeBsNodes.clear();
    initHistoryMap();
    checkDistances = checkD;

    BsNodeIdentifier bsNId;
    for (unsigned int i = 1; i < BST->V.size(); i++) {
        bsNId.setNodeId(BST->V[i].getNodeId(), BST->V[i].getNodeDir());
        bsNId.BSNode = i;
        historyMap.insert(bsNId);
        G->nodesTab[BST->V[i].getNodeId()].setVisited(BST->V[i].getNodeDir());

    }
    return 0;
}

void BackwardsWalker::initHistoryMap() {
    int OGnode;
    multiset<BsNodeIdentifier>::iterator it;
    for (multiset<BsNodeIdentifier>::iterator it = historyMap.begin(); it != historyMap.end(); it++) {
        OGnode = it->nIdDir;
        if (OGnode < 0)
            OGnode = -OGnode;
        G->nodesTab[OGnode].unsetVisited();
    }
    historyMap.clear();

}

int BackwardsWalker::acceptChoice(unsigned int chosenBsNode) {
    if (chosenBsNode == 0) {
        return 0;
    }

    bool ok = false;

    unsigned int node = BST->V[headBSNode].getSon();

    while (node != 0) {
        if (chosenBsNode == node) {
            ok = true;
            break;
        }
        node = BST->V[node].getBrother();
    }

    if (ok == false) {
        return 0;
    }
    for (unsigned int i = 0; i < BSNode::branchingNodeIds.size(); i++)
        BST->V[BSNode::branchingNodeIds[i]].setPEMatchesNull();

    BST->cleanTree(chosenBsNode, headBSNode);

    updateHistoryMap(chosenBsNode);
    headBSNode = chosenBsNode;
    BST->headBSNode = headBSNode; //required for dot visualization
    return 1;
}

int BackwardsWalker::elongate(vector<unsigned int> &path) {
    unsigned int ogNodeId;
    unsigned int maxJump = P->getGlobalMaxAllowedDistance();
    unsigned int previousHead = 1;
    unsigned int lastNew = 1;
    unsigned int BSNode;
    unsigned int chosenBsNode;


    bool PRINT_TREE = false;
    unsigned int val;
    unsigned int maxD;
    unsigned int nStep = 0;
    unsigned int nStepLimit = 50;
    unsigned int explorationDepth;

    //    bool reDo = true; //debug purpose
    //    while (reDo) {
    //    reDo = false;

//    if (maxJump == 0)
//        continue;

    maxD = BST->V[headBSNode].getDistance() + maxJump;
    
    do {
        previousHead = headBSNode;

        OverlapsGraph::AP.printProgress(cout);
        val = stepElongate(maxD, chosenBsNode); //could update headBSNode
        //val==2: dead end
        //val==4: search canceled
        //val==5: non resolvable bubble , hopeless
        //val==0: tree updated
        //chosenBsNode = 0: undetermined
        //chosenBsNode = uintmax: libraries conflict

        //nucleotide distance
        explorationDepth = BST->getUsableDistance() - BST->V[headBSNode].getDistance();

        //reject chosen node when exploration depth smaller 15 nt.
        if (explorationDepth < 15 && val != 5) //val==5 means a bubble -> hopeless to be resolved later
        {
            chosenBsNode = 0; //forces a deeper exploration to ensure a decent sampling
        }

        if (chosenBsNode == 0 || chosenBsNode == numeric_limits<unsigned int>::max())
            nStep++;
        else
            nStep = 0;

        if (PRINT_TREE) {
            dot("test");
        }

        if (chosenBsNode == numeric_limits<unsigned int>::max()) //conflict
        {
            if (PRINT_TREE)
                dot("test");
        } else if (chosenBsNode > 0) //elongation: update headBSNode
        {//prune tree, update history map, update graph flags

            //snapshot
            //                st_count++;
            //                ostringstream oss;
            //                oss << "snapshot/" << setfill('0') << setw(4) << st_count << "_nstep" << nStep;
            //                dot (oss.str(),"png");
            //                nStep=0;

            for (unsigned int i = 0; i < BSNode::branchingNodeIds.size(); i++)
                BST->V[BSNode::branchingNodeIds[i]].setPEMatchesNull();

            BST->cleanTree(chosenBsNode, headBSNode);

            updateHistoryMap(chosenBsNode);
            headBSNode = chosenBsNode;
            BST->headBSNode = headBSNode; //required for dot visualization

            NodeIt nIt;
            nIt.initNodeIt(BST->V[headBSNode].getNodeId(), BST->V[headBSNode].getNodeDir());
        }

        maxD = BST->V[headBSNode].getDistance() + maxJump;

        if (PRINT_TREE)
            dot("test");

        if (headBSNode != previousHead) {//has been elongated

            BSNode = previousHead;

            do {
                ogNodeId = BST->V[BSNode].getOGNodeId();
                // ogNodeDir = BST->V[BSNode].getOGNodeDir();

                if (G->nodesTab[ogNodeId].isAlreadyUsed() == false) {
                    unsigned int add = BST->V[BSNode].getDistance() - BST->V[lastNew].getDistance() + 1;
                    OverlapsGraph::AP.addKb(add);
                    OverlapsGraph::AP.printProgress(cerr, true);

                    lastNew = BSNode;
                    G->nodesTab[ogNodeId].setAlreadyUsed();
                }

                BSNode = BST->V[BSNode].getSon();

            } while (BSNode != 0);
        } else
            if (BST->getUsableDistance() > maxD) {
            //out of paired-end range
            //non-resolvable ambiguity
            val = 3;
        }

        //redundancy check
        if (BST->V[headBSNode].getDistance() - BST->V[lastNew].getDistance() > maxJump) {
            //collision
            //enter an already assembled path
            val = 1;
        } else
            if (nStep >= nStepLimit) {
            nStep = 0;
            val = 4;
        }
    } while (val == 0);

    //        if (PRINT_TREE)
    //            dot("test");

    BST->cleanTree(headBSNode);

    // } //end redo

    
   
    if (lastNew == headBSNode)
        BST->getPath(lastNew, path);
    else {
         //include redundant ending path up to 'maxRedundany'
        unsigned int bsNode = lastNew;
        unsigned int nextNode = BST->V[bsNode].getSon();

        while (nextNode != 0) {
            if (BST->V[nextNode].getDistance() - BST->V[lastNew].getDistance() + 1 > maxRedundancy)
                break;

            bsNode = nextNode;
            nextNode = BST->V[bsNode].getSon();
        }

        BST->getPath(bsNode, path);
    }

    //    if (lastNew > BST->V[headBSNode].getFather())
    //        BST->getPath(lastNew, path);
    //    else
    //        BST->getPath(BST->V[headBSNode].getFather(), path);

    //1 collision
    //2 search tree reaches a dead-end
    //3 out of paired-end range
    //4 step limit reached
    //5 non resolvable bubble

    return val;
}

int BackwardsWalker::stepElongate(unsigned int maxD, unsigned int &BSNode) {
    BSNode = 0;
    unsigned int maxTreeD = maxD; //!!
    unsigned int ret;
    bool PRINT_TREE = false;
    bool firstUpdate = BST->getQueueAt(0) == headBSNode;

    //OverlapsGraph::AP.printProgress(cout);

    ret = BST->stepUpdateTree(maxTreeD);
    //ret=0 : canceled (max tree size reached)
    //ret=uintmax : dead-end

    if (PRINT_TREE) {
        dot("test");
    }

    if (firstUpdate) //first update
    {
        nPotentialElongation = BST->getQueueSize();
        BSNode::branchingNodeIds.clear();
        BSNode::nReadsChecked.clear();
        BST->initPEMatchesStorage();

        //init the PEM that sums up all pe matching information
        BST->summedPEM = BST->PEMS->getPPEMatches(BST->getQueueSize());

        for (size_t i = 0; i < BST->getQueueSize(); i++) {
            unsigned int bsn = BST->getQueueAt(i);
            BST->V[bsn].setBranchingNodeIndex(i);
            BSNode::branchingNodeIds.push_back(bsn);
            BSNode::nReadsChecked.push_back(0);

            if (BST->V[bsn].PEMatchesAllocated() != 0) {
                cout << "problem";
                exit(0);
            }
            BST->V[bsn].initPEMatches(1);
        }

        //init activeBsNode
        for (size_t i = 0; i < activeBsNodes.size(); i++) {
            BST->V[ activeBsNodes[i] ].setInactive();
            BST->V[ activeBsNodes[i] ].setPEMatchesNull();
        }
        activeBsNodes.clear();

    } else {
        unsigned int id;
        bool dir;

        unsigned int bsn;

        for (size_t i = 0; i < BST->getQueueSize(); i++) {
            bsn = BST->getQueueAt(i);
            unsigned int father = BST->V[bsn].getFather();

            if (father != headBSNode) {
                BST->V[bsn].setBranchingNodeIndex(BST->V[father].getBranchingNodeIndex());
            }
        }

        //check for unresolved bubble
        // ->hopeless to try past the bubble 
        unsigned int c = 0;
        for (size_t i = 0; i < BST->getQueueSize(); i++) {
            bsn = BST->getQueueAt(i);
            if (i == 0) {
                id = BST->V[bsn].getNodeId();
                dir = BST->V[bsn].getNodeDir();
            } else {
                if (BST->V[bsn].getNodeId() == id && BST->V[bsn].getNodeDir() == dir)
                    c++;
            }

        }
        if (c == BST->getQueueSize() - 1 && BST->getQueueSize() >= 2) {
            //checkDistances

            unsigned int d, maxd = 0, mind = 1E8;

            for (size_t i = 0; i < BST->getQueueSize(); i++) {
                bsn = BST->getQueueAt(i);
                d = BST->V[bsn].getDistance();
                if (d > maxd)
                    maxd = d;
                if (d < mind)
                    mind = d;
            }
            if (maxd - mind < 5) //supposed non-resolvable
            {
                BSNode = 0;
                return 5; //non resolvable bubble
            }
        }
    }

    //map the current interval to history
    backwardsMap();

    BSNode = getDecision();
    // BSNode = BST->getDecision(activeBsNodes, minNPair, minRatio);
    BST->chosen = BSNode;

    //return value: 1=collision, 2=dead-end, 3=ambiguity, 4=canceled
    if (ret == numeric_limits<unsigned int>::max()) { //dead-end

        return 2;
    }


    if (ret == 0) {//canceled : max beam width reached

        return 4;
    }

    //return value: 2=dead-end, 4=canceled, 5=bubble ambiguity
    return 0;
}

unsigned int BackwardsWalker::getDecision() {
    PEMatches *pem;
    vector<BSNode>::iterator pNode;
    unsigned int chosen = 0;
    unsigned int nHit, sumHit, maxHit;
    size_t maxIndex;
    unsigned int nChoice = BSNode::branchingNodeIds.size();
    static vector<unsigned int> sumUp;
    static vector<unsigned int> sumUpTmp;
    static vector<unsigned int> hits;
    static vector<double> expNHits;
    sumUp.assign(nChoice, 0);
    hits.assign(nChoice, 0);

    for (unsigned int i = 0; i < activeBsNodes.size(); i++) {
        pNode = BST->V.begin() + activeBsNodes[i];
        pem = pNode->ppeMatches;
        pNode->flag = false;

        sumUpTmp.assign(nChoice, 0);

        for (unsigned int lib = 1; lib <= PEMatches::nLibrary; lib++) {
            //currentChoice = 0;
            BST->computeExpectedNHits(lib, activeBsNodes[i], expNHits);

            sumHit = maxHit = 0.0;
            for (unsigned int index = 0; index < nChoice; index++) {
                nHit = pem->getCount(lib, index);

                //safeguard: limit nHit to 2*expected bHits. Expected nHits is
                //an estimation of the number of hits expected given sampling range and a
                //coverage of 1. This is done to avoid bias caused by sequences repeated many times
                //which can cause mis-assemblies

                if (nHit > expNHits[index]*2.0)
                    nHit = expNHits[index]*2.0;

                hits[index] = nHit;

                if (nHit >= maxHit) {
                    maxHit = nHit;
                    maxIndex = index;
                }
                sumHit += nHit;
            }

            if (sumHit == 0)
                continue;

            if (maxHit >= minNPair && double(maxHit) / sumHit >= minRatio) {
                //mark node for visualization
                pNode->flag = true;

                for (unsigned int index = 0; index < nChoice; index++) {
                    sumUpTmp[index] += hits[index];
                }

            }
        }

        for (unsigned int index = 0; index < nChoice; index++)
            sumUp[index] += sumUpTmp[index];
    }

    maxHit = 0;
    sumHit = 0;

    for (unsigned int i = 0; i < nChoice; i++) {
        if (sumUp[i] > maxHit) {
            maxHit = sumUp[i];
            maxIndex = i;
        }
        sumHit += sumUp[i];
    }

    if (sumHit == 0)
        return 0; //undetermined

    if ((double) maxHit / sumHit >= minRatio) {
        chosen = BSNode::branchingNodeIds[maxIndex];
    } else
        return numeric_limits<unsigned int>::max(); //libraries conflict 

    return chosen;
}

void BackwardsWalker::CountPEMatches(vector<unsigned int> path, size_t p1, size_t p2) {
    //check path
    vector<unsigned int> P1P2;
    P1P2.insert(P1P2.begin(), path.begin(), path.begin() + p2);

    //init beamSearchTree with P1-P2
    BeamSearchTree T(P);
    T.init(P1P2);


    //index P1
    updateHistoryMap(p1);
    //init PEmatch

    //go through p3 and map reads

    //output PEMatch
}

int BackwardsWalker::updateHistoryMap(unsigned int BSNode) {
    BsNodeIdentifier bsNId;
    unsigned int OGNodeId;
    bool OGNodeDir;
    unsigned int currentBsNode = headBSNode;

    do {
        currentBsNode = BST->V[currentBsNode].getSon();

        OGNodeId = BST->V[currentBsNode].getNodeId();
        OGNodeDir = BST->V[currentBsNode].getNodeDir();

        bsNId.setNodeId(OGNodeId, OGNodeDir);
        bsNId.BSNode = currentBsNode;
        historyMap.insert(bsNId);
        G->nodesTab[OGNodeId].setVisited(OGNodeDir);
    } while (currentBsNode != BSNode);
    return 0;
}

void BackwardsWalker::backwardsMap() {

    bool noStepLimit = false;
    int lib;

    unsigned int maxSearchD = BST->getUsableDistance(); //step limit
    unsigned int maxJump = P->getGlobalMaxAllowedDistance();
    int mateOrientation = 1;

    if (maxSearchD > BST->V[headBSNode].getDistance() + maxJump - rL)
        maxSearchD = BST->V[headBSNode].getDistance() + maxJump - rL;

    BsNodeIdentifier bsNId;
    size_t BSNodeInHistory;
    int mo;
    unsigned int min, max;
    bool requiredSourceReadDir;
    bool pairedNodeDir;
    unsigned int (ReadsLayout::*getNext)(size_t) const;
    unsigned int d1, d2, actualDistance;
    unsigned int searchD;
    unsigned int backwardNodeDistance;

    //nPotentialElongation is the number of nodes under headBsNode



    //                 <<<<<<<<<<<<<<<<<<<<<<< backwards search   root <<<< leaves            
    //     pairedNode                               sourceNode
    // (          -->       )...............(      <--           )
    //            |..d1.....|                         |....d2....|

    //for all enqueued node


    for (size_t queueIndex = 0; queueIndex < BST->getQueueSize(); queueIndex++) {
        size_t sourceBSNode = BST->getQueueAt(queueIndex);
        unsigned int bIndex = BST->V[sourceBSNode].getBranchingNodeIndex(); //bIndex equal queueIndex??
        size_t sourceNodeId = BST->V[sourceBSNode].getNodeId();
        bool sourceNodeDir = BST->V[sourceBSNode].getNodeDir();
        unsigned int sourceNodeDistance = BST->V[sourceBSNode].getDistance();
        unsigned int sourceNodeLength = G->nodesTab[sourceNodeId].getSequenceLength();



        unsigned int currentRead = BST->V[sourceBSNode].getPLayout();

        if (sourceNodeDir)
            getNext = &ReadsLayout::getNext;
        else
            getNext = &ReadsLayout::getPrevious;

        unsigned int nReadsChecked = 0; //used to correct bias due to difference in nodes coverage
        unsigned int nMatchPE = 0;

        //        if (BST->V[sourceBSNode].PEMatchesAllocated() == 0)
        //            BST->V[sourceBSNode].initPEMatches(1);

        for (; currentRead != 0; currentRead = (L->*getNext)(currentRead)) {
            //            if (currentRead%1000==0)
            //             OverlapsGraph::AP.printProgress(cout);
            //first check distance
            if (sourceNodeDir)
                d2 = sourceNodeLength - (L->getPosition(currentRead) + rL - 1);
            else
                d2 = L->getPosition(currentRead) - 1;

            searchD = sourceNodeDistance - d2;

            if (noStepLimit == false) {
                if (searchD > maxSearchD)
                    break;
            }

            nReadsChecked++;

            //then check whether read is paired
            size_t pairedRead = P->getPairing(currentRead);
            if (pairedRead == 0)
                continue;


            mateOrientation = P->getMateOrientation(currentRead);


            if (mateOrientation == 1)
                requiredSourceReadDir = !sourceNodeDir;
            else
                requiredSourceReadDir = sourceNodeDir;

            if (L->getDirection(currentRead) != requiredSourceReadDir)
                continue;


            size_t pairedNode = L->getNodeId(pairedRead);

            //todo: keep distances locally
            P->getDistanceRange(currentRead, min, max, mo);

            //            if (mo != mateOrientation)
            //                continue;



            bool pairedDir = L->getDirection(pairedRead);

            if (mateOrientation == 1)
                pairedNodeDir = pairedDir;
            else
                pairedNodeDir = !pairedDir;

            //Does paired node belong to history ?
            if (G->nodesTab[pairedNode].isVisited(pairedNodeDir) == 0)
                continue;

            //get the corresponding(s) BSnodes
            pair<multiset<BsNodeIdentifier>::iterator, multiset<BsNodeIdentifier>::iterator> ret;
            multiset<BsNodeIdentifier>::iterator it;

            bsNId.setNodeId(pairedNode, pairedNodeDir);

            ret = historyMap.equal_range(bsNId); //should contain a single instance most of the time

            if (ret.first == ret.second) {
                cout << "void BackwardsWalker::backwardsMap(int mateOrientation)\n";
                cout << "This should not happen\n";
                exit(0);
            }

            for (it = ret.first; it != ret.second; ++it) {//most of the time only a single iteration

                BSNodeInHistory = it->getBSNode();
                backwardNodeDistance = BST->V[BSNodeInHistory].getDistance();
                //check distance for BSNode it

                if (pairedNodeDir)
                    d1 = G->nodesTab[pairedNode].getSequenceLength() - L->getPosition(pairedRead) + 1;
                else
                    d1 = L->getPosition(pairedRead) + rL - 1;

                actualDistance = sourceNodeDistance - backwardNodeDistance + d1 - d2;


                if ((actualDistance >= min && actualDistance <= max) || !checkDistances) {//pair connected within allowed range
                    nMatchPE++;

                    if (BST->V[BSNodeInHistory].isActive() == false) {
                        activeBsNodes.push_back(BSNodeInHistory);
                        BST->V[BSNodeInHistory].setActive();
                        BST->V[BSNodeInHistory].initPEMatches(nPotentialElongation);
                    }
                    lib = P->getPeLibraryID(currentRead);
                    //    unsigned int bIndex=BST->V[sourceBSNode].getBranchingNodeIndex();
                    BST->V[BSNodeInHistory].addDistanceSample(lib, bIndex, actualDistance);
                    unsigned int n = BSNode::branchingNodeIds.at(bIndex);
                    BST->V[n].addDistanceSample(lib, 0, actualDistance);
                    BST->summedPEM->addMatch(lib, bIndex, actualDistance);
                }
            }

        } // end for current read
        BSNode::nReadsChecked[bIndex] += nReadsChecked;

        BST->V[sourceBSNode].setPLayout(currentRead); //ready for next round
    } //end for queue
}

void BackwardsWalker::getPath(unsigned int treeNode, vector<unsigned int> &p) {
    BST->getPath(treeNode, p);
}

void BackwardsWalker::getPath(vector<unsigned int> &p) {
    BST->getPath(headBSNode, p);
}

