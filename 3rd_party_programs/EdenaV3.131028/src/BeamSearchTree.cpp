/* 
 * File:   BeamSearchTree.cpp
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

#include "BeamSearchTree.h"
#include "globalFunc.h"
#include <limits>
using namespace std;

unsigned int BeamSearchTree::maxAchievedBeamWidth = 0;
vector<unsigned int> BSNode::markedNode;
vector<unsigned int> BSNode::branchingNodeIds;
vector<unsigned int> BSNode::nReadsChecked;
PEMatchesStorage * BSNode::PEMS = 0x0;
extern ofstream outGlob;

BeamSearchTree::BeamSearchTree(Pairing *p) {

    P = p;
    PEMS = new PEMatchesStorage(p);
    BSNode::PEMS = PEMS;
}

BeamSearchTree::BeamSearchTree(const BeamSearchTree& orig) {
}

BeamSearchTree::~BeamSearchTree() {

    delete PEMS;
}

void BeamSearchTree::init(NodeIt n) {
    BSNode bsn;
    V.clear();
    nodeQueue.clear();
    V.push_back(bsn);

    bsn.init(n);
    V.push_back(bsn);
    nodeQueue.push_back(1); //root
    usableDistance = 0;
    headBSNode = 1;
}

void BeamSearchTree::init(vector<unsigned int> &path) {
    BSNode bsn;
    NodeIt n;
    bool initOrientation;
    V.clear();
    nodeQueue.clear();
    V.push_back(bsn);
    headBSNode = 0;

    if (path.size() < 2)
        return;

    if (path[0] > 0)
        initOrientation = true;
    else
        initOrientation = false;

    n.initNodeIt(path[1], initOrientation);
    bsn.init(n);
    V.push_back(bsn); //root

    size_t p = 2;

    unsigned int d = 0;

    while (p < path.size()) {
        //n.getNeighbor(path[p],n);
        n.getNeighborGivenAbsoluteIndex(path[p], n);
        bsn.init(n);
        d += bsn.nodeLength - bsn.ovSize;
        bsn.distance = d;
        V.push_back(bsn);

        V[V.size() - 2].son = V.size() - 1;
        V[V.size() - 1].father = V.size() - 2;

        p++;
    }

    usableDistance = 0;
    headBSNode = V.size() - 1;
    nodeQueue.push_back(headBSNode);
    // dot("test","test");
}

unsigned int BeamSearchTree::stepUpdateTree(unsigned int maxTreeDepth) {
    size_t n, last;
    size_t qSize;
    NodeIt nIt, nextNIt;
    BSNode bsNode;
    nodeToMap.clear();

    unsigned int newUsableDistance = numeric_limits<unsigned int>::max();

    //possibly assign newUsableDistance here
    //when the next step does not require to add
    //a new node.
    if (nodeQueue.size() > 1)//when some node in queue have been pruned
    {
        for (unsigned int i = 0; i < nodeQueue.size(); i++) {
            n = nodeQueue[i];

            if (V[n].distance > usableDistance) {
                if (V[n].distance < newUsableDistance)
                    newUsableDistance = V[n].distance;
            }
        }
    }

    qSize = nodeQueue.size();

    for (unsigned int i = 0; i < qSize; i++) {
        if (nodeQueue.size() > maxAchievedBeamWidth)
            maxAchievedBeamWidth = nodeQueue.size();

        if (nodeQueue.size() > maxBeamWidth) {
            //give up
            return 0;
            //keeps usableDistance to previous state
        }

        //get node from queue
        n = nodeQueue.front();
        nodeQueue.pop_front();

        if (V[n].distance >= newUsableDistance ||
                (maxTreeDepth != 0 && V[n].distance > maxTreeDepth)) {
            nodeQueue.push_back(n); //keep for next round
        } else { //add one level below V[n]

            nIt.initNodeIt(V[n].nodeId, V[n].nodeDir);
            nIt.getNext(nextNIt); //get first son

            if (nextNIt.isNull() == false) {//connect father to first son
                V[n].son = V.size(); //first son
            }

            last = 0;
            while (nextNIt.isNull() == false) { //add remaining sons
                bsNode.init(nextNIt);
                bsNode.distance = V[n].distance + nextNIt.getNodeLength() - nextNIt.getArrivalEdgeSize();
                V.push_back(bsNode);

                if (bsNode.distance < newUsableDistance)
                    newUsableDistance = bsNode.distance;

                if (last != 0) {
                    V[last].brother = last + 1;
                    last++;
                } else
                    last = V.size() - 1;

                V[last].father = n;
                nodeQueue.push_back(last);

                nIt.getNext(nextNIt);
            }
        }
    } // end for all i in queue


    if (newUsableDistance != numeric_limits<unsigned int>::max()) {
        usableDistance = newUsableDistance;
    }

    //return 0 if canceled (maxQueue size reached)
    //uintmax if cannot be updated (queue empty)

    return newUsableDistance;
}

void BeamSearchTree::getPath(unsigned int leaf, vector<unsigned int> &pp) {

    static vector<unsigned int> tmpPath;
    tmpPath.clear();

    pp.clear();
    if (V[1].getNodeDir() == true)
        pp.push_back(1);
    else
        pp.push_back(0);


    while (leaf != 0) {
        tmpPath.push_back(V[leaf].edgeIndex);
        leaf = V[leaf].father;
    }

    for (size_t i = tmpPath.size(); i > 0; i--)
        pp.push_back(tmpPath[i - 1]);
}

void BeamSearchTree::updateQueue(unsigned int node) {
    nodeQueue.clear();
    updateQueue_(node);
}

void BeamSearchTree::updateQueue_(unsigned int node) {
    if (V[node].son == 0)
        nodeQueue.push_back(node);
    else
        updateQueue_(V[node].son);

    if (V[node].brother != 0)
        updateQueue_(V[node].brother);
}

void BeamSearchTree::cleanTree(unsigned int currentHead) {
    V[currentHead].son = 0;
    nodeQueue.clear();
    nodeQueue.push_back(currentHead);
}

int BeamSearchTree::cleanTree(unsigned int n, unsigned int currentHead) {

    unsigned int chosenNode = n;
    if (chosenNode <= currentHead) {
        cout << "int BeamSearchTree::cleanTree(...) problem\n";
        exit(0);
    }

    unsigned int father;


    V[chosenNode].son = 0;


    while (chosenNode != currentHead) { //kill all brothers
        V[chosenNode].brother = 0;
        father = V[chosenNode].father;
        V[father].son = chosenNode;
        chosenNode = father;
    }

    nodeQueue.clear();
    nodeQueue.push_back(n);

    return 0;
}

void BeamSearchTree::dot(string fileName, string title) {
    dot(fileName, title, "ps");
}

void BeamSearchTree::dot(string fileName, string title, string format) {
    
    ofstream out((fileName + ".dot").c_str());
    out << "digraph " << title << " {\n";
    //out << "size=\"8,10\";" << endl;
    out << "node [shape = ellipse,fontsize = 10];\n";

    if (V.size() >= 2)
        dot(out, 1);

    out << "}\n";
    out.close();
    out.clear();

    string command = "dot -T" + format + " -o ";
    // string command = "dot -Tps -o ";
    command += fileName + '.' + format + " " + fileName + ".dot";
    //cout << "executing " << command << " ... " << flush;
    int unused __attribute__((unused)); //get rid of gcc warning
    unused = system(command.c_str());
    cerr << "File " << fileName + '.' + format << " written" << endl;

}

void BeamSearchTree::dot(ostream &out, size_t node) {
    unsigned int son = V[node].son;
    NodeIt nIt;
    nIt.initNodeIt(V[node].nodeId, V[node].nodeDir);
    // static vector<double> meanShifts;
    static vector<double> summedProb;
    static vector<double> eMean;

    string color = "black";
    string shape = "ellipse";
    string style = "solid";


    if (V[node].active == true) {
        if (V[node].flag == true) {
            color = "orange";
            style = "filled";
        } else {
            color = "cadetblue1";
            style = "filled";
        }
    }

    if (node == headBSNode)
        shape = "hexagon";

    if (node == chosen) {
        color = "green";
    }

    out << node << "[style=" << style << ", color=" << color << ", shape=" << shape << ",";

    out << "label=\"";
    out << node << " d=" << V[node].distance << "\\n";


    if (V[node].nodeDir == false)
        out << "(!";
    else
        out << '(';
    out << V[node].nodeId << ')'
            << " c=" << nIt.getCoverage()
            << " l=" << V[node].nodeLength << "\\n";

    if (V[node].ppeMatches != 0x0 && node <= headBSNode) {
        for (size_t lib = 1; lib <= PEMatches::nLibrary; lib++) {
            computeExpectedNHits(lib, node, summedProb);
            computeExpectedMeans(lib, node, eMean);

            if (V[node].ppeMatches->getSumCount(lib) > 0) {
                out << "lib" << lib << ": ";
                //counts
                for (size_t i = 0; i < V[node].ppeMatches->getNChoice(); i++) {
                    out << V[node].ppeMatches->getCount(lib, i);
                    out << '(' << setprecision(4) << summedProb[i] << ')';

                    if (V[node].ppeMatches->getNChoice() - i > 1)
                        out << ", ";
                    else
                        out << "\\n";
                    //out << " | ";
                }
                //                for (size_t i = 0; i < V[node].ppeMatches->getNChoice(); i++)
                //                {
                //                    out << V[node].ppeMatches->getMean(lib, i);
                //                    out << '(' << setprecision(4) << eMean[i] << ')';
                //                    
                //                    if (V[node].ppeMatches->getNChoice() - i > 1)
                //                        out << ", ";
                //                    else
                //                        out << "\\n";
                //                        //out << " * ";
                //                }
            }
        }
    } else if (V[node].getFather() == headBSNode) {
        for (size_t lib = 1; lib <= PEMatches::nLibrary; lib++) {

            //has to be adapted
            // computeExpectedNHits(lib,node,summedProb);
            computeExpectedMeans(lib, eMean);

            //   if (V[node].ppeMatches->getSumCount(lib) > 0)
            if (summedPEM->getSumCount(lib) > 0) {
                out << "lib" << lib << ": ";
                //counts
                for (size_t i = 0; i < summedPEM->getNChoice(); i++) {
                    out << summedPEM->getCount(lib, i);
                    out << "()";
                    //  out << '(' << setprecision(4) << summedProb[i] << ')';

                    if (summedPEM->getNChoice() - i > 1)
                        out << ", ";
                    else
                        //out << "\\n";
                        out << " | ";
                }
                for (size_t i = 0; i < summedPEM->getNChoice(); i++) {
                    out << summedPEM->getMean(lib, i);
                    //out << "()";
                    out << '(' << setprecision(4) << eMean[i] << ')';

                    if (summedPEM->getNChoice() - i > 1)
                        out << ", ";
                    else
                        out << "\\n";
                    //out << " * ";

                }
            }
        }

    }

    //mean distance
    //                for (size_t i = 0; i < V[node].ppeMatches->getNChoice(); i++)
    //                {
    //                    if (V[node].ppeMatches->getCount(lib, i) == 0)
    //                        out << "---";
    //                    else
    //                        out << V[node].ppeMatches->getMean(lib, i);
    //                    if (V[node].ppeMatches->getNChoice() - i > 1)
    //                        out << ", ";
    //                    else
    //                        out << "\\n";
    //                }

    //                unsigned int maxIndex;
    //
    //                if (V[node].active == true)
    //                { //in history
    //
    //                    //relative mean shifts
    //                    V[node].ppeMatches->computeRelativeMeanShifts(lib, maxIndex, meanShifts);
    //                    for (size_t i = 0; i < V[node].ppeMatches->getNChoice(); i++)
    //                    {
    //                        if (V[node].ppeMatches->getCount(lib, i) == 0)
    //                            out << "---";
    //                        else
    //                            out << meanShifts.at(i);
    //                        if (V[node].ppeMatches->getNChoice() - i > 1)
    //                            out << ", ";
    //                        else
    //                            out << "\\n";
    //                    }
    //                }
    //                
    //                else
    //                { //in decision zone
    //                     summedPEM->computeRelativeMeanShifts(lib,maxIndex, meanShifts);
    //                     out << meanShifts.at(V[node].branchingNodeIndex) << "\\n";
    //                }
    //                
    //absolute mean shifts
    //                double s_mean, s_sd, s_nEchant;
    //                
    //                int nPosA=getUsableDistance()-V[headBSNode].getDistance();
    //                int nPosB=V[node].getNodeLength();
    //                int spacer=V[headBSNode].getDistance()-V[node].getDistance()-53;
    //                V[node].ppeMatches->computeAbsoluteMeanShifts(lib,
    //                                                              meanShifts,
    //                                                              s_mean,
    //                                                              s_sd,
    //                                                              s_nEchant,
    //                                                              nPosA,
    //                                                              nPosB,
    //                                                              spacer);
    //                for (size_t i = 0; i < V[node].ppeMatches->getNChoice(); i++)
    //                {
    //                    if (V[node].ppeMatches->getCount(lib, i) == 0)
    //                        out << "---";
    //                    else
    //                        out << meanShifts.at(i) << "(" << s_mean<< ")";
    //                    if (V[node].ppeMatches->getNChoice() - i > 1)
    //                        out << ", ";
    //                    else
    //                        out << "\\n";
    //                }

    //summedProb
    //                static vector<double> summedProb;
    //                 static vector<double> eMean;
    //               computeSummedProbs(lib,node,summedProb);
    //                computeExpectedMeans(lib,node,eMean);
    //                for (size_t i = 0; i < V[node].ppeMatches->getNChoice(); i++)
    //                {
    //                   
    //                        out << "p" << summedProb.at(i) << "e" << eMean.at(i);
    //                    if (V[node].ppeMatches->getNChoice() - i > 1)
    //                        out << ", ";
    //                    else
    //                        out << "\\n";
    //                }

    //                for (size_t i = 0; i < V[node].ppeMatches->getNChoice(); i++)
    //                {
    //                    if (V[node].ppeMatches->getCount(lib, i) == 0)
    //                        out << "---";
    //                    else
    //                        out << V[node].ppeMatches->getSDMeanShift(lib, i);
    //                    if (V[node].ppeMatches->getNChoice() - i > 1)
    //                        out << ", ";
    //                    else
    //                        out << "\\n";
    //                }
    //            }
    //        }
    //    }

    out << "\"]\n";

    while (son != 0) {
        out << node << " -> " << son;
        out << " [arrowsize = 0.7,";
        out << "arrowtail = inv,";
        out << "arrowhead = none];\n";

        dot(out, son);
        son = V[son].brother;
    }
}

//estimate the expected number of PE hits in the node given
//the search zone depth (shift)
//and assuming a mate coverage of 1.0

void BeamSearchTree::computeExpectedNHits(unsigned int lib,
        unsigned int bsNode,
        vector<double> &summedProbs) {
    unsigned int nChoice = BSNode::branchingNodeIds.size();
    summedProbs.clear();
    summedProbs.assign(nChoice, 0.0);
    double summed_p;

    for (unsigned int i = 0; i < nChoice; i++) {
        unsigned a = V[headBSNode].getDistance() - V[bsNode].getDistance() + 1;
        a += 2 * (PELibrary::readLength - 1); //need to 
        //substract ovlp at headbsnode
        a -= V[BSNode::branchingNodeIds[i]].ovSize;

        unsigned b = a + V[bsNode].getNodeLength() - PELibrary::readLength + 1;
        unsigned int shift = getUsableDistance() - V[headBSNode].getDistance();
        // shift-= PELibrary::readLength-1;

        // summed_p = P->getLibraryIt(lib-1)->getSummedProb(a, b, shift);
        summed_p = P->getLibraryIt(lib - 1)->lengthDistr.sumCdf(a - 1, b, shift);
        summedProbs[i] = summed_p;
        summedProbs[i] = summedProbs[i] * P->getLibraryIt(lib - 1)->getExpectedMateCoverage();
    }
}

void BeamSearchTree::computeExpectedMeans(unsigned int lib,
        unsigned int bsNode,
        vector<double> &eMean) {
    unsigned int nChoice = BSNode::branchingNodeIds.size();
    eMean.clear();
    eMean.assign(nChoice, 0.0);
    double eM;

    for (unsigned int i = 0; i < nChoice; i++) {
        unsigned a = V[headBSNode].getDistance() - V[bsNode].getDistance() + 1;
        a += 2 * (PELibrary::readLength - 1); //need to 
        //substract ovlp at headbsnode
        a -= V[BSNode::branchingNodeIds[i]].ovSize;

        unsigned b = a + V[bsNode].getNodeLength();
        unsigned int shift = getUsableDistance() - V[headBSNode].getDistance();
        // shift-= PELibrary::readLength-1;

        eM = P->getLibraryIt(lib - 1)->lengthDistr.mu(a, b, shift);
        eMean[i] = eM;
    }
}

//estimate the expected mean distance for the total mates
//of a given lib mapped in the tree

void BeamSearchTree::computeExpectedMeans(unsigned int lib,
        vector<double> &eMean) {
    unsigned int nChoice = BSNode::branchingNodeIds.size();
    eMean.clear();
    eMean.assign(nChoice, 0.0);
    double eM;

    for (unsigned int i = 0; i < nChoice; i++) {
        //unsigned a = V[headBSNode].getDistance() - V[bsNode].getDistance()  +1;

        unsigned int a = 0;
        a += 2 * (PELibrary::readLength - 1); //need to 
        //substract ovlp at headbsnode
        a -= V[BSNode::branchingNodeIds[i]].ovSize;

        //unsigned b = a + V[bsNode].getNodeLength();

        unsigned int b = a + V[headBSNode].getDistance();
        b += V[1].getNodeLength();
        b -= PELibrary::readLength - 1; //must get real overlap here

        unsigned int shift = getUsableDistance() - V[headBSNode].getDistance();
        // shift-= PELibrary::readLength-1;

        eM = P->getLibraryIt(lib - 1)->lengthDistr.mu(a, b, shift);
        eMean[i] = eM;
    }
}

BSNode::BSNode() {
    init();
}

void BSNode::init() {

    edgeIndex = 0;
    nodeId = 0;
    nodeDir = true;
    distance = 0; //sequence length from the root
    son = 0;
    brother = 0;
    father = 0;
    flag = true;
    pLayout = 0;
    ppeMatches = 0x0;
    active = false;
}

void BSNode::init(NodeIt &n) {
    init();
    edgeIndex = n.getAbsArrivalEdgeIndex();
    ovSize = n.getArrivalEdgeSize();
    nodeId = n.getNodeId();
    nodeDir = n.getDirection();
    nodeLength = n.getNodeLength();
    if (nodeDir)
        pLayout = n.getFirstReadInLayout();
    else
        pLayout = n.getLastReadInLayout();

    ppeMatches = 0x0;
}

void BSNode::getNodeIt(NodeIt& nit) const {
    nit.initNodeIt(nodeId, nodeDir);
}

void BSNode::addDistanceSample(unsigned int lib, unsigned int branching, unsigned int d) {
    if (ppeMatches == 0x0) {
        cout << "void BSNode::addDistanceSample(...) problem\n";
        sendBugReportPlease(cout);
    }
    ppeMatches->addMatch(lib, branching, d);
}

