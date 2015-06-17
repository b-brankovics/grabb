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

#include <algorithm>
using namespace std;

#include "node.h"
#include "readsStorage.h"
//#include "layouts.h"
#include "globalFunc.h"
#include "NodeIt.h"
#include "stat.h"
#include "overlapGraph.h"

extern ofstream outGlob;
extern bool DEV_FLAG;
extern unsigned long global_count;
Node* Node::N = 0x0;
OverlapsGraph* Node::G = 0x0;
ReadsStorage * Node::R = 0x0;
ReadsLayout * Node::L = 0x0;
Pairing * Node::P = 0x0;
unsigned int Node::S_counter;
unsigned int Node::nNodes;
vector<unsigned int> Node::path;
double Node::pathLengthMean = 0;
double Node::pathLengthSD = 0;
unsigned int Node::nPathFound = 0;
vector<unsigned int> Node::nodeList;
vector<unsigned int> Node::dotNodeList;

deque<NodeMem> Node::stack;
bool Node::edgeSorted = false;

unsigned int Node::maxLog = 0;

Node::Node() {
    nOv = 0;
    nOvRight = 0;
    layout = 0;

    ovId = 0x0;
    ovSize = 0x0;
    edgeValue = 0x0;
    edgeBitsFlag = 0x0;
    value = 0;
    flagsBit &= 0; //unset all flags
}

Node::~Node() {
    freeMemory();
}

void Node::allocate(int n) {
    if (ovId != 0x0) {
        free(ovId);
        //delete[] ovId;
        ovId = 0x0;
    }
    if (ovSize != 0x0) {
        free(ovSize);
        //delete[] ovSize;
        ovSize = 0x0;
    }

    if (n > 0) {
        nOv = n;

        ovId = (unsigned int*) malloc(n * sizeof (unsigned int));
        ovSize = (short*) malloc(n * sizeof (short));
    }

    allocateEdgeFlag();
}

void Node::allocateEdgeValue() {
    if (edgeValue != 0x0) {
        free(edgeValue);
        //delete[] edgeValue;
        edgeValue = 0x0;
    }

    if (nOv > 0) {

        edgeValue = (float*) malloc(nOv * sizeof (float));

        for (size_t i = 0; i < nOv; i++)
            edgeValue[i] = 0;
    }
}

void Node::allocateEdgeFlag() {
    if (edgeBitsFlag != 0x0) {
        free(edgeBitsFlag);
        //  delete[] edgeBitsFlag;
        edgeBitsFlag = 0x0;
    }

    if (nOv > 0) {
        edgeBitsFlag = (unsigned char*) malloc(nOv * sizeof (unsigned char));

        for (size_t i = 0; i < nOv; i++)
            edgeBitsFlag[i] = false;
    }
}

void Node::reallocMemory() {
    if (nOv == 0) {
        freeMemory();
        return;
    }
    void * p;

    p = realloc((void*) ovId, nOv * sizeof (unsigned int));
    ovId = (unsigned int*) p;

    p = realloc((void*) ovSize, nOv * sizeof (short));
    ovSize = (short*) p;

    if (edgeValue != 0x0) {
        p = realloc((void*) edgeValue, nOv * sizeof (float));
        edgeValue = (float*) p;
    }

    if (edgeBitsFlag != 0x0) {
        p = realloc((void*) edgeBitsFlag, nOv * sizeof (unsigned char));
        edgeBitsFlag = (unsigned char*) p;
    }
}

void Node::freeMemory() {
    if (ovId != 0x0) {
        free(ovId);
        //  delete[] ovId;
        ovId = 0x0;
    }
    if (ovSize != 0x0) {
        free(ovSize);
        //delete[] ovSize;
        ovSize = 0x0;
    }
    if (edgeValue != 0x0) {
        free(edgeValue);
        //delete[] edgeValue;
        edgeValue = 0x0;
    }
    if (edgeBitsFlag != 0x0) {
        free(edgeBitsFlag);
        //delete[] edgeBitsFlag;
        edgeBitsFlag = 0x0;
    }

    nOv = 0;
    nOvRight = 0;
    //  layout=0;
}

void Node::setEdgeSortedFlag(bool v) {
    edgeSorted = v;
}

void Node::computeOverlaps() {
    static vector<unsigned int>ID;
    static vector<short> SIZE;
    multiset<_OV, orderSet> ovSet;

    nOvRight = R->determineOverlaps2(getThisId(), ID, SIZE, ovSet);

    nOv = ID.size();
    allocate(nOv);

    for (unsigned int ii = 0; ii < nOv; ii++) {
        ovId[ii] = ID[ii];
        ovSize[ii] = SIZE[ii];
    }

    setOvComputed();

    // OverlapsGraph::nodeQueue.push_front(getThisId());

    OverlapsGraph::g_count++;
    if (OverlapsGraph::g_count % 1000 == 0) {
        cout << "Overlapped: " << OverlapsGraph::g_count++ << " Queued for reduction: " << OverlapsGraph::nodeQueue.size() << " "
                << 100 * (double) OverlapsGraph::nodeQueue.size() / R->getN_nrReads() << "%\r" << flush;
    }
}

unsigned int Node::getSequenceLength() {
    return L->getSequenceLength(layout);
};

double Node::getCoverage() {

    double nBases = L->getNReads(layout);
    unsigned int rL = R->getReadsLength();

    nBases *= rL; //n nucl.

    return nBases / (getSequenceLength() - rL + 1);
}

double Node::getUnbiaisedCoverage(double &nodeLength) {
    unsigned int nBases = L->getNReads(layout);

    // if (nBases < 3) //(==nReads at this point)
    //     return 0;

    nBases *= R->getReadsLength();
    double meanOv = L->getMeanOverlapSize(getLayout());

    if (nodeLength == R->getReadsLength()) {//empirical: consider adjacent overlaps
        double meanAdjOv = 0.0;

        for (unsigned int i = 0; i < getNOverlap(); i++)
            meanAdjOv += getOverlapSize(i);

        if (getNOverlap() != 0) {
            if (hasSingleSideOV())//
            {
                meanAdjOv /= (getNOverlap() + 1);
            } else {
                meanAdjOv /= getNOverlap();
            }
        }

        nodeLength -= meanAdjOv;
    } else
        nodeLength -= meanOv;

    return (double) nBases / nodeLength;
}

void Node::getVCoverage(vector<unsigned int> &cov, bool direction) {
    L->getVcoverage(getLayout(), direction, cov);
}

unsigned int Node::getNReads() {
    return L->getNReads(layout);
}

string Node::getDirectSequence() {
    return L->getSequence(layout);
}

string Node::getReverseSequence() {
    return L->getReverseSequence(layout);
}

void Node::getNeighbor(bool direction, size_t index, NodeIt& nIt) {
    if (!direction)
        index += nOvRight;

    nIt.pNode = N + (ovId[index]);
    nIt.absArrivalEdgeIndex = index;

    int s = (int) ovSize[index];
    if (s > 0) {
        nIt.arrivalEdgeSize = s;
        nIt.dir = direction;
    } else {
        nIt.arrivalEdgeSize = -s;
        nIt.dir = !direction;
    }

    nIt.initIterator();
}

unsigned int Node::getRightOverlap(unsigned int index, unsigned int &id, unsigned int &size, bool &direction) {
    if (index < 0 || index >= getNRightOverlap()) {
        id = 0;
        size = 0;
    } else {
        id = ovId[index];
        int s = (int) ovSize[index];
        if (s > 0) {
            direction = true;
            size = s;
        } else {
            size = -s;
            direction = false;
        }
    }
    return size;
}

unsigned int Node::getLeftOverlap(unsigned int index, unsigned int &id, unsigned int &size, bool &direction) {
    if (index < 0 || index >= getNLeftOverlap()) {
        id = 0;
        size = 0;
    } else {
        index += nOvRight;

        id = ovId[index];
        int s = (int) ovSize[index];
        if (s > 0) {
            direction = true;
            size = s;
        } else {
            size = -s;
            direction = false;
        }
    }
    return size;
}

unsigned int Node::getOverlap(unsigned int index, unsigned int &id, unsigned int&size, bool &direction) {
    if (index < 0 || index >= getNOverlap()) {
        id = 0;
        size = 0;
    } else {
        id = ovId[index];
        int s = (int) ovSize[index];
        if (s > 0) {
            direction = true;
            size = s;
        } else {
            size = -s;
            direction = false;
        }
    }
    return size;
}

unsigned int Node::getRightOverlapSize(unsigned int index) {
    if (index < 0)
        return 0;

    if (index > nOvRight - 1)
        cout << "Node::getRightOverlapSize() Problem" << endl;

    int s = (int) ovSize[index];
    if (s < 0)
        s = -s;

    return s;
}

unsigned int Node::getLeftOverlapSize(unsigned int index) {
    if (index > nOv - 1)
        cout << "Node::getLeftOverlapSize() Problem" << endl;

    index += nOvRight;
    int s = (int) ovSize[index];
    if (s < 0)
        s = -s;

    return s;
}

unsigned int Node::getOverlapSize(unsigned int index) {
    if (index > nOv - 1)
        cout << "Node::getOverlapSize() Problem" << endl;

    int s = (int) ovSize[index];
    if (s < 0)
        s = -s;

    return s;
}

bool Node::getOverlapDirection(unsigned int index) {
    return ovSize[index] > 0;
}

unsigned int Node::getOverlapId(unsigned int index) {
    return ovId[index];
};

unsigned int Node::getLeftOverlapId(unsigned int index) {
    if (index < nOv - nOvRight)
        return ovId[index + nOvRight];
    else
        return 0;
}

unsigned int Node::getRightOverlapId(unsigned int index) {
    if (index < nOvRight)
        return ovId[index];
    else
        return 0;
}

unsigned int Node::getNext(bool right, unsigned int index, unsigned int& id, unsigned int& size, bool& direction) {
    if (right) {
        return getRightOverlap(index, id, size, direction);
    } else {
        return getRightOverlap(index, id, size, direction);
    }
}

unsigned int Node::getEdgeIndex(bool right, unsigned int targetId, unsigned int targetSize, bool targetDirection) {
    //binary search on the edge tab, firstly based on the overlap size, then on the ID
    //edge size and IDs are sorted from the larger to the smaller

    if (edgeSorted == false) {
        cerr << "\nNode::getEdgeIndex(...) problem\n";
        sendBugReportPlease(cerr);
    }

    int high, low, i;

    if (right) {
        low = -1;
        high = nOvRight - 1;
    } else {
        low = nOvRight - 1;
        high = nOv - 1;
    }

    while (high - low > 1) {
        i = low + (high - low) / 2;

        if (targetSize > getOverlapSize(i))
            high = i;
        else if (targetSize < getOverlapSize(i))
            low = i;
        else //==
        {
            if (targetId >= getOverlapId(i))
                high = i;
            else
                low = i;
        }
    }

    if (high == -1)
        return nOv;
    if (getOverlapSize(high) != targetSize)
        return nOv;

    if (right) {
        while (getOverlapId(high) == targetId && (unsigned int) high < nOvRight) {
            if (getOverlapDirection(high) == targetDirection)
                return high;
            high++;
        }
    } else {
        while (getOverlapId(high) == targetId && (unsigned int) high < nOv) {
            if (getOverlapDirection(high) == targetDirection)
                return high;
            high++;
        }
    }

    return nOv; //not found
}

unsigned int Node::getReciprocal(unsigned int index) {
    unsigned int id, size;
    bool dir = true, right = false;

    if (index < nOvRight)
        right = true;

    getOverlap(index, id, size, dir);

    return N[id].getEdgeIndex(right != dir, getThisId(), size, dir);
}

void Node::removeEdge(unsigned int index) {
    //get reciprocal edge
    Node *pN = &N[ovId[index]];
    unsigned int rec = getReciprocal(index);

    if (pN == this) {//special case, self overlap
        if (rec != index) {//not palindromic => two edges are to be removed
            if (rec > index) {
                removeEdgeNR(rec);
                removeEdgeNR(index);
            } else {
                removeEdgeNR(index);
                removeEdgeNR(rec);
            }
            //            edgeBitsFlag[rec]=true;
            //            edgeBitsFlag[index]=true;
            //            removeMarkedEdge();
            return;
        }
    }

    if (rec == N[ovId[index]].nOv) {
        cerr << "Node::removeEdge problem\n";
        sendBugReportPlease(cerr);
    }

    for (unsigned int i = index; i < nOv - 1; i++) {
        ovId[i] = ovId[i + 1];
        ovSize[i] = ovSize[i + 1];
        if (edgeValue != 0x0)
            edgeValue[i] = edgeValue[i + 1];
        if (edgeBitsFlag != 0x0)
            edgeBitsFlag[i] = edgeBitsFlag[i + 1];
    }

    if (index < nOvRight)
        nOvRight--;
    nOv--;

    if (pN != this)//not a self palindromic overlap
    {
        for (unsigned int i = rec; i < pN->nOv - 1; i++) {
            pN->ovId[i] = pN->ovId[i + 1];
            pN->ovSize[i] = pN->ovSize[i + 1];
            if (pN->edgeValue != 0x0)
                pN->edgeValue[i] = pN->edgeValue[i + 1];
            if (pN->edgeBitsFlag != 0x0)
                pN->edgeBitsFlag[i] = pN->edgeBitsFlag[i + 1];
        }

        if (rec < pN->nOvRight)
            pN->nOvRight--;
        pN->nOv--;

    }
}

void Node::removeEdgeNR(unsigned int index)//do not care of the reciprocal
{
    for (unsigned int i = index; i < nOv - 1; i++) {
        ovId[i] = ovId[i + 1];
        ovSize[i] = ovSize[i + 1];
        if (edgeValue != 0x0)
            edgeValue[i] = edgeValue[i + 1];
        if (edgeBitsFlag != 0x0)
            edgeBitsFlag[i] = edgeBitsFlag[i + 1];
    }

    if (index < nOvRight)
        nOvRight--;
    nOv--;
}

void Node::removeEdges() {
    while (nOv > 0) {
        removeEdge(nOv - 1);
    }
}

void Node::removeAllEdges() {
    //do not take care of the reciprocal edges!
    nOv = nOvRight = 0;
}

void Node::flagReads(char state) {
    L->flagReads(layout, state);
}

void Node::initializeEdgeValues(float v) {
    for (size_t i = 0; i < getNOverlap(); i++)
        edgeValue[i] = v;
}

void Node::incrementEdgeValue(size_t i) {
    edgeValue[i]++;

    Node *pN = &N[ovId[i]];
    int rec = getReciprocal(i);
    pN->edgeValue[rec]++;
}

void Node::decrementEdgeValue(size_t i) {
    edgeValue[i]--;

    Node *pN = &N[ovId[i]];
    int rec = getReciprocal(i);
    pN->edgeValue[rec]--;
}

float Node::getEdgeValue(size_t i) {
    return edgeValue[i];
}

void Node::setEdgeValue(unsigned int index, float value) {
    Node * pN = &N[ovId[index]];
    int rec = getReciprocal(index);
    edgeValue[index] = value;
    pN->edgeValue[rec] = value;
}

void Node::setEdgeFlag(unsigned int index, bool state) {
    //mark edge i and its reciprocal
    Node * pN = &N[ovId[index]];
    unsigned int rec = getReciprocal(index);
    if (rec == pN->getNOverlap())
        cout << " Node::setEdgeFlag(unsigned int index, bool state) problem" << endl;
    edgeBitsFlag[index] = state;
    pN->edgeBitsFlag[rec] = state;

}

void Node::setEdgeFlagUnrec(unsigned int index, bool state) {
    //do not take care of the reciprocal
    edgeBitsFlag[index] = state;
}

void Node::initEdgeFlags(bool state) {
    for (unsigned int i = 0; i < nOv; i++) {
        edgeBitsFlag[i] = state;
        // setEdgeFlag(i,state);
    }
}

int Node::removeMarkedEdge() {
    //!does not care of the reciprocal edge!
    //!Must be called on all nodes

    int R = nOvRight;
    int N = nOv;
    int p = 0;
    for (unsigned int i = 0; i < nOv; i++) {
        if (getEdgeFlag(i) == false) //=edge to keep
        {
            ovId[p] = ovId[i];
            ovSize[p] = ovSize[i];
            edgeBitsFlag[p] = edgeBitsFlag[i];

            //not necessarily allocated
            //not allocated during the overlapping mode
            if (edgeValue != 0x0)
                edgeValue[p] = edgeValue[i];
            p++;
        } else {
            if (i < nOvRight)
                R--;
            N--;
        }
    }

    int nRemoved = nOv - N;
    nOvRight = R;
    nOv = N;
    return nRemoved;
}

int Node::removeEdgesByValue(float suspectCutoff, float GIncoherentCutoff) {
    //!does not care of the reciprocal edge!
    //!must be called on all nodes

    int R = nOvRight;
    int N = nOv;
    int p = 0;

    //    string seq;
    //    path.clear();
    //    path.push_back(0);
    //    path.push_back(getThisId());
    //    path.push_back(0);
    //    unsigned int rl=this->R->getReadsLength();
    string right, left, seq1, seq2;
    unsigned int endRead = L->getEnd(getLayout());
    if (L->getDirection(endRead))
        right = L->getDirectRead(endRead);
    else
        right = L->getReverseRead(endRead);

    endRead = L->getBegin(getLayout());
    if (L->getDirection(endRead))
        left = L->getReverseRead(endRead);
    else
        left = L->getDirectRead(endRead);

    NodeIt nIt, nIt2;

    for (unsigned int i = 0; i < nOv; i++) {
        //        outGlob << ">" << getThisId() << "->" << this->ovId[i] << "_" << getOverlapSize(i) << '_';
        //        if ((getEdgeValue(i) < suspectCutoff && getEdgeFlag(i) == false) || getEdgeValue(i) < GIncoherentCutoff)
        //            outGlob << "Gincoherent\n";
        //        else
        //            outGlob << "OK\n";
        //        if (i < nOvRight) {
        //            nIt.initNodeIt(getThisId(), true);
        //            seq1 = right;
        //        } else {
        //            nIt.initNodeIt(getThisId(), false);
        //            seq1 = left;
        //        }
        //
        //        nIt.getNeighborGivenAbsoluteIndex(i, nIt2);
        //        nIt2.getSequence(seq2);
        //        seq2 = seq2.substr(this->getOverlapSize(i));
        //        seq1 += seq2[0];

        //outGlob << seq1 << '\n';


        // if ( (getEdgeValue(i) < suspectCutoff && getEdgeFlag(i)==false) || (getEdgeValue(i) < GIncoherentCutoff && getEdgeFlag(i)==false))
        if ((getEdgeValue(i) < suspectCutoff
                && getEdgeFlag(i) == false)
                || getEdgeValue(i) < GIncoherentCutoff) {//false positive edge

            //            outGlob << ">node" << getThisId() << "to" << this->getOverlapId(i) << " " << this->getOverlapSize(i) << "\n";
            //            path[2]=i;
            //            seq=G->pathToSeq(path);
            //            seq=seq.substr(this->getSequenceLength()-getOverlapSize(i)-1);
            //            seq = seq.substr(0,rl+1);
            //            outGlob << seq << '\n';

            if (i < nOvRight)
                R--;
            N--;
        } else {
            ovId[p] = ovId[i];
            ovSize[p] = ovSize[i];
            // edgeBitsFlag[p]=edgeBitsFlag[i];
            edgeBitsFlag[p] = false; //set back to default value
            edgeValue[p] = edgeValue[i];
            p++;
        }
    }

    int nRemoved = nOv - N;
    nOvRight = R;
    nOv = N;
    return nRemoved;
}

//compute edges prob by sampling maxD bp around

void Node::computeEdgesProb(unsigned int maxD, unsigned int maxN, double reliableCutoff) {
    //reliableCutoff: used to flag edges that are, or have at least one reliable brother

    unsigned int distr1[512];
    unsigned int distr2[512];
    unsigned int nOhSampled1=0;
    unsigned int nOhSampled2=0;
    unsigned int oh;
    unsigned int totSample1 = 0, sum1 = 0, totSample2 = 0, sum2 = 0;
    double prob, mean;
    unsigned int rl = R->getReadsLength();
    bool done = false;
    NodeIt nIt, nextIt;

    bool sampleDir = true;

    for (int i = 0; i < 2; i++) {

        nIt.initNodeIt(getThisId(), sampleDir);
        bool atLeastOneReliable = false; //means at least one reliable on each side of the edge
        // a G-incoherent edge should be accompanied by G-coherent edges

        while (nIt.getNext(nextIt) != 0) {
            if (getEdgeValue(nextIt.getAbsArrivalEdgeIndex()) >= reliableCutoff)
                atLeastOneReliable = true;

            if (getEdgeValue(nextIt.getAbsArrivalEdgeIndex()) != -1.0)
                continue; //edge already examined

            oh = rl - nextIt.getArrivalEdgeSize();
            if (oh == 1) //P(oh=1) set to 1.0 by definition
            {
                setEdgeValue(nextIt.getAbsArrivalEdgeIndex(), 1.0);
                atLeastOneReliable = true;
                continue;
            }

            if (!done) {
                //sample this
                memset(distr1, 0, rl * sizeof (unsigned int));
                nOhSampled1 = sampleOverHangs(!sampleDir, maxD, maxN, distr1);
                totSample1 = sum1 = 0;
                for (unsigned int i = 0; i < rl; i++) {
                    totSample1 += distr1[i];
                    sum1 += i * distr1[i];
                }

                done = true;
            }
            memset(distr2, 0, rl * sizeof (unsigned int));
            nOhSampled2 = nextIt.pNode->sampleOverHangs(nextIt.dir, maxD, maxN, distr2);

            totSample2 = sum2 = 0;
            for (unsigned int i = 0; i < rl; i++) {
                totSample2 += distr2[i];
                sum2 += i * distr2[i];
            }

            totSample2 += totSample1;
            sum2 += sum1;
            nOhSampled2 += nOhSampled1;
            oh = rl - nextIt.getArrivalEdgeSize();

            if (totSample2 < 25) {
                setEdgeValue(nextIt.getAbsArrivalEdgeIndex(), 2.0); //canceled  
            } else {
                mean = (double) sum2 / totSample2;
                prob = cdfOH(1.0 / (mean + 1), oh, (double) totSample2 / nOhSampled2);
                //  prob = cdfOH(1.0 / (mean + 1), oh, 1.0);
                if (prob >= reliableCutoff)
                    atLeastOneReliable = true;
                setEdgeValue(nextIt.getAbsArrivalEdgeIndex(), prob);
            }
        }

        if (atLeastOneReliable == false) {
            nIt.initNodeIt(getThisId(), sampleDir);
            while (nIt.getNext(nextIt) != 0) {
                setEdgeFlag(nextIt.getAbsArrivalEdgeIndex(), true);
            }
        }
        sampleDir = !sampleDir;
    }
}

unsigned int Node::bubble(bool dir, double minCov) {
    NodeIt nIt, bub1, bub2, end1, end2;
    unsigned int maxLength = 2 * R->getReadsLength() - 1;
    string s1, s2;
    double cov1, cov2, ratio;

    nIt.initNodeIt(getThisId(), dir);

    if (nIt.getNNeighbor() == 2) {
        nIt.getNext(bub1);
        nIt.getNext(bub2);
        if (bub1.getNNeighbor() == 1 &&
                bub2.getNNeighbor() == 1 &&
                bub1.getNodeLength() <= maxLength &&
                bub2.getNodeLength() <= maxLength) {
            bub1.getNext(end1);
            bub2.getNext(end2);

            if (end1 == end2) {//potential bubble
                // dotLocalGraph(5, "debug");
                bub1.getSequence(s1);
                bub2.getSequence(s2);
                s1 = s1.substr(bub1.getArrivalEdgeSize());
                s2 = s2.substr(bub2.getArrivalEdgeSize());
                s1 = s1.substr(0, s1.length() - end1.getArrivalEdgeSize());
                s2 = s2.substr(0, s2.length() - end2.getArrivalEdgeSize());

                //                int dSize=s1.size()-s2.size();
                //                if (dSize<0)
                //                    dSize=-dSize;

                if (//(dSize <=1 && s1.size() <=1) || 
                        hamming(s1.c_str(), s2.c_str()) == 1) {
                    cov1 = bub1.getCoverage();
                    cov2 = bub2.getCoverage();

                    if (cov1 > cov2)
                        ratio = cov2 / cov1;
                    else
                        ratio = cov1 / cov2;

                    if (ratio <= 0.4 && (cov1 <= minCov || cov2 <= minCov)) {
                        //cout << "bub:" << getThisId() << endl;
                        if (cov1 > cov2)
                            bub2.discard();
                        else
                            bub1.discard();

                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

unsigned int Node::sampleOverHangs(bool direction,
        unsigned int maxD,
        unsigned int maxN,
        unsigned int *distr) {
    nodeList.clear();

    unsigned int v = sampleOverHangsRec(direction, maxD, maxN, 0, 0, distr);

    for (size_t i = 0; i < nodeList.size(); i++)
        N[nodeList.at(i)].unsetVisited();

    return v;
}

unsigned int Node::sampleOverHangsRec(bool direction,
        unsigned int maxD,
        unsigned int maxN,
        unsigned int currentD,
        unsigned int currentN,
        unsigned int *distr) {
    unsigned int nOH;
    unsigned int nodeLength = L->getSequenceLength(getLayout());

    nodeList.push_back(getThisId()); //memorize marked nodes
    setVisited();

    nOH = L->sampleOH(getLayout(), direction, maxD - currentD, maxN - currentN, distr);

    currentD += nodeLength - R->getReadsLength() + 1; //only if node fully sampled, but not an issue
    currentN += nOH;

    if (currentD < maxD && currentN < maxN) {//sample further

        NodeIt nIt, nextIt;
        Node* nodeMem = 0x0;
        nIt.initNodeIt(getThisId(), direction);

        double cov, maxCov = 0.0;
        while (nIt.getNext(nextIt)) { //chose the one to go one
            // -> avoid recursion on all neighbors

            if (nextIt.pNode->isVisited())
                continue;

            cov = nextIt.getCoverage();
            if (cov > maxCov) {
                maxCov = cov;
                nodeMem = nextIt.pNode;
            }
        }

        if (nodeMem != 0x0) {
            nodeList.push_back(nodeMem->getThisId());
            nodeMem->setVisited();
            nOH += nodeMem->sampleOverHangsRec(nextIt.dir,
                    maxD,
                    maxN,
                    // currentD-nextIt.getArrivalEdgeSize(),
                    currentD,
                    currentN,
                    distr);
        }

        //        //recursion loop
        //        while (nIt.getNext(nextIt))
        //        {
        //            if (nextIt.pNode->isTraversed())
        //                continue;
        //            else
        //            {
        //                nodeList.push_back(nextIt.getNodeId());
        //                nextIt.pNode->setTraversed();
        //            }
        //
        //            //Sampling interNode overHanging
        //            //may introduce spurious oh
        //            //oh=R->getReadsLength()-nextIt.getArrivalEdgeSize();
        //            //distr[oh]++;
        //            //nOH++;
        //            //currentN++;
        //            //currentD+=oh;
        //            
        //            
        //            nOH+=nextIt.pNode->sampleOverHangsRec(nextIt.dir,
        //                                                  maxD,
        //                                                  maxN,
        //                                                 // currentD-nextIt.getArrivalEdgeSize(),
        //                                                  currentD,
        //                                                  currentN,
        //                                                  distr);
        //        }
    }


    return nOH;
}

bool Node::testReciprocal() {
    unsigned int index;
    bool passed = true;

    for (unsigned int i = 0; i < nOv; i++) {
        index = getReciprocal(i);
        if (index == N[ovId[i]].getNOverlap()) {
            cout << "node: " << getThisId() << " edgeIndex:" << i << " to node: " << ovId[i] << " not reciprocal" << endl;
            cout << "size: " << getOverlapSize(i);
            if (i < nOvRight)
                cout << "-->";
            else
                cout << "<--";
            if (this->getOverlapDirection(i))
                cout << " +" << endl;
            else
                cout << " -" << endl;

            passed = false;
        } else {
            if (edgeValue != 0x0) {
                if (edgeValue[i] != N[ovId[i]].getEdgeValue(index)) {
                    cout << "node: " << getThisId() << " edgeIndex:" << i << " to node: " << ovId[i] << " edge value problem" << endl;
                    cout << "size: " << this->getOverlapSize(i);
                    if (i < nOvRight)
                        cout << "-->";
                    else
                        cout << "<--";
                    if (this->getOverlapDirection(i))
                        cout << " +" << endl;
                    else
                        cout << " -" << endl;
                    passed = false;
                }
            }
            if (edgeBitsFlag[i] != N[ovId[i]].getEdgeFlag(index)) {
                cout << "node: " << getThisId() << " edgeIndex:" << i << " to node: " << ovId[i] << " edge flag problem" << endl;
                cout << "size: " << this->getOverlapSize(i);
                if (i < nOvRight)
                    cout << "-->";
                else
                    cout << "<--";
                if (this->getOverlapDirection(i))
                    cout << " +" << endl;
                else
                    cout << " -" << endl;
                passed = false;
            }
        }
    }

    return passed;
}

//an other slow method to check the edges reciprocity
//but that do not require the edges to be sorted.

bool Node::testReciprocal2() {
    unsigned int id = 0, id2 = 0, size = 0, size2 = 0;
    bool direction = true, direction2 = true;
    bool ok;

    for (unsigned int i = 0; i < getNRightOverlap(); i++) {
        getRightOverlap(i, id, size, direction);
        ok = false;
        if (direction) {
            for (unsigned int j = 0; j < N[id].getNLeftOverlap(); j++) {
                N[id].getLeftOverlap(j, id2, size2, direction2);
                if (id2 == getThisId() && size2 == size && direction == direction2) {
                    ok = true;
                    break;
                }
            }
        } else {
            for (unsigned int j = 0; j < N[id].getNRightOverlap(); j++) {
                N[id].getRightOverlap(j, id2, size2, direction2);
                if (id2 == getThisId() && size2 == size && direction == direction2) {
                    ok = true;
                    break;
                }
            }
        }
        //
        //        if (getThisId() == id)
        //        {
        //            cout << "self " << id << endl;
        //        }

        if (ok == false) {
            cout << "read " << getThisId() << " ovRight " << id << " not reciproqual " << endl;
            cout << id2 << " " << size2 << " " << direction2 << '\n';
            cout << id << " " << size << " " << direction << '\n';

            return false;
        }

    }

    for (unsigned int i = 0; i < getNLeftOverlap(); i++) {
        getLeftOverlap(i, id, size, direction);
        ok = false;
        if (!direction) {
            for (unsigned int j = 0; j < N[id].getNLeftOverlap(); j++) {
                N[id].getLeftOverlap(j, id2, size2, direction2);
                if (id2 == getThisId() && size2 == size && direction == direction2) {
                    ok = true;
                    break;
                }
            }
        } else {
            for (unsigned int j = 0; j < N[id].getNRightOverlap(); j++) {
                N[id].getRightOverlap(j, id2, size2, direction2);
                if (id2 == getThisId() && size2 == size && direction == direction2) {
                    ok = true;
                    break;
                }

            }
        }

        //        if (getThisId() == id)
        //        {
        //            cout << "self " << id << endl;
        //        }

        if (ok == false) {
            cout << "read " << getThisId() << " ovLeft " << id << " not reciproqual " << endl;
            cout << id2 << " " << size2 << " " << direction2 << '\n';
            cout << id << " " << size << " " << direction << '\n';

            return false;
        }


    }

    return true;
}

void Node::overlapSizeCutoff(unsigned int cutoff) {
    //must be called on for ALL nodes

    for (size_t i = 0; i < nOv; i++)
        edgeBitsFlag[i] = false;

    //right
    for (unsigned int i = 0; i < getNRightOverlap(); i++) {
        if (getRightOverlapSize(i) < cutoff) {
            edgeBitsFlag[i] = true;
        }
    }

    //left
    for (unsigned int i = 0; i < getNLeftOverlap(); i++) {
        if (getLeftOverlapSize(i) < cutoff) {
            edgeBitsFlag[i + nOvRight] = true;
        }
    }

    removeMarkedEdge();
}

void Node::removeShortBranchingOverlaps(unsigned int cutoff) {
    unsigned int size;

    if (isBranching(true)) {
        for (unsigned int i = 1; i <= getNRightOverlap(); i++) {
            size = getOverlapSize(i - 1);
            if (size < cutoff) {
                removeEdge(i - 1);
                i--;
            }
        }
    }
    if (isBranching(false)) {
        for (unsigned int i = 1; i <= getNLeftOverlap(); i++) {
            size = getOverlapSize(i - 1);
            if (size < cutoff) {
                removeEdge((i - 1) + nOvRight);
                i--;
            }
        }
    }
}

void Node::markTransitiveEdges() {

    unsigned int id, idB, size, sizeB;
    bool direction, directionB;
    unsigned int longestOH;
    unsigned int overhang;
    unsigned int nEdges;

    int targetSize;
    unsigned int edgeIndex;

    //Linear expected algorithm adapted from:
    //E.W. Myers, The fragment assembly graph, Bioinformatics vol 21 Suppl.2 2005 pages ii79-ii85

    for (unsigned int i = 0; i < getNRightOverlap(); i++) {
        getNeighbor(true, i, id, size, direction);

        if (N[id].isVisited(direction))
            N[id].setMultipleEdges(direction); //indicate multiple incoming edges
        else
            N[id].setVisited(direction);
    }

    longestOH = R->getReadsLength() - size; //longest overhanging

    for (unsigned int i = 0; i < getNRightOverlap(); i++) {
        if (edgeBitsFlag[i] == true) {
            continue;
        }
        getNeighbor(true, i, id, size, direction);

        if (N[id].isEliminated(direction) == false) {
            overhang = R->getReadsLength() - size;

            nEdges = N[id].getNEdges(direction);

            for (unsigned int ii = 0; ii < nEdges; ii++) {
                global_count++;
                N[id].getNeighbor(direction, ii, idB, sizeB, directionB);
                if ((R->getReadsLength() - sizeB) + overhang > longestOH)
                    break;

                if (N[idB].hasMultipleEdges(directionB)) {//in this case, the "good" transitive edge from *this has to be searched
                    //given current distance and orientation
                    //this point as been added from Myer's algorithm which do not support
                    //multiple edges between two nodes

                    targetSize = size + sizeB - R->getReadsLength();
                    edgeIndex = getEdgeIndex(true, idB, targetSize, directionB);
                    if (edgeIndex == getNOverlap())
                        cout << "";

                    edgeBitsFlag[edgeIndex] = true;
                } else if (N[idB].isVisited(directionB))
                    N[idB].setEliminated(directionB);
            }
        }
    }

    for (unsigned int i = 0; i < getNRightOverlap(); i++) {
        getNeighbor(true, i, id, size, direction);

        if (N[id].isEliminated(direction))
            edgeBitsFlag[i] = true;
    }

    for (unsigned int i = 0; i < getNRightOverlap(); i++) {
        getNeighbor(true, i, id, size, direction);
        //  N[id].flagsBit&=0;
        N[id].flagsBit &= ~254;
    }


    for (unsigned int i = 0; i < getNLeftOverlap(); i++) {
        getNeighbor(false, i, id, size, direction);

        if (N[id].isVisited(direction))
            N[id].setMultipleEdges(direction); //indicate multiple incoming edges
        else
            N[id].setVisited(direction);
    }

    longestOH = R->getReadsLength() - size; //longest overhanging

    for (unsigned int i = 0; i < getNLeftOverlap(); i++) {
        if (edgeBitsFlag[i + nOvRight] == true)
            continue;

        getNeighbor(false, i, id, size, direction);

        if (N[id].isEliminated(direction) == false) {
            overhang = R->getReadsLength() - size;

            nEdges = N[id].getNEdges(direction);

            for (unsigned int ii = 0; ii < nEdges; ii++) {
                global_count++;
                N[id].getNeighbor(direction, ii, idB, sizeB, directionB);
                if ((R->getReadsLength() - sizeB) + overhang > longestOH)
                    break;

                if (N[idB].hasMultipleEdges(directionB)) {//in this case, the "good" transitive edge from *this has to be searched
                    //given current distance and orientation
                    //this point as been added from Myer's algoithm which do not take account for
                    //multiple edges between two nodes

                    targetSize = size + sizeB - R->getReadsLength();
                    edgeIndex = getEdgeIndex(false, idB, targetSize, !directionB);
                    if (edgeIndex == getNOverlap())
                        cout << "";

                    edgeBitsFlag[edgeIndex] = true;
                } else if (N[idB].isVisited(directionB))
                    N[idB].setEliminated(directionB);
            }
        }
    }

    for (unsigned int i = 0; i < getNLeftOverlap(); i++) {
        getNeighbor(false, i, id, size, direction);

        if (N[id].isEliminated(direction))
            edgeBitsFlag[i + nOvRight] = true;
    }

    for (unsigned int i = 0; i < getNLeftOverlap(); i++) {
        getNeighbor(false, i, id, size, direction);
        // N[id].flagsBit&=0;
        N[id].flagsBit &= ~254;
    }
}


//void Node::markTransitiveEdges()
//{ //naive algorithm
//    //mark transitive edges
//    unsigned int i,j;
//
//    //right overlaps
//    i=0;
//
//    while ( i+1 < getNRightOverlap() )
//    {
//        if (getEdgeFlag(i)==false)
//        {
//            j=i+1;
//
//            while (j < getNRightOverlap() )
//            {
//                if (getEdgeFlag(i)==false && getOverlapSize(i) != getOverlapSize(j))
//                {
//                    if (areConsistentR(i,j))
//                    {
//                       // cout << "node: " << getThisId() << " edge: " << j << " flagged" << endl;
//                        edgeBitsFlag[j]=true;
//                    }
//                }
//                j++;
//            }
//        }
//        i++;
//    }
//
//    //left overlaps
//    i=0;
//
//    while ( i+1 < getNLeftOverlap() )
//    {
//        if (getEdgeFlag(i+nOvRight)==false)
//        {
//            j=i+1;
//            while (j < getNLeftOverlap() )
//            {
//                if (getEdgeFlag(j+nOvRight)==false && getOverlapSize(i+nOvRight) != getOverlapSize(j+nOvRight))
//                {
//                    if (areConsistentL(i,j))
//                    {
//                       // cout << "node: " << getThisId() << " edge: " << j << " flagged" << endl;
//                        edgeBitsFlag[j+nOvRight]=true;
//                    }
//                }
//                j++;
//            }
//        }
//        i++;
//    }
//}

void Node::initDot(ostream &out) {
    for (unsigned int i = 0; i <= nNodes; i++) {
        N[i].unsetVisited();
        N[i].unsetTargeted();
        N[i].initEdgeFlags(false);
    }

    Node::dotNodeList.clear();
    out << "digraph overlapGraph {\n";
    out << "node [shape = ellipse,fontsize = 7];\n";
    out << "edge [dir=both];\n";
}

void Node::dotAround(ostream &out, int depthLimit, int currentDepth) {
    unsigned int id, size;
    bool direction=true;

    //avoid cycling  
    if (isVisited() == true) {
        return;
    }

    Node::dotNodeList.push_back(getThisId());

    if (currentDepth == depthLimit)
        return;

    setVisited();
    currentDepth++;
    if (isTargeted())
        currentDepth = 0;

    for (size_t i = 0; i < this->getNOverlap(); i++) {
        getOverlap(i, id, size, direction);

        //already drawn
        if (edgeBitsFlag[i] == true)
            continue;
        setEdgeFlag(i, true);

        out << getThisId() << " -> " << id;
        out << " [label = \"" << size << "\\n";
        out << setprecision(2) << getEdgeValue(i) << "\\n";

        //        if (getEdgeFlag(i))
        //            out << "T";
        //        else
        //            out << "F";
        out << "\", ";
        out << " arrowsize = 0.7,";

        if (isTargeted() && N[id].isTargeted())
            out << " penwidth = 2,";
        else
            out << " penwidth = 1,";

        if (i >= nOvRight) {
            direction = !direction;
            out << "arrowtail = normal,";
        } else
            out << "arrowtail = inv,";

        if (direction)
            out << "arrowhead = normal];\n";
        else
            out << "arrowhead = inv];\n";




    }

    for (size_t i = 0; i < this->getNOverlap(); i++) {
        getOverlap(i, id, size, direction);
        //  N[id].writeDotGraph(out, depthLimit, currentDepth);
        N[id].dotAround(out, depthLimit, currentDepth);
    }
}

void Node::closeDot(ostream& out) {
    for (vector<unsigned int>::iterator it = Node::dotNodeList.begin(); it != Node::dotNodeList.end(); it++) {
        bool missingEdges = false;
        for (unsigned int i = 0; i < N[*it].getNOverlap(); i++) {
            if (N[*it].edgeBitsFlag[i] == false) {
                missingEdges = true;
                break;
            }
        }

        out << *it << " [";
        out << "label = \""
                << *it << "\\nl="
                << N[*it].getSequenceLength()
                << " c=" << N[*it].getCoverage()
                // << " " << (int)N[*it].value 
                << "\"";
        if (missingEdges)
            out << " style = dashed";
        else
            out << " style = solid";
        if (N[*it].isTargeted())
            out << " ,style= filled, fillcolor = orange";
        //        else
        //            out << " ,fillcolor = white";
        out << "]\n";
    }

    for (vector<unsigned int>::iterator it = Node::dotNodeList.begin(); it != Node::dotNodeList.end(); it++) {
        N[*it].unsetVisited();
        N[*it].initEdgeFlags(false);
    }

    out << "}\n";
}

void Node::dotLocalGraph(int depth, string fileName) {

    for (unsigned int i = 0; i <= nNodes; i++) {
        N[i].unsetVisited();
        N[i].initEdgeFlags(false);
    }

    Node::dotNodeList.clear();

    ofstream out((fileName + ".dot").c_str());

    out << "digraph overlapGraph {\n";
    //out << "size=\"8,10\";" << endl;
    out << "node [shape = ellipse,fontsize = 7];\n";
    out << "edge [dir=both];\n";

    writeDotGraph(out, depth, 0);

    for (vector<unsigned int>::iterator it = Node::dotNodeList.begin(); it != Node::dotNodeList.end(); it++) {

        bool missingEdges = false;
        for (unsigned int i = 0; i < N[*it].getNOverlap(); i++) {
            if (N[*it].edgeBitsFlag[i] == false) {
                missingEdges = true;
                break;
            }
        }

        out << *it << " [";
        out << "label = \""
                << *it << "\\nl="
                << N[*it].getSequenceLength()
                << " c=" << N[*it].getCoverage()
                // << " " << (int)N[*it].value 
                << "\"";
        if (missingEdges)
            out << " style = dashed";
        out << "]\n";

    }

    for (vector<unsigned int>::iterator it = Node::dotNodeList.begin(); it != Node::dotNodeList.end(); it++) {
        N[*it].unsetVisited();
        N[*it].initEdgeFlags(false);
    }

    out << "}\n";
    out.close();
    out.clear();

    string command = "dot -Tps -o ";
    command += fileName + ".ps " + fileName + ".dot";
    cout << "executing " << command << " ... " << flush;
    int unused __attribute__((unused)); //get rid of gcc warning
    unused = system(command.c_str());
    cout << "done" << endl;

}

void Node::writeDotGraph(ostream &out, int &depthLimit, int currentDepth) {
    unsigned int id, size;
    bool direction=true;

    //avoid cycling  
    if (isVisited() == true) {
        return;
    }

    Node::dotNodeList.push_back(getThisId());

    if (currentDepth == depthLimit)
        return;

    setVisited();
    currentDepth++;

    for (size_t i = 0; i < this->getNOverlap(); i++) {
        getOverlap(i, id, size, direction);

        //already drawn
        if (edgeBitsFlag[i] == true)
            continue;
        setEdgeFlag(i, true);

        out << getThisId() << " -> " << id;
        out << " [label = \"" << size << "\\n";
        out << setprecision(2) << getEdgeValue(i) << "\\n";

        //        if (getEdgeFlag(i))
        //            out << "T";
        //        else
        //            out << "F";
        out << "\", ";
        out << " arrowsize = 0.7,";

        if (i >= nOvRight) {
            direction = !direction;
            out << "arrowtail = normal,";
        } else
            out << "arrowtail = inv,";

        if (direction)
            out << "arrowhead = normal];\n";
        else
            out << "arrowhead = inv];\n";
    }

    for (size_t i = 0; i < this->getNOverlap(); i++) {
        getOverlap(i, id, size, direction);
        N[id].writeDotGraph(out, depthLimit, currentDepth);
    }
}

unsigned int Node::getLongestNonAmbiguousPath(string &s, vector<unsigned int> &coverage, unsigned int &leftId, bool &leftDir, unsigned int &rightId, bool &rightDir) {
    nodeList.clear();
    s.clear();
    unsigned int id, size;
    bool direction = true;
    unsigned int currentId = getThisId();
    coverage.clear();

    //left
    bool state = false;
    N[currentId].setVisited(state);
    nodeList.push_back(currentId);

    do {
        if (!N[currentId].hasSingleEdge(state)) {
            break; // currentId is the left most read of the Contig
        }

        if (state)
            N[currentId].getRightOverlap(0, id, size, direction);
        else
            N[currentId].getLeftOverlap(0, id, size, direction);

        if (!N[id].isBranching(direction ? !state : state)) {
            //if no incoming path

            //since a node having a single self reverse edge is considered
            //as a branching node, the following two possibilities should be equivalent

            //if (N[id].isTraversed(direction ? state : !state))
            if (N[id].isVisited())
                break;

            //id is added to the Contig
            if (!direction)
                state = !state;

            N[id].setVisited(state);
            nodeList.push_back(id);
            currentId = id;
        } else
            break; //id is not included in the Contig

    } while (true);

    leftId = currentId;
    leftDir = !state;

    for (vector<unsigned int>::iterator it = nodeList.begin(); it != nodeList.end(); it++)
        N[*it].unsetVisited();

    nodeList.clear();
    state = !state; // U-turn

    unsigned int ContigLength = N[currentId].getSequenceLength();
    unsigned int nConcatenatedNode = 1;

    unsigned int myLayout;

    if (state)
        myLayout = N[currentId].getLayout();
    else
        myLayout = L->reverseComplement(N[currentId].getLayout());

    N[currentId].setLayout(0);
    S_counter++;

    N[currentId].setVisited(state);
    nodeList.push_back(currentId);

    do {
        if (!N[currentId].hasSingleEdge(state)) {
            break; //
        }

        if (state)
            N[currentId].getRightOverlap(0, id, size, direction);
        else
            N[currentId].getLeftOverlap(0, id, size, direction);

        if (!N[id].isBranching(direction ? !state : state)) {
            //if no incoming path


            //if (N[id].isTraversed(direction ? state : !state))
            if (N[id].isVisited())
                break;

            //id is added to the Contig

            if (!direction)
                state = !state;

            N[id].setVisited(state);
            nodeList.push_back(id);

            if (state) {
                myLayout = L->merge(myLayout, N[id].getLayout(), true, L->getSequenceLength(myLayout) - size);
            } else {
                myLayout = L->merge(myLayout, N[id].getLayout(), false, L->getSequenceLength(myLayout) - size);
            }
            N[id].setLayout(0);

            ContigLength += N[id].getSequenceLength() - size;
            nConcatenatedNode++;
            currentId = id;
            S_counter++;
        } else
            break;

    } while (true);

    rightId = currentId;
    rightDir = state;

    setLayout(myLayout);
    L->setLayoutNodeId(myLayout, getThisId());

    //return ContigLength;
    return nConcatenatedNode;
}

void Node::sortByIds() {
    unsigned int startIndex, endIndex, index;
    unsigned int size;

    if (getNRightOverlap() >= 2) {
        index = 0;
        startIndex = 0;
        size = getOverlapSize(0);
        index++;

        while (index < getNRightOverlap()) {
            if (getOverlapSize(index) != size) {
                endIndex = index - 1;
                sortByIds(startIndex, endIndex);
                startIndex = index;
                size = getOverlapSize(index);
            }
            index++;
        }
        endIndex = getNRightOverlap() - 1;
        sortByIds(startIndex, endIndex);
    }

    if (getNLeftOverlap() >= 2) {
        index = getNRightOverlap();
        startIndex = index;
        size = getOverlapSize(index);
        index++;

        while (index < getNOverlap()) {
            if (getOverlapSize(index) != size) {
                endIndex = index - 1;
                sortByIds(startIndex, endIndex);
                startIndex = index;
                size = getOverlapSize(index);
            }
            index++;
        }
        endIndex = getNOverlap() - 1;
        sortByIds(startIndex, endIndex);
    }
}

bool orderf(Edge i, Edge j) {
    return (i.ovId > j.ovId);
    //sort from the larger to the smaller
}

void Node::sortByIds(unsigned int startIndex, unsigned int endIndex) {

    if (startIndex < 0 || startIndex > nOv - 1 || endIndex < 0 || endIndex > nOv - 1 || startIndex > endIndex) {
        cout << "Node::sortByIds(...) problem\n";
        exit(0);
    }

    Edge myEdge;
    static vector<Edge> edgesTmp;
    edgesTmp.clear();

    //copy edges in Edge vector
    for (unsigned int i = startIndex; i <= endIndex; i++) {
        myEdge.ovId = ovId[i];
        myEdge.ovSize = ovSize[i];
        myEdge.edgeBitFlag = edgeBitsFlag[i];
        if (edgeValue != 0x0)
            myEdge.edgeValue = edgeValue[i];

        edgesTmp.push_back(myEdge);
    }

    sort(edgesTmp.begin(), edgesTmp.end(), orderf);
    vector<Edge>::iterator eIt = edgesTmp.begin();

    //copy back vector to individual edge fields
    for (unsigned int i = startIndex; i <= endIndex; i++) {

        ovId[i] = eIt->ovId;
        ovSize[i] = eIt->ovSize;
        edgeBitsFlag[i] = eIt->edgeBitFlag;
        if (edgeValue != 0x0)
            edgeValue[i] = eIt->edgeValue;
        eIt++;
    }
}

int Node::searchPathToPaired(bool currentDir,
        unsigned int targetNode,
        bool targetDir,
        int minDistance,
        int maxDistance,
        int shift,
        double maxLogNLeaves) {
    double cancel = 0;
    path.clear();
    //VP.clear();

    if (targetNode == 0)
        return 0;
    //   vPath.clear();
    path.push_back(0); //index 0 stores the length of the path
    path.push_back(getThisId());

    pathLengthMean = pathLengthSD = 0.0;
    nPathFound = 0;

    searchPathToPaired(currentDir, targetNode, targetDir, minDistance, maxDistance, shift, 0, cancel, maxLogNLeaves);

    if (cancel == -1) {
        // VP.clear();
        pathLengthMean = pathLengthSD = 0.0;
        path.clear();
        return -1; //canceled
    }

    if (nPathFound > 1) {
        pathLengthMean /= nPathFound;
        pathLengthSD = (nPathFound - 1) *
                (
                pathLengthSD -
                nPathFound * pathLengthMean * pathLengthMean
                );
    } else
        pathLengthSD = 0.0;

    // return VP.getNPath();
    return nPathFound;
}

//note: node1->searchPathToPaired(...,node2,...) could return a different result than node2->searchPathToPaired(...,node1,...)
//when one of the search is canceled

int Node::searchPathToPaired(bool currentDir,
        unsigned int targetNode,
        bool targetDir,
        int minDistance,
        int maxDistance,
        int shift,
        int currentDistance,
        double &cancel,
        double maxLogNLeaves) {
    unsigned int id, size;
    bool direction;
    unsigned int distance;
    int retVal = 0;

    if (currentDistance > maxDistance)
        return 0;

    if (cancel == -1.0) {
        return 0;
    }

    if (getThisId() == targetNode &&
            currentDir == targetDir &&
            currentDistance >= minDistance &&
            currentDistance <= maxDistance) {
        nPathFound++;
        path[0] = currentDistance - shift;
        pathLengthMean += currentDistance - shift;
        pathLengthSD += (currentDistance - shift)*(currentDistance - shift);
        //  VP.insert(path);

        //search continue since another valid path could be found
    }

    //cancel the search if too much branching
    //
    //ln(0)=nan
    //ln(1)=0
    if (currentDir) {
        if (getNRightOverlap() > 1)
            cancel += fastLog(getNRightOverlap());
    } else {
        if (getNLeftOverlap() > 1)
            cancel += fastLog(getNLeftOverlap());
    }

    //    if (cancel > maxLog)
    //        maxLog=cancel;

    if (cancel > maxLogNLeaves) {
        cancel = -1.0; //stop
        return 0;
    }

    if (currentDir) {
        for (unsigned int i = 0; i < getNRightOverlap(); i++) {
            getRightOverlap(i, id, size, direction);

            //            if (N[id].value == 0)
            //                continue;

            //distance = currentDistance + N[id].getSequenceLength()-size;
            distance = currentDistance + getSequenceLength() - size;

            path.push_back(i);

            if (direction) {
                retVal = N[id].searchPathToPaired(currentDir, targetNode, targetDir, minDistance, maxDistance, shift, distance, cancel, maxLogNLeaves);
            } else {
                retVal = N[id].searchPathToPaired(!currentDir, targetNode, targetDir, minDistance, maxDistance, shift, distance, cancel, maxLogNLeaves);
            }

            path.pop_back();
        }
    } else {
        for (unsigned int i = 0; i < getNLeftOverlap(); i++) {
            getLeftOverlap(i, id, size, direction);
            //distance = currentDistance + N[id].getSequenceLength()-size;
            distance = currentDistance + getSequenceLength() - size;

            path.push_back(i + nOvRight);

            if (direction) {
                retVal = N[id].searchPathToPaired(currentDir, targetNode, targetDir, minDistance, maxDistance, shift, distance, cancel, maxLogNLeaves);
            } else {
                retVal = N[id].searchPathToPaired(!currentDir, targetNode, targetDir, minDistance, maxDistance, shift, distance, cancel, maxLogNLeaves);
            }

            path.pop_back();
        }
    }

    if (cancel != -1.0) {
        if (currentDir) {
            if (getNRightOverlap() > 1)
                cancel -= fastLog(getNRightOverlap());

            if (cancel < 0.0)
                cancel = 0.0;
        } else {
            if (getNLeftOverlap() > 1)
                cancel -= fastLog(getNLeftOverlap());

            if (cancel < 0.0)
                cancel = 0.0;
        }
    }

    return retVal;
}

bool Node::isBranching()//true if branching on the right or left hand side
{
    //false if internal or ending or alone

    if (getNRightOverlap() > 1)
        return true;
    if (getNLeftOverlap() > 1)
        return true;

    //Nodes with a single self reverse overlaps is a branching node
    if (getNRightOverlap() == 1) {
        if (ovId[0] == getThisId() && ovSize[0] < 1)
            return true;
    }
    if (getNLeftOverlap() == 1) {
        if (ovId[nOvRight] == getThisId() && ovSize[nOvRight] < 1)
            return true;
    }

    return false;
}

bool Node::isBranching(bool front) {
    if (front) {
        if (getNRightOverlap() > 1)
            return true;

        //self overlaps
        if (getNRightOverlap() == 1) {
            if (ovId[0] == getThisId() && ovSize[0] < 1)
                return true;
        }
    } else {
        if (getNLeftOverlap() > 1)
            return true;

        //self overlaps
        if (getNLeftOverlap() == 1) {
            if (ovId[nOvRight] == getThisId() && ovSize[nOvRight] < 1)
                return true;
        }
    }

    return false;
}

bool Node::isInternal() //a single OV on both side (=not branching and not ending)
{
    if (!isBranching() && !isEnding())
        return true;
    else
        return false;
}

bool Node::hasSingleEdge(bool front) {
    if (front) {
        if (getNRightOverlap() == 1) {
            if (ovId[0] == getThisId() && ovSize[0] < 1)
                return false;
            else
                return true;
        }
    } else {
        if (getNLeftOverlap() == 1) {
            if (ovId[nOvRight] == getThisId() && ovSize[nOvRight] < 1)
                return false;
            else
                return true;
        }
    }
    return false;
}

bool Node::isEnding() {
    if (getNRightOverlap() == 0 || getNLeftOverlap() == 0)
        return true;
    else
        return false;
}

bool Node::isEnding(bool dir) {
    if (dir) {
        if (getNRightOverlap() == 0) {
            return true;
        }
    } else {
        if (getNLeftOverlap() == 0) {
            return true;
        }
    }
    return false;
}

bool Node::hasNoEdge() {
    return (getNLeftOverlap() == 0 && getNRightOverlap() == 0);
}

bool Node::hasNeighbors() {//other than itself (self overlap)
    if (hasNoEdge())
        return false;

    for (size_t i = 0; i < getNOverlap(); i++) {
        if (ovId[i] != getThisId())
            return true;
    }

    return false;
}

bool Node::hasSingleSideOV()//Ov on a single side
{
    return ( (getNLeftOverlap() == 0 && getNRightOverlap() > 0) ||
            (getNLeftOverlap() > 0 && getNRightOverlap() == 0));
}

bool Node::hasBothSidesOV() {
    return (getNLeftOverlap() > 0 && getNRightOverlap() > 0);
}


