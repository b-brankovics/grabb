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

#include "overlapGraph.h"
#include "readsLayout.h"
#include "globalFunc.h"
#include "stat.h"
#include "NodeIt.h"
#include "crc.h"
#include "logWriter.h"
#include <cmath>
#include <algorithm>
#include <set>
#include <sstream>
#include <limits>
#include <ctime>
#include <pthread.h>

extern bool DEV_INFO;
extern logWriter LOG;
bool PRINT_TREE;

Node* OverlapsGraph::nodesTab = 0x0;
AssemblyProgress OverlapsGraph::AP;
deque<unsigned int> OverlapsGraph::nodeQueue;
unsigned int OverlapsGraph::g_count;

OverlapsGraph::OverlapsGraph() {
    nodesTab = 0x0;
    nNodes = 0;
    nEdges = 0;
    Node::N = 0x0;
    Node::nNodes = nNodes;
    depthLimit = 10; //default
    targetSizeEstimation = 0;
}

OverlapsGraph::~OverlapsGraph() {
    freeMemory();
}

void OverlapsGraph::init(ReadsStorage *r, Pairing *p, ReadsLayout *l) {
    R = r;
    P = p;
    L = l;
    Node::R = R;
    Node::L = L;
    Node::P = P;
}

void OverlapsGraph::allocateNodeTab(unsigned int maxN) {
    maxNodes = maxN;
    try {
        nodesTab = new Node[maxN + 1];
        Node::N = nodesTab;
        Node::G = this;
        Node::R = R;
        Node::L = L;
        Node::P = P;
        NodeIt::N = nodesTab;
    } catch (bad_alloc ex) {
        cout << ex.what() << "\not enough memory. OverlapsGraph::allocateNodeTab(unsigned int)" << '\n';
        exit(0);
    }
    //  Node::VP.allocate((size_t)R->getEffectiveNReads());//should be ok
}

void OverlapsGraph::freeMemory() {
    if (nodesTab != 0x0) {
        for (unsigned int i = 0; i < maxNodes + 1; i++) {
            nodesTab[i].freeMemory();
        }

        delete[] nodesTab;
    }

    //  Node::VP.cleanMemory();

    nodesTab = 0x0;
    Node::N = 0x0;
    Node::G = 0x0;
    Node::R = 0x0;
    Node::L = 0x0;
    NodeIt::N = 0x0;
    nNodes = 0;
}

bool OverlapsGraph::loadData(string file) {
    bool checkFailed = false;

    ifstream in_bin;
    in_bin.open(file.c_str(), ios_base::binary);
    if (!in_bin) {
        cout << "Cannot open " << file << endl;
        return false;
    }

    char *tag = new char[13];
    in_bin.read(tag, 12 * sizeof (char));
    tag[12] = '\0';
    string tmp = tag;

    bool edena302 = false;

    if (tmp == "EDENA302####") {
        edena302 = true;
    } else
        if (tmp != "EDENA303####") {
        if (tmp.substr(0, 5) == "EDENA") {
            cout << "The file \"" << file << "\" is not compatible with current version.\n";
            cout << "Please, generate a new ovl file\n";
            return false;
        } else {
            cout << "The file " << file << " does not appear to be an EDENA file\n";
            return false;
        }
    }

    ovlVersion = tmp;

    freeMemory();

    cout << "Loading file \"" << file << "\"" << "... " << flush;

    if (checkFailed == false)
        if (R->load(in_bin) == false)
            checkFailed = true;

    if (checkFailed == false)
        if (P->load(in_bin) == false)
            checkFailed = true;

    if (checkFailed == false)
        if (L->load(in_bin, R) == false)
            checkFailed = true;

    if (checkFailed == false) {
        Crc32 CRC;
        unsigned int crcCheck;

        in_bin.read((char*) (&nNodes), sizeof (unsigned int));
        in_bin.read((char*) (&minOverlap), sizeof (unsigned int));

        CRC.AddData((uint8_t*) & nNodes, sizeof (unsigned int));
        CRC.AddData((uint8_t*) & minOverlap, sizeof (unsigned int));

        unsigned int nOverlaps = 0;
        unsigned int nOv, nOvr;
        unsigned int mylayout;

        allocateNodeTab(nNodes);
        Node::nNodes = nNodes;

        unsigned int crc;

        for (unsigned int i = 1; i <= nNodes; i++) {

            in_bin.read((char*) &mylayout, sizeof (unsigned int));
            in_bin.read((char*) (&nOv), sizeof (unsigned int));
            in_bin.read((char*) (&nOvr), sizeof (unsigned int));
            //allocate memory for nOv overlaps
            nOverlaps += nOv;
            nodesTab[i].allocate(nOv);
            nodesTab[i].allocateEdgeValue();
            nodesTab[i].layout = mylayout;
            nodesTab[i].nOvRight = nOvr;

            if (nOv > 0) {
                in_bin.read((char*) nodesTab[i].ovId, nOv * sizeof (unsigned int));
                in_bin.read((char*) nodesTab[i].ovSize, nOv * sizeof (short int));
                in_bin.read((char*) nodesTab[i].edgeBitsFlag, nOv * sizeof (char));
                CRC.AddData((uint8_t*) nodesTab[i].ovId, nOv * sizeof (unsigned int));
                CRC.AddData((uint8_t*) nodesTab[i].ovSize, nOv * sizeof (short int));
                CRC.AddData((uint8_t*) nodesTab[i].edgeBitsFlag, nOv * sizeof (char));
            }

            if (i % 10000 == 0) {
                if (edena302)
                    continue;

                crc = CRC.GetCrc32();
                in_bin.read((char*) &crcCheck, sizeof (unsigned int));
                if (crc != crcCheck) {
                    checkFailed = true;
                    break;
                }

            }
        }

        in_bin.read((char*) &crcCheck, sizeof (unsigned int));

        crc = CRC.GetCrc32();
        if (crc != crcCheck)
            checkFailed = true;

        nEdges = nOverlaps / 2;
        in_bin.close();
        delete [] tag;

    }

    if (checkFailed) {
        cout << "The file is corrupted, CRC check failed\n"
                << "Please generate a new ovl file\n";
        return false;
    }

    Node::edgeSorted = true;

    cout << " done\n";
    return true;
}

bool OverlapsGraph::saveData(string file) {
    ofstream out_bin;
    out_bin.open(file.c_str(), ios_base::binary);
    if (!out_bin)
        return false;
    cout << "writing the overlaps graph to the file " << file << "... " << flush;
    out_bin << "EDENA303####"; //12 bytes

    R->save(out_bin);
    P->save(out_bin);
    L->save(out_bin);

    Crc32 CRC;
    unsigned int crcCheck;

    out_bin.write((char*) &nNodes, sizeof (unsigned int));
    out_bin.write((char*) &minOverlap, sizeof (unsigned int));

    CRC.AddData((uint8_t*) & nNodes, sizeof (unsigned int));
    CRC.AddData((uint8_t*) & minOverlap, sizeof (unsigned int));

    unsigned int nOv, nOvr;
    for (unsigned int i = 1; i <= nNodes; i++) {

        out_bin.write((char*) &nodesTab[i].layout, sizeof (unsigned int));
        nOv = nodesTab[i].nOv;
        out_bin.write((char*) (&nOv), sizeof (unsigned int));
        nOvr = nodesTab[i].nOvRight;
        out_bin.write((char*) (&nOvr), sizeof (unsigned int));

        if (nOv > 0) {
            out_bin.write((char*) (nodesTab[i].ovId), nOv * sizeof (unsigned int));
            out_bin.write((char*) (nodesTab[i].ovSize), nOv * sizeof (short int));
            out_bin.write((char*) (nodesTab[i].edgeBitsFlag), nOv * sizeof (char));
            CRC.AddData((uint8_t*) nodesTab[i].ovId, nOv * sizeof (unsigned int));
            CRC.AddData((uint8_t*) nodesTab[i].ovSize, nOv * sizeof (short int));
            CRC.AddData((uint8_t*) nodesTab[i].edgeBitsFlag, nOv * sizeof (char));
        }

        if (i % 10000 == 0) {
            crcCheck = CRC.GetCrc32();
            out_bin.write((char*) &crcCheck, sizeof (unsigned int));
        }
    }

    crcCheck = CRC.GetCrc32();
    out_bin.write((char*) &crcCheck, sizeof (unsigned int));

    cout << "done" << endl;

    out_bin.close();
    return true;
}

unsigned int OverlapsGraph::getNOverlap(unsigned int id) {
    return nodesTab[id].nOv;
}

unsigned int OverlapsGraph::getNRightOverlap(unsigned int id) {
    return nodesTab[id].nOvRight;
}

unsigned int OverlapsGraph::getNLeftOverlap(unsigned int id) {
    return nodesTab[id].nOv - nodesTab[id].nOvRight;
}

unsigned int OverlapsGraph::getReadLength() {
    return R->readsLength;
}

void* OverlapsGraph::thread_maker(void* dat_s) {
    Thread_param* v = static_cast<Thread_param*> (dat_s);
    OverlapsGraph* e = v->_this;
    void ((OverlapsGraph::*f)(void*)) = v->funcPTR;
    (e->*f)(dat_s);
    return NULL;
}

void OverlapsGraph::computeOverlaps(unsigned int nThreads) {
    minOverlap = R->minOvSize;
    nNodes = R->getN_nrReads();
    Node::nNodes = nNodes;

    allocateNodeTab(nNodes);
    g_count++;

    //Set layout in nodes (one read per node at this point)
    for (unsigned int i = 1; i <= nNodes; i++)
        nodesTab[i].layout = L->getLastIdentical(i);

    if (nThreads > nNodes)
        nThreads = nNodes;

    int chunkSize = nNodes / nThreads;
    if (chunkSize == 0)
        chunkSize = 1;

    pthread_t *threads = new pthread_t[nThreads];
    Thread_param *param = new Thread_param[nThreads];

    unsigned int start, end;
    for (unsigned int nn = 0; nn < nThreads; nn++) {
        start = (nn * chunkSize + 1);
        if (nn == nThreads - 1)
            end = nNodes;
        else
            end = (nn * chunkSize + 1) + chunkSize - 1;

        param[nn]._this = this;
        param[nn].funcPTR = &OverlapsGraph::overlapNodeRange;
        param[nn].start = start;
        param[nn].end = end;
        pthread_create(&threads[nn], NULL, thread_maker, (void *) &param[nn]);
    }

    cout << "Computing overlaps >=" << minOverlap << "... " << flush;
    for (unsigned int nn = 0; nn < nThreads; nn++) {
        pthread_join(threads[nn], NULL);
    }
    cout << "\rComputing overlaps >=" << minOverlap << " done        \n" << flush;
    R->cleanPrefixTables();
    Node::edgeSorted = true;
}

void OverlapsGraph::overlapNodeRange(void* ptr) {
    static pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
    static pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
    Thread_param *param = (Thread_param*) ptr;
    vector<unsigned int> ID, ID2;
    vector<short> SIZE, SIZE2;
    multiset<_OV, orderSet> ovSet;
    static clock_t t;
    unsigned int nOv = 0;
    unsigned int nOvRight = 0;
    static unsigned int count;
    static int activeThreads;

    pthread_mutex_lock(&mutex1);
    activeThreads++;
    pthread_mutex_unlock(&mutex1);

    for (unsigned i = param->start; i <= param->end; i++) {
        if (nodesTab[i].isAllocated())
            continue;

        pthread_mutex_lock(&mutex2);

        //~every 1/4 seconds if all threads run at 100%
        cerr << setprecision(3) << fixed;
        if ((clock() - t) >= activeThreads * CLOCKS_PER_SEC / 4) {
            //if (count % 10000 == 0) {
            t = clock();
            cerr << "\rComputing overlaps >=" << minOverlap << "...";
            cerr << count << " (" << 100.0 * count / this->getNNodes() << "%)"
                    << "                              " << flush;
        }

        count++;
        pthread_mutex_unlock(&mutex2);

        // nOvRight=R->determineOverlaps(i,ID,SIZE,ovSet);

        nOvRight = R->determineOverlaps2(i, ID, SIZE, ovSet);

        //        if (ID.size() != ID2.size())
        //            cout << "";
        //        else
        //        for (unsigned int j=0; j<ID.size(); j++)
        //        {
        //            if (ID[j]!=ID2[j])
        //                cout << "";
        //            if (SIZE[j]!=SIZE2[j])
        //                cout << "";
        //        }

        nOv = ID.size();

        nodesTab[i].allocate(nOv);
        nodesTab[i].nOv = nOv;
        nodesTab[i].nOvRight = nOvRight;

        for (unsigned int ii = 0; ii < nOv; ii++) {
            nodesTab[i].ovId[ii] = ID[ii];
            nodesTab[i].ovSize[ii] = SIZE[ii];
        }

    }
    pthread_mutex_lock(&mutex1);
    activeThreads--;
    pthread_mutex_unlock(&mutex1);
}

void OverlapsGraph::condense(bool verbose, bool sortEdges) {

    string concatSeq;
    vector<unsigned int> coverage;
    unsigned int leftId, rightId;
    bool leftDir, rightDir;

    unsigned int uCounter = 0;
    unsigned int nOv, nOvRight;

    vector<unsigned int> ovID;
    vector<short> ovSize;
    vector<float> ovValue;
    vector<unsigned char> ovFlags;

    unsigned int* newId = new unsigned int [nNodes + 1];
    bool* newDir = new bool [nNodes + 1];
    memset(newId, 0, (nNodes + 1) * sizeof (unsigned int));

    Node::S_counter = 0;
    Node::edgeSorted = false;

    for (unsigned int i = 1; i <= nNodes; i++) {
        nodesTab[i].unsetExtended();
        nodesTab[i].unsetVisited();
    }

    unsigned int inter = nNodes / 100;
    unsigned int prev = 0;

    for (unsigned int i = 1; i <= nNodes; i++) {

        if (verbose && Node::S_counter - prev > inter) {
            prev = Node::S_counter;
            cerr << "Condensing overlaps graph..." << (int) (((float) Node::S_counter / nNodes)*100) << "%\r" << flush;
        }

        if (nodesTab[i].isDiscarded()) {
            //discarded nodes must be orphan
            if (nodesTab[i].nOv != 0) {
                cerr << "[bug] OverlapsGraph::condense(...)\n";
                sendBugReportPlease(cerr);
            }

            L->setLayoutNodeId(nodesTab[i].layout, 0);
            continue;
        }

        if (!nodesTab[i].isVisited()) {
            nodesTab[i].getLongestNonAmbiguousPath(
                    concatSeq,
                    coverage,
                    leftId,
                    leftDir,
                    rightId,
                    rightDir);
            uCounter++;

            nodesTab[i].setExtended();

            newId[i] = i;
            newId[leftId] = i;
            newId[rightId] = i;
            newDir[i] = true;
            newDir[leftId] = leftDir;
            newDir[rightId] = rightDir;

            ovID.clear();
            ovSize.clear();
            ovValue.clear();
            ovFlags.clear();

            if (rightDir == true) {
                nOvRight = nodesTab[rightId].getNRightOverlap();
                for (unsigned int ii = 0; ii < nOvRight; ii++) {
                    ovID.push_back(nodesTab[rightId].ovId[ii]);
                    ovSize.push_back(nodesTab[rightId].ovSize[ii]);
                    ovValue.push_back(nodesTab[rightId].edgeValue[ii]);
                    ovFlags.push_back(nodesTab[rightId].edgeBitsFlag[ii]);
                }
            } else {
                nOvRight = nodesTab[rightId].getNLeftOverlap();
                for (unsigned int ii = nodesTab[rightId].getNRightOverlap(); ii < nodesTab[rightId].getNOverlap(); ii++) {
                    ovID.push_back(nodesTab[rightId].ovId[ii]);
                    ovSize.push_back(-(nodesTab[rightId].ovSize[ii]));
                    ovValue.push_back(nodesTab[rightId].edgeValue[ii]);
                    ovFlags.push_back(nodesTab[rightId].edgeBitsFlag[ii]);
                }
            }

            if (leftDir == true) {
                nOv = nOvRight + nodesTab[leftId].getNLeftOverlap();

                for (unsigned int ii = nodesTab[leftId].getNRightOverlap(); ii < nodesTab[leftId].getNOverlap(); ii++) {
                    ovID.push_back(nodesTab[leftId].ovId[ii]);
                    ovSize.push_back(nodesTab[leftId].ovSize[ii]);
                    ovValue.push_back(nodesTab[leftId].edgeValue[ii]);
                    ovFlags.push_back(nodesTab[leftId].edgeBitsFlag[ii]);
                }
            } else {
                nOv = nOvRight + nodesTab[leftId].getNRightOverlap();
                for (unsigned int ii = 0; ii < nodesTab[leftId].getNRightOverlap(); ii++) {
                    ovID.push_back(nodesTab[leftId].ovId[ii]);
                    ovSize.push_back(-(nodesTab[leftId].ovSize[ii]));
                    ovValue.push_back(nodesTab[leftId].edgeValue[ii]);
                    ovFlags.push_back(nodesTab[leftId].edgeBitsFlag[ii]);
                }
            }

            //reconnect extended node (i) to new lest and right neighbors)
            //Do not take care of the reciprocal edges at the moment

            nodesTab[i].freeMemory();

            if (nOv > 0) {
                nodesTab[i].allocate(nOv);
                nodesTab[i].allocateEdgeValue();
                nodesTab[i].nOvRight = nOvRight;

                for (size_t ii = 0; ii < nOv; ii++) {
                    nodesTab[i].ovId[ii] = ovID[ii];
                    nodesTab[i].ovSize[ii] = ovSize[ii];
                    nodesTab[i].edgeValue[ii] = ovValue[ii];
                    nodesTab[i].edgeBitsFlag[ii] = ovFlags[ii];
                }
            }
        }
    }


    if (verbose) {
        cout << "Condensing overlaps graph... done          " << endl;
        cout << "   Updating node IDs...\r" << flush;
    }

    //update reciprocal edges to correct IDs
    for (unsigned int i = 1; i <= nNodes; i++) {
        if (nodesTab[i].isExtended()) {
            for (unsigned int ii = 0; ii < nodesTab[i].getNOverlap(); ii++) {
                //update edges to correct IDs
                if (newDir[ nodesTab[i].ovId[ii] ] == false) {
                    nodesTab[i].ovSize[ii] = -nodesTab[i].ovSize[ii];
                }
                nodesTab[i].ovId[ii] = newId[ nodesTab[i].ovId[ii] ];
            }
        }
    }

    if (verbose)
        cout << "   Updating node IDs... done" << endl;

    //required
    renumber(verbose, sortEdges);

    delete [] newId;
    delete [] newDir;
}

void OverlapsGraph::renumber(bool verbose, bool sortEdges) {
    //nodes marked as Extended are kept and renumbered
    //others are discarded
    if (verbose)
        cout << "   Renumbering nodes..." << flush;
    unsigned int newIdCounter = 0;

    unsigned int* newId = new unsigned int [nNodes + 1];
    memset(newId, 0, (nNodes + 1) * sizeof (unsigned int));

    for (unsigned int i = 1; i <= nNodes; i++) {
        if (nodesTab[i].isExtended()) {

            newIdCounter++;

            nodesTab[i].unsetExtended();
            newId[i] = newIdCounter;

            if (i != newIdCounter) {
                //nodesTab[newIdCounter].freeMemory();//not required
                memcpy(&nodesTab[newIdCounter], &nodesTab[i], sizeof (Node));
                L->setLayoutNodeId(nodesTab[newIdCounter].layout, newIdCounter);
                //do not call node::freeMemory() here!!
                nodesTab[i].ovId = 0x0;
                nodesTab[i].ovSize = 0x0;
                nodesTab[i].edgeValue = 0x0;
                nodesTab[i].edgeBitsFlag = 0x0;
                nodesTab[i].nOv = nodesTab[i].nOvRight = 0;
                nodesTab[i].layout = 0;
            }
        } else //not extended so do not take part to new condensed graph
        {
            nodesTab[i].freeMemory();
        }
    }

    nNodes = newIdCounter;
    Node::nNodes = newIdCounter;
    nEdges = 0;
    unsigned int nE;

    for (unsigned int i = 1; i <= nNodes; i++) {
        nE = nodesTab[i].getNOverlap();
        nEdges += nE;
        for (unsigned int ii = 0; ii < nE; ii++) {
            nodesTab[i].ovId[ii] = newId[ nodesTab[i].ovId[ii] ];
        }

        if (sortEdges)
            nodesTab[i].sortByIds();
    }

    nEdges /= 2;

    if (sortEdges)
        Node::edgeSorted = true;

    if (verbose) {
        cout << " done" << endl;
    }

    delete [] newId;
}

void OverlapsGraph::printSummary(unsigned int id, bool direct, bool right) {
    string seq;
    unsigned int formatLength = 70;
    string fill;

    seq = nodesTab[id].getSequence(direct);


    if (seq.length() < formatLength) {
        fill = "";
        for (size_t i = 0; i < formatLength - seq.length(); i++)
            fill += ' ';
        if (right)
            seq = fill + seq;
        else
            seq = seq + fill;
    } else if (seq.length() > formatLength) {
        if (right) {
            seq = seq.substr((seq.length() - formatLength) + 3);
            seq = "..." + seq;

        } else {
            seq = seq.substr(0, formatLength - 3);
            seq += "...";
        }
    }

    if (direct) {
        cout << seq << " (+)";
        cout << " id=" << id;
        cout << " l=" << nodesTab[id].getSequenceLength();
        //     cout << " " << nodesTab[id].getNLeftOverlap() << "<->";
        //        cout << nodesTab[id].getNRightOverlap();
    } else {
        cout << seq << " (-)";
        cout << " id=!" << id;
        cout << " l=" << nodesTab[id].getSequenceLength();
        //       cout << " " << nodesTab[id].getNRightOverlap() << "<->";
        //        cout << nodesTab[id].getNLeftOverlap();
    }

}

void OverlapsGraph::nicePrint(unsigned int id, bool direct) {

    unsigned int neighbor, size;
    bool direction;


    cout << "Node ";
    if (!direct)
        cout << '!';
    cout << id << endl;
    cout << "   length: " << nodesTab[id].getSequenceLength() << endl;
    cout << "   cov:    " << nodesTab[id].getCoverage() << endl;
    cout << "   nReads: " << nodesTab[id].getNReads() << endl;
    cout << "   nEdges: " << nodesTab[id].getNOverlap() << endl;

    for (unsigned int i = 0; i < nodesTab[id].getNOverlap(); i++) {
        if (i >= nodesTab[id].getNLeftOverlap())
            cout << "      > ";
        else
            cout << "      < ";

        nodesTab[id].getOverlap(i, neighbor, size, direction);
        if (direction == false)
            cout << "!";

        cout << neighbor << " ovlp=" << size << "\n";
    }
}

unsigned int OverlapsGraph::countEdges() {
    nEdges = 0;
    for (unsigned int i = 1; i <= nNodes; i++) {
        nEdges += nodesTab[i].nOv;
    }
    nEdges /= 2;
    return nEdges;
}

void OverlapsGraph::computeEdgesProb(double reliableCutoff) {
    cout << "Contextual cleaning: step1...\r" << flush;

    for (size_t i = 1; i <= nNodes; i++) {
        nodesTab[i].initializeEdgeValues(-1.0);
        //(default value)  nodesTab[i].initEdgeFlags(false); //used to tell whether the edge has at least one
        //reliable brother
        nodesTab[i].unsetVisited();
    }

    for (size_t i = 1; i <= nNodes; i++) {
        if (i % 50 == 0)
            cerr << setprecision(2) << "Contextual cleaning: step1..."
            << (double) 100 * i / nNodes
                << "%           \r" << flush;


        //  nodesTab[i].computeEdgesProb(200,50, reliableCutoff);
        //  nodesTab[i].computeEdgesProb(200,100, reliableCutoff);
        nodesTab[i].computeEdgesProb(100, 100, reliableCutoff);
    }
    cout << "Contextual cleaning: step1... done                    " << endl;

}

unsigned int OverlapsGraph::removeEdgesByValue(float suspectCutoff, float GIncoherentCutoff) {

    unsigned int nRemoved = 0;
    cout << "Contextual cleaning: step2...\r" << flush;

    for (size_t i = 1; i <= nNodes; i++) {
        //        if (i % 50 == 0)
        //            cout << setprecision(2) << "Contextual cleaning: step2..."
        //            << (double) 100 * i / nNodes
        //                << "%           \r" << flush;

        nRemoved += nodesTab[i].removeEdgesByValue(suspectCutoff, GIncoherentCutoff);
    }

    cout << "Contextual cleaning: step2... done                 " << endl;
    return nRemoved / 2; //self overlap do not have reciprocal edges !!
}

void OverlapsGraph::overlapSizeCutoff(unsigned int cutoff) {
    for (unsigned int i = 1; i <= nNodes; i++) {
        nodesTab[i].overlapSizeCutoff(cutoff);
    }
    minOverlap = cutoff;
}

unsigned int OverlapsGraph::coverageCutoff(unsigned int cutoff) {
    unsigned int count = 0;
    for (unsigned int i = 1; i <= nNodes; i++) {
        if (nodesTab[i].getCoverage() <= cutoff) {
            nodesTab[i].isolate();
            count++;
        }
    }

    return count;
}

void OverlapsGraph::removeTransitiveEdges() {

    for (unsigned int i = 1; i <= nNodes; i++)
        memset(nodesTab[i].edgeBitsFlag, 0, nodesTab[i].getNOverlap() * sizeof (unsigned char));

    cout << "Flagging transitive edges...\r" << flush;

    ofstream out_bin;
    ifstream in_bin;
    // out_bin.open("debugDATA", ios_base::binary);
    //  in_bin.open("debugDATA", ios_base::binary);
    vector<bool> flagged;

    for (unsigned int i = 1; i <= nNodes; i++) {
        if (i % 1000 == 0) {
            cerr << "Flagging transitive edges..." << (int) (((float) i / nNodes)*100) << "%\r" << flush;
        }

        //     flagged.clear();


        nodesTab[i].markTransitiveEdges();
        //      
        //     unsigned int nnov=nodesTab[i].getNOverlap();
        //     out_bin.write((char*)&nnov,sizeof(unsigned int) );
        //       if (nodesTab[i].getNOverlap() > 0)
        //       out_bin.write((char*) ((nodesTab[i].edgeBitsFlag)), nodesTab[i].getNOverlap()*sizeof (unsigned char));


        //        unsigned char p;
        //      flagged.clear();
        //      unsigned int nnov;
        //      in_bin.read( (char*) &nnov, sizeof(unsigned int));
        //      for (unsigned int j=0; j<nnov; j++)
        //      {
        //          in_bin.read( (char*)(&p), sizeof(unsigned char));
        //          flagged.push_back(p);
        //      }

        // nodesTab[i].anotherMarkTransitiveEdges2(flagged);

    }
    cout << "Flagging transitive edges... done    " << endl;
    //out_bin.close();
    for (unsigned int i = 1; i <= nNodes; i++) {
        if (i % 1000 == 0) {
            cerr << "Removing flagged edges..." << (int) (((float) i / nNodes)*100) << "%\r" << flush;
        }
        nodesTab[i].removeMarkedEdge();
    }
    cout << "Removing flagged edges... done" << endl;
}

bool OverlapsGraph::checkConsistency() {
    cout << "Performing a consistency check...\r" << flush;

    for (unsigned int i = 1; i <= nNodes; i++) {
        if (i % 1000 == 0) {
            cerr << "Performing a consistency check..." << (int) (((float) i / nNodes)*100) << "%\r" << flush;
        }

        if (!nodesTab[i].testReciprocal()) {
            cout << "\nfailed\n";
            return false;
        }
    }
    cout << "Performing a consistency check... successful        \n";
    return true;
}

unsigned int OverlapsGraph::discardShortOrphanNodes(unsigned int &nDiscardedReads) {
    unsigned int nDiscardedNodes = 0;
    nDiscardedReads = 0;

    unsigned int maxSize = R->getReadsLength() + R->getReadsLength() / 2;

    for (size_t i = 1; i <= getNNodes(); i++) {
        if (nodesTab[i].getSequenceLength() < maxSize &&
                nodesTab[i].hasNeighbors() == false) {
            nodesTab[i].removeAllEdges(); //should be only self overlaps
            //   nodesTab[i].isolate(); //required edges to be sorted
            nodesTab[i].setDiscarded();

            nDiscardedNodes++;
            nDiscardedReads += nodesTab[i].getNReads();
        }
    }

    return nDiscardedNodes;
}

unsigned int OverlapsGraph::discardSuspiciousNode(double minCov) {
    //suspicious nodes are defined as follows:
    //node coverage < minCov
    //connected, with in- and out-degree no more that one
    //neighbor node(s) coverage at least 20x the suspicious node
    //neighbor node as another neighbor (=no disruption)

    //such nodes are mainly dead-ends that were too long to be removed
    //by the conventional dead-end removal procedure.

    unsigned int neighbor1, neighbor2;
    bool ok1, ok2;
    unsigned int size;
    bool dir;
    double cov;
    unsigned int nSuspicious = 0;

    for (unsigned int i = 1; i <= getNNodes(); i++) {
        cov = nodesTab[i].getCoverage();
        if (cov >= minCov)
            continue;

        if (nodesTab[i].getSequenceLength() > 3 * getReadLength())
            continue;

        neighbor1 = 0;
        neighbor2 = 0;
        ok1 = false;
        ok2 = false;

        if (nodesTab[i].getNRightOverlap() == 1) {
            nodesTab[i].getRightOverlap(0, neighbor1, size, dir);

            if (nodesTab[neighbor1].getNEdges(!dir) < 2)
                continue;

            if (nodesTab[neighbor1].getCoverage() > 20 * cov)
                ok1 = true;
        }

        if (nodesTab[i].getNLeftOverlap() == 1) {
            nodesTab[i].getLeftOverlap(0, neighbor2, size, dir);

            if (nodesTab[neighbor2].getNEdges(dir) < 2)
                continue;

            if (nodesTab[neighbor2].getCoverage() > 20 * cov)
                ok2 = true;
        }

        if (neighbor1 != 0 && ok1 == false)
            continue;
        if (neighbor2 != 0 && ok2 == false)
            continue;
        if (neighbor1 == 0 && neighbor2 == 0)
            continue;

        nSuspicious++;

        //print out removed nodes
        //        ostringstream oss;
        //        oss << "suspicious_" << i;         
        //        nodesTab[i].dotLocalGraph(4,oss.str());

        nodesTab[i].setDiscarded();
        nodesTab[i].isolate();
    }

    return nSuspicious;
}

unsigned int OverlapsGraph::identifyDeadEnd(unsigned int &dLimit, unsigned int &nrreads) {

    unsigned int nDeadEnds;
    unsigned int nTotalDeadEnds = 0;
    nrreads = 0;

    if (dLimit == 0)
        dLimit = 2 * R->getReadsLength() - 1;

    //     do {
    //                nDeadEnd = G.identifyDeadEnd(param.maxDeadEndLength, nReads);
    //                sum += nDeadEnd;
    //            } while (nDeadEnd != 0);

    do { //the whole process is iterated until no dead-end remains
        Node::nodeList.clear();
        nDeadEnds = 0;

        for (unsigned int i = 1; i <= getNNodes(); i++) {

            if (nodesTab[i].isEnding(true)) {
                identifyDeadEnd(i, false, nodesTab[i].getSequenceLength(), nDeadEnds, dLimit);
            }

            if (nodesTab[i].isEnding(false)) {

                identifyDeadEnd(i, true, nodesTab[i].getSequenceLength(), nDeadEnds, dLimit);
            }
        }

        unsigned int nodeId;

        // for all edges
        for (unsigned int i = 1; i <= getNNodes(); i++) {
            for (unsigned int j = 0; j < nodesTab[i].getNOverlap(); j++) {

                nodeId = nodesTab[i].getOverlapId(j);
                if (nodesTab[nodeId].isDiscarded())
                    nodesTab[i].setEdgeFlagUnrec(j, true);
            }
            nodesTab[i].removeMarkedEdge();

        }

        for (vector<unsigned int>::iterator it = Node::nodeList.begin();
                it != Node::nodeList.end();
                it++) {
            // if (nodesTab[*it].getNOverlap() != 0)
            nrreads += nodesTab[*it].getNReads();

            nodesTab[*it].removeAllEdges(); //no need to remove the reciprocal as well
            //already removed during previous step
            // nodesTab[*it].isolate();

        }
        nTotalDeadEnds += nDeadEnds;

    } while (nDeadEnds != 0);

    return nTotalDeadEnds;
}

bool OverlapsGraph::identifyDeadEnd(unsigned int nodeId,
        bool dir,
        unsigned int distance,
        unsigned int &nDeadEnd,
        unsigned int dLimit) {
    //check dead-end length and mark nodes 
    unsigned int id, size;
    bool ovDir;
    bool markBack = false;

    if (distance > dLimit)
        return false;

    if (dir) {
        for (unsigned int i = 0; i < nodesTab[nodeId].getNRightOverlap(); i++) {
            nodesTab[nodeId].getRightOverlap(i, id, size, ovDir);

            if (nodesTab[id].isBranching(dir != ovDir)) {
                markBack = true;
                nDeadEnd++;
            } else {
                if (identifyDeadEnd(id,
                        dir == ovDir,
                        distance + nodesTab[id].getSequenceLength() - size,
                        nDeadEnd,
                        dLimit))
                    markBack = true;
            }
            break; //explore only the first edge is enough to remove all DE path
            //since the procedure is iterated until no DE remains
        }
    } else {
        for (unsigned int i = 0; i < nodesTab[nodeId].getNLeftOverlap(); i++) {
            nodesTab[nodeId].getLeftOverlap(i, id, size, ovDir);

            if (nodesTab[id].isBranching(dir != ovDir)) {
                markBack = true;
                nDeadEnd++;
            } else {
                if (identifyDeadEnd(id,
                        dir == ovDir,
                        distance + nodesTab[id].getSequenceLength() - size,
                        nDeadEnd,
                        dLimit))
                    markBack = true;
            }
            break; //explore only the first edge is enough to remove all DE path
            //since the procedure is iterated until no DE remains
        }
    }

    if (markBack) {
        if (nodesTab[nodeId].isDiscarded() == false) {
            Node::nodeList.push_back(nodeId);
            nodesTab[nodeId].setDiscarded();
        }
    }
    return markBack;
}

unsigned int OverlapsGraph::bubbles(double minCov) {
    unsigned int nBub = 0;

    cout << "Resolving bubbles..." << flush;
    for (size_t i = 1; i <= getNNodes(); i++) {
        nBub += nodesTab[i].bubble(true, minCov);
        nBub += nodesTab[i].bubble(false, minCov);
    }
    cout << " done\n";
    return nBub;

}

int OverlapsGraph::getPathBetweenTwoReads(
        unsigned int id1,
        bool d1,
        unsigned int id2,
        bool d2,
        int minLength,
        int maxLength,
        double maxLogNLeaves) {
    unsigned int node1, node2, p1, p2;
    bool dirInNode1, dirInNode2;
    int minD, maxD, retValue;

    //get reads position in the graph:
    node1 = L->getNodeId(id1);
    node2 = L->getNodeId(id2);

    //read position in target node
    p1 = L->getPosition(id1);
    p2 = L->getPosition(id2);

    //read direction in target node
    dirInNode1 = L->getDirection(id1);
    dirInNode2 = L->getDirection(id2);

    int deltaDist1, deltaDist2;

    if (dirInNode1 == d1)
        deltaDist1 = p1 - 1;
    else
        deltaDist1 = nodesTab[node1].getSequenceLength() - p1 - R->getReadsLength() + 1;

    if (dirInNode2 != d2)
        deltaDist2 = -nodesTab[node2].getSequenceLength() + p2 - 1;
    else
        deltaDist2 = -p2 - R->getReadsLength() + 1;


    maxD = maxLength + deltaDist1 + deltaDist2;
    minD = minLength + deltaDist1 + deltaDist2;

    if (!dirInNode1)
        d1 = !d1;
    if (!dirInNode2)
        d2 = !d2;

    int shift = deltaDist1 + deltaDist2;

    retValue = nodesTab[node1].searchPathToPaired(dirInNode1, node2, !dirInNode2, minD, maxD, shift, maxLogNLeaves);


    //
    //    for (size_t i = 0; i < Node::VP.getNPath(); i++)
    //    {
    //
    //        pathLength = *(Node::VP.getPath(i));
    //        pathLength -= (deltaDist1 + deltaDist2);
    //        Node::VP.setPathLength(i,pathLength);
    //        Node::pathLengthMean += pathLength;
    //        Node::pathLengthSD += (pathLength * pathLength);
    //    }
    //
    //    Node::pathLengthMean /= Node::VP.getNPath();
    //
    //     if (Node::VP.getNPath() > 1)
    //    {
    //         Node::pathLengthSD =
    //                sqrt(1.0 / (Node::VP.getNPath() - 1)) *
    //                (
    //                Node::pathLengthSD -
    //                Node::VP.getNPath() * Node::pathLengthMean * Node::pathLengthMean
    //                );
    //    }
    //    else
    //        Node::pathLengthSD = 0.0;

    return retValue;
}

void OverlapsGraph::estimatePairedDistance(unsigned int PEhorizon, unsigned int MPhorizon, double nsd, string prefix) {
    int retValue;
    unsigned int id1, id2;

    unsigned int pathLength;
    unsigned int nUnique = 0, nMultiple = 0, nCanceled = 0, nZero = 0;
    double maxLogNLeaves = 3.0;
    unsigned int nCounted = 0;
    unsigned int globalCounter = 0;
    bool mate1Dir = true, mate2Dir = false;
    double meanLength, sdLength;
    unsigned int usableMateEstimate;

    if (P->getNPairing() == 0)
        return;

    //    unsigned int PEhorizon = 1000;
    //    unsigned int MPhorizon = 15000;

    //here reverse direct pairing are supposed to be long Illumina mate-pair.

    //    if (horizon != 0) {
    //        PEhorizon = horizon;
    //        MPhorizon = horizon;
    //        LOG.oss << "PE horizon user setting: " << horizon;
    //    } else {
    //        LOG.oss << "Short-range paired-end horizon: " << PEhorizon << endl;
    //        LOG.oss << "Long-range paired-end horizon:  " << MPhorizon << endl;
    //
    //    }

    bool PE = false, MP = false;
    for (unsigned int i = 0; i < P->getNLibrary(); i++) {
        if (P->getLibraryIt(i)->getMateOrientation() == 1)
            PE = true;
        if (P->getLibraryIt(i)->getMateOrientation() == 2)
            MP = true;
    }
    if (PE)
        LOG.oss << "Short-range paired-end horizon: " << PEhorizon << endl;
    if (MP)
        LOG.oss << "Long-range paired-end horizon:  " << MPhorizon << endl;
    LOG.flushStream(TOSTDOUT);

    vector<unsigned int> distr;

    size_t start, end;
    ofstream out;
    ostringstream oss;
    string outFile;

    for (unsigned int lib = 0; lib < P->getNLibrary(); lib++) {
        meanLength = 0.0;
        sdLength = 0.0;
        oss.clear();
        oss.str("");
        oss << prefix;
        oss << "_peDistSample_";
        oss << (lib + 1);
        outFile = oss.str();
        out.open(outFile.c_str());

        start = P->getLibraryIt(lib)->getStartIndex();
        end = P->getLibraryIt(lib)->getEndIndex();
        nCounted = 0;
        nUnique = 0;
        nMultiple = 0;
        nCanceled = 0;
        nZero = 0;
        usableMateEstimate = 0;

        unsigned int horizon = PEhorizon;

        if (P->getLibraryIt(lib)->getMateOrientation() == 1) {
            mate1Dir = true;
            mate2Dir = false;
            horizon = PEhorizon;
        } else if (P->getLibraryIt(lib)->getMateOrientation() == 2) {
            mate1Dir = false;
            mate2Dir = true;
            horizon = MPhorizon;
        } else {
            cerr << "[bug] OverlapsGraph::estimatePairedDistance(...)\n";
            sendBugReportPlease(cerr);
        }

        distr.assign(horizon, 0);


        for (unsigned int i = start; i <= end; i++) {
            //for all pairing

            if (globalCounter % 1000 == 0) {
                cerr << "Estimating pairing distances..." << (int) (((float) globalCounter / P->getNPairing())*100) << "%\r" << flush;
            }
            globalCounter++;
            //get IDs of the paired reads
            id1 = P->getR1(i);
            id2 = P->getR2(i);

            retValue = getPathBetweenTwoReads(id2, mate1Dir, id1, mate2Dir, 0, horizon, maxLogNLeaves);

            if (retValue == -1) {
                nCanceled++;
            } else if (retValue == 1) {
                //out << *(Node::VP.getPath(0)) << " ";
                //   out << ">" << nUnique << "_1\n" << L->getDirectRead(id2) << '\n';
                //     out << ">" << nUnique << "_2\n" << L->getDirectRead(id1) << '\n';
                nUnique++;
                pathLength = (size_t) Node::pathLengthMean;

                distr[pathLength]++;
                meanLength += pathLength;
                sdLength += pathLength*pathLength;

                out << pathLength << ' ';
                usableMateEstimate++;
            } else if (retValue > 1) {
                nMultiple++;
                //
                //                if (Node::pathLengthSD < 5)
                //                {
                //                    meanLength += Node::pathLengthMean;
                //                    sdLength += Node::pathLengthMean * Node::pathLengthMean;
                //                    nCounted++;
                //                }
                usableMateEstimate++;
            } else
                nZero++;
        }

        out.close();

        nCounted += nUnique;

        meanLength /= nCounted;
        sdLength = sqrt((1.0 / (nCounted - 1))*(sdLength - nCounted * meanLength * meanLength));

        P->getLibraryIt(lib)->setMeanCloneSize(meanLength);
        P->getLibraryIt(lib)->setSDCloneSize(sdLength);
        P->getLibraryIt(lib)->setNEchant(nCounted);

        int min = (int) (meanLength - (nsd * sdLength));
        int max = (int) (meanLength + (nsd * sdLength));
        if (min < 0)
            min = 0;

        if (nCounted < 20)
            min = max = 0;

        //bounded sample (mean +- nsd*sd)
        for (int i = 0; i < max - min + 1; i++)
            distr[i] = distr[min + i];

        distr.resize(max - min + 1, 0);
        P->getLibraryIt(lib)->lengthDistr.setDistribution(distr, min, max);
        //  P->getLibraryIt(lib)->setVecDistr(distr);

        usableMateEstimate += (nCanceled * ((double) usableMateEstimate / (usableMateEstimate + nZero)));
        P->getLibraryIt(lib)->setNUsableMates(usableMateEstimate);
        P->getLibraryIt(lib)->setExpectedMateCoverage((double) usableMateEstimate / targetSizeEstimation);
    }
    cout << "Estimating pairing distance... done               \n";
    P->updatePERange(nsd);
    PELibrary::readLength = getReadLength();
}

void OverlapsGraph::testChimera(unsigned int maxLength) {

    int retValue;
    unsigned int id1, id2;

    unsigned int nUnique = 0, nMultiple = 0, nCancelled = 0, nZero = 0;
    double maxLogNLeaves = 6.0;

    for (unsigned int i = 0; i < P->getNPairing(); i++) {
        //for all pairing

        if (i % 1000 == 0) {
            cout << "Testing for chimera..." << (int) (((float) i / P->getNPairing())*100) << "%\r" << flush;
        }

        //get IDs of the paired reads
        id1 = P->getR1(i);
        id2 = P->getR2(i);

        retValue = getPathBetweenTwoReads(id2, true, id1, true, 0, maxLength, maxLogNLeaves);

        if (retValue == -1) {
            nCancelled++;
        } else if (retValue == 1) {
            nUnique++;
        } else if (retValue > 1) {
            nMultiple++;
        } else {
            nZero++;
        }

        retValue = getPathBetweenTwoReads(id2, false, id1, false, 0, maxLength, maxLogNLeaves);

        if (retValue == -1) {
            nCancelled++;
        } else if (retValue == 1) {
            nUnique++;
        } else if (retValue > 1) {
            nMultiple++;
        } else {
            nZero++;
        }
    }

    cout << "Testing for chimera... done                \n";
    cout << "   Total pairing  " << P->getNPairing() << endl;
    cout << "   Unique bad path    " << nUnique << " (" << ((float) nUnique / P->getNPairing())*100 << "%)\n";
    cout << "   Multiple bad paths " << nMultiple << " (" << ((float) nMultiple / P->getNPairing())*100 << "%)\n";
    cout << "   No paths       " << nZero << " (" << ((float) nZero / P->getNPairing())*100 << "%)\n";
    cout << "   Canceled      " << nCancelled << " (" << ((float) nCancelled / P->getNPairing())*100 << "%)\n";
}

void OverlapsGraph::sortNodesByLength() {
    rank.clear();
    nodeRank myRank;
    for (unsigned int i = 1; i <= getNNodes(); i++) {
        myRank.nodeIndex = i;
        rank.push_back(myRank);
    }

    cerr << "Sorting nodes..." << flush;
    sort(rank.begin(), rank.end());
    cerr << "done" << endl;
}

void OverlapsGraph::assemble(
        string prefix,
        unsigned int minContigSize,
        double minCoverage,
        int trim,
        int minNPair,
        double minRatio,
        unsigned int maxRedundancy) {

    AP.init(); //assembly progress status
    unsigned int totSize = 0;
    sortNodesByLength();

    if (maxRedundancy > 0) {
        LOG.oss << "[warn] maxRedundancy: " << maxRedundancy << '\n';
        LOG.oss << "[warn] Assembly may contain redundant information and be artificially \n";
        LOG.oss << "[warn] increased in size\n";
        LOG.flushStream(TOSTDOUT);
    }
    ofstream out;
    ofstream out_cov;
    ofstream out_lay;
    ofstream out_tab;
    ofstream out_nodeInfo;
    out.open((prefix + "_contigs.fasta").c_str());
    out_cov.open((prefix + "_contigs.cov").c_str());

    // if (DEV_INFO)
    // {
    out_lay.open((prefix + "_contigs.lay").c_str());
    out_tab.open((prefix + "_nodesPosition").c_str());
    out_nodeInfo.open((prefix + "_nodesInfo").c_str());
    // }

    unsigned int contigCount = 0;
    string contig;
    bool endingLeft, endingRight;
    unsigned int nBranching = 0, nNoOv = 0;
    ostringstream oss;
    unsigned int nReadInContig;
    int breakCauseRight = 0, breakCauseLeft = 0;
    int nCollision = 0, nDE = 0, nAmb = 0, nCanceled = 0, nBubbles = 0;

    vector<unsigned int> mainPath, nodePath, tmpPath, coverage;
    vector< vector<unsigned int> > assembly;
    vector<bool> nodePathDir;
    vector<int> echant;
    Node::S_counter = 0;


    //flagB: has been traversed

    cout << "Building contigs...\n" << flush;

    for (unsigned int i = 1; i <= getNNodes(); i++)
        nodesTab[i].unsetAlreadyUsed();
    //   nodesTab[i].unsetFlagB();

    for (unsigned int node = 1; node <= getNNodes(); node++) {

        // if (DEV_INFO)
        // {
        out_nodeInfo << "node" << node << " "
                << nodesTab[node].getSequenceLength() << " "
                << nodesTab[node].getCoverage() << '\n';
        // }

        unsigned int i = rank.at(node - 1).nodeIndex;

        if (nodesTab[i].isAlreadyUsed() == true)
            continue;

        nodesTab[i].isAlreadyUsed(); // traversed
        Node::S_counter++;

        //1 collision
        //2 search tree reaches a dead-end
        //3 out of paired-end range
        //4 canceled, step limit reached or max tree size reached
        //5 non resolvable bubble

        breakCauseRight = PEelongate(tmpPath, i, true, minNPair, minRatio, maxRedundancy);

        if (P->getNPairing() != 0) {
            breakCauseLeft = PEelongate(mainPath, i, false, minNPair, minRatio, maxRedundancy);
            reverseEPath(mainPath);
            mainPath.insert(mainPath.end(), tmpPath.begin() + 2, tmpPath.end());
        } else {
            mainPath = tmpPath;
        }

        //  assembly.push_back(mainPath);

        ePathToNodePath(mainPath, nodePath, nodePathDir);

        if (P->getNPairing() == 0) {
            //unpaired data: only two possible elongation break causes
            if (nodesTab[nodePath.at(1)].isEnding(!nodePathDir.at(1)))
                breakCauseLeft = 2; //dead end
            else
                breakCauseLeft = 3; //branching 

            if (nodesTab[*(nodePath.end() - 1)].isEnding(*(nodePathDir.end() - 1)))
                breakCauseRight = 2; //dead end
            else
                breakCauseRight = 3; //branching
        }

        //not used anymore
        endingLeft = nodesTab[nodePath.at(1)].isEnding(!nodePathDir.at(1));
        endingRight = nodesTab[*(nodePath.end() - 1)].isEnding(*(nodePathDir.end() - 1));

        //Basic rule to determine the circularity of a contig
        //may miss some cases, to improve.
        bool potentialyCircular = false;
        bool circular = false;
        NodeIt nIt, nIt2, nextNIt;

        nIt.initNodeIt(*(nodePath.end() - 1), *(nodePathDir.end() - 1));
        nIt2.initNodeIt(*(nodePath.begin() + 1), *(nodePathDir.begin() + 1));

        if (nIt.getNNeighbor() == 1) {
            nIt.getNext(nextNIt);
            if (nextNIt.getNodeId() == nodePath.at(1) && nextNIt.getDirection() == nodePathDir.at(1)) {
                circular = true;
            }
        } else {
            nIt.initIterator();
            while (nIt.getNext(nextNIt)) {
                if (nextNIt == nIt2) {
                    potentialyCircular = true;
                    break;
                }
            }
        }

        oss.clear();
        oss.str("");

        //return value: 1=collision, 2=dead-end, 3=ambiguity, 4=canceled
        switch (breakCauseLeft) {
            case 1: //out of paired-end range in a collision context
                oss << ':';
                nCollision++;
                break;
            case 2: //dead-end
                oss << '|';
                nDE++;
                break;
            case 3: //out of paired-end range
                oss << '~';
                nAmb++;
                break;
            case 4: //search canceled
                oss << 'X';
                nCanceled++;
                break;
            case 5: //non resolved bubble
                oss << 'B';
                nBubbles++;
                break;
            default:
                cout << "[bug] This should not happen\n";
                sendBugReportPlease(cerr);
        }
        switch (breakCauseRight) {
            case 1: //out of paired-end range in a collision context
                oss << ':';
                nCollision++;
                break;
            case 2: //dead-end
                oss << '|';
                nDE++;
                break;
            case 3: //out of paired-end range
                oss << '~';
                nAmb++;
                break;
            case 4: //search canceled
                oss << 'X';
                nCanceled++;
                break;
            case 5: //non resolved bubble
                oss << 'B';
                nBubbles++;
                break;
            default:
                cout << "[bug] This should not happen\n";
                sendBugReportPlease(cerr);
        }

        oss << ' ';

        //        if (endingLeft)
        //            oss << '[';
        //        if (endingRight)
        //            oss << ']';

        if (potentialyCircular)
            oss << "==";
        else if (circular)
            oss << "OO";

        oss << "\n";

        for (size_t ii = 1; ii < nodePath.size(); ii++) {
            if (!nodePathDir[ii])
                oss << '!';
            oss << nodePath.at(ii);

            if (ii < nodePath.size() - 1)
                oss << '-';
        }

        contig = ePathToSeq(mainPath);
        ePathToCov(mainPath, coverage);
        nReadInContig = getNReadsInEPath(mainPath);

        double contigCoverage = (double) nReadInContig * R->getReadsLength() / (contig.length() - R->getReadsLength() + 1);
        if (nReadInContig < 3)
            contigCoverage = 0.0;

        if (contigCoverage < minCoverage)
            continue;

        int tLeft = 0, tRight = contig.size() - 1;
        if (trim > 1) {
            while (tLeft < (int) coverage.size()) {
                if ((int) coverage[tLeft] >= trim)
                    break;
                tLeft++;
            }

            while (tRight >= 0) {
                if ((int) coverage[tRight] >= trim)
                    break;
                tRight--;
            }

            if (tLeft <= tRight) {
                contig = contig.substr(tLeft, tRight - tLeft + 1);
            } else {
                contig = "";
                coverage.clear();
            }
        }

        if (contig.length() >= minContigSize) {
            echant.push_back(contig.length());
            contigCount++;
            totSize += contig.length();

            AP.setNContigs(contigCount);
            AP.setKb(totSize);
            AP.printProgress(cerr, true);

            out_lay << ">" << prefix << "_" << contigCount;
            out_lay << " size=" << contig.length();
            out_lay << " cov=" << contigCoverage;
            out_lay << " node:" << i << " ";
            out_lay << oss.str();
            out_lay << '\n';
            out_lay << flush;
            out_tab << ">" << prefix << "_" << contigCount;
            out_tab << " size=" << contig.length();
            out_tab << " cov=" << contigCoverage << "\n";
            ePathTabularInfo(mainPath, out_tab);

            oss.clear();
            oss.str("");

            oss << ">" << prefix << "_" << contigCount;
            oss << " size=" << contig.length();
            oss << " cov=" << contigCoverage;
            if (potentialyCircular)
                oss << " possibly_circular";
            if (circular)
                oss << " circular";
            oss << '\n';

            out_cov << oss.str();

            for (int ii = tLeft; ii <= tRight; ii++)
                out_cov << coverage.at(ii) << ' ';
            out_cov << '\n';

            out << oss.str();

            lineWrap(out, contig, 70);

            if (endingRight)
                nNoOv++;
            else
                nBranching++;
            //            else if (ambRight)
            //                nBranching++;

            if (endingLeft)
                nNoOv++;
            else
                nBranching++;
            //            else if (ambLeft)
            //                nBranching++;
        }
    }

    // cerr << endl;

    cout << "\rBuilding contigs... done                               " << endl;

    LOG.oss << "Number of contigs:  " << contigCount << endl;
    stats(echant.begin(), echant.end(), LOG.oss);
    LOG.oss << "Assembly breaks occurred due to:\n";

    // if (P->getNPairing() == 0) {
    LOG.oss << "   non-resolved ambiguity: " << nBranching << "\n";
    LOG.oss << "   dead-end: " << nNoOv << "\n";
    //    }
    //    else
    //    {
    //        LOG.oss << " ~  ambiguity: out of paired end range: " << nAmb << endl;
    //        LOG.oss << " B  ambiguity: unresolved bubble(s): " << nBubbles << endl;
    //        LOG.oss << " :  ambiguity: lost in redundancy (out of paired-end range): " << nCollision << endl;
    //        LOG.oss << " |  dead-end: " << nDE << endl;
    //        LOG.oss << " X  canceled search: " << nCanceled << endl;
    //    }
    LOG.flushStream(TOSTDOUT);

    out.close();
    out_cov.close();



    //   if (DEV_INFO)
    //   {  
    out_lay.close();
    out_tab.close();
    out_nodeInfo.close();
    //   }
}

int OverlapsGraph::PEelongate(vector<unsigned int>& path,
        unsigned int node,
        bool direction,
        int minNPair,
        double minRatio,
        unsigned int maxRedundancy) {

    if (P->getNPairing() == 0 || nodesTab[node].isEnding(direction)) {
        path.clear();
        //path.push_back(nodesTab[node].getSequenceLength());
        if (direction)
            path.push_back(1); //first node is direct
        else
            path.push_back(0);
        path.push_back(node);

        if (nodesTab[node].isEnding(direction))
            return 2; //dead-end
        else
            return 3; //unresolved ambiguity
    }

    static BackwardsWalker BK;
    if (!BK.isAllocated())
        BK.init(this);

    BK.initWalker(node, direction, minNPair, minRatio, maxRedundancy, 500);

    return BK.elongate(path);
}

string OverlapsGraph::ePathToSeq(const vector<unsigned int> &path) {

    //A path is a vector v where:
    //v[0] : first node direction: 1 direct, 0 reverse, used when v.size()==2)
    //v[1] : starting node
    //v[2] .. v[n] : edges index forming the path.

    string s;
    string tmp;
    bool checkDir = true; //to check path validity

    if (path.size() < 2)
        return "";
    else if (path.size() == 2) {
        if (path[0] == 1)
            return nodesTab[path.at(1)].getDirectSequence();
        else if (path[0] == 0)
            return nodesTab[path.at(1)].getReverseSequence();
        else {
            cout << "OverlapsGraph::ePathToSeq(...) problem\n";
            sendBugReportPlease(cout);
        }
    }

    unsigned int id, size;
    bool dir, currentDir;
    unsigned int currentNode;

    currentNode = path.at(1);
    if (path.at(2) < nodesTab[currentNode].getNRightOverlap())
        currentDir = true;
    else
        currentDir = false;

    s = nodesTab[currentNode].getSequence(currentDir);

    for (unsigned int i = 2; i < path.size(); i++) {
        //check path validity
        if (path.at(i) < nodesTab[currentNode].getNRightOverlap())
            checkDir = true;
        else
            checkDir = false;

        if (checkDir != currentDir) {
            return "";
        }

        nodesTab[currentNode].getOverlap(path.at(i), id, size, dir);
        if (!dir)
            currentDir = !currentDir;
        tmp = nodesTab[id].getSequence(currentDir);
        tmp = tmp.substr(size);
        s += tmp;
        //s+=nodesTab[id].getSequence(dir).substr(size);
        currentNode = id;
    }
    return s;
}

unsigned int OverlapsGraph::getNReadsInEPath(const vector<unsigned int> &path) {
    if (path.size() < 2)
        return 0;
    else if (path.size() == 2)
        return nodesTab[path.at(1)].getNReads();

    unsigned int id, size;
    bool dir, currentDir, checkDir;
    unsigned int currentNode;
    unsigned int nReads = 0;

    currentNode = path.at(1);
    if (path.at(2) < nodesTab[currentNode].getNRightOverlap())
        currentDir = true;
    else
        currentDir = false;

    nReads += nodesTab[currentNode].getNReads();

    for (unsigned int i = 2; i < path.size(); i++) {
        //check path validity
        if (path.at(i) < nodesTab[currentNode].getNRightOverlap())
            checkDir = true;
        else
            checkDir = false;

        if (checkDir != currentDir) {
            return 0;
        }

        nodesTab[currentNode].getOverlap(path.at(i), id, size, dir);
        if (!dir)
            currentDir = !currentDir;

        nReads += nodesTab[id].getNReads();
        currentNode = id;
    }
    return nReads;
}

void OverlapsGraph::nodeListToEPath(string args, vector<unsigned int> &ePath)//for DEV mode
{
    //return an empty ePath is the stringPath is not valid in G

    ePath.clear();
    vector<string> nodesStr;
    vector<unsigned int> nodeList;
    vector<bool> nodeDir;
    deque<unsigned int> nl;
    deque<bool> nd;
    unsigned int value;
    string buffer;
    istringstream iss;
    bool reversePath = false;

    if (args[0] == 'R') {
        reversePath = true;
        args = args.substr(1);
    }

    iss.str(args);

    do {
        getline(iss, buffer, '-');
        if (!iss.fail()) {
            nodesStr.push_back(buffer);
        }

    } while (!iss.eof());

    for (size_t i = 0; i < nodesStr.size(); i++) {
        if (nodesStr.at(i)[0] == '!') {
            nodeDir.push_back(false);
            nodesStr.at(i) = nodesStr.at(i).substr(1);
        } else
            nodeDir.push_back(true);

        iss.clear();
        iss.str(nodesStr.at(i));
        iss >> value;
        if (iss.fail())
            return;
        nodeList.push_back(value);
    }

    if (reversePath == true) {
        for (size_t i = 0; i < nodesStr.size(); i++) {
            nl.push_front(nodeList[i]);
            nd.push_front(!nodeDir[i]);
        }
        for (size_t i = 0; i < nodesStr.size(); i++) {
            nodeList[i] = nl[i];
            nodeDir[i] = nd[i];
        }
    }

    if (nodeList.size() == 0)
        return;

    nodePathToEPath(nodeList, nodeDir, ePath);
}

void OverlapsGraph::ePathToNodePath(
        const vector<unsigned int> &path,
        vector<unsigned int> &nodePath,
        vector<bool> &vDir) {
    unsigned int id, size;
    bool dir, currentDir;
    unsigned int edgeIndex;
    nodePath.clear();
    vDir.clear();

    if (path.size() < 2)
        return;

    if (path.size() == 2) {

        nodePath.push_back(nodesTab[path[1]].getSequenceLength());
        nodePath.push_back(path.at(1));
        vDir.push_back(true); // not used

        if (path.at(0) == 1)
            vDir.push_back(true);
        else
            vDir.push_back(false);
        return;
    }

    nodePath.push_back(0);
    vDir.push_back(true); // not used

    unsigned int currentNode = path.at(1);
    unsigned int pathSize = nodesTab[currentNode].getSequenceLength();

    if (path.at(2) < nodesTab[currentNode].getNRightOverlap())
        currentDir = true;
    else
        currentDir = false;

    nodePath.push_back(currentNode);
    vDir.push_back(currentDir);

    for (size_t i = 2; i < path.size(); i++) {
        edgeIndex = path.at(i);
        nodesTab[currentNode].getOverlap(edgeIndex, id, size, dir);
        if (!dir)
            currentDir = !currentDir;
        currentNode = id;
        nodePath.push_back(currentNode);
        vDir.push_back(currentDir);
        pathSize += nodesTab[currentNode].getSequenceLength() - size;
    }
    nodePath.at(0) = pathSize;
}

void OverlapsGraph::nodePathToEPath(const vector<unsigned int> &nodePath, const vector<bool> &vDir, vector<unsigned int> &path) {
    path.clear();

    if (nodePath.size() != vDir.size()) {
        cout << "OverlapsGraph::nodePathToPath(...) problem\n";
        sendBugReportPlease(cout);
    }

    if (nodePath.size() == 0)
        return;

    //store the starting direction
    if (vDir[0] == true)
        path.push_back(1);
    else
        path.push_back(0);

    path.push_back(nodePath.at(0));

    bool found = false;
    unsigned int memIndex;
    NodeIt nIt1, nIt2, nItNext;

    for (size_t n = 0; n < nodePath.size() - 1; n++) {
        nIt1.initNodeIt(nodePath.at(n), vDir.at(n));
        nIt2.initNodeIt(nodePath.at(n + 1), vDir.at(n + 1));

        found = false;
        while (nIt1.getNext(nItNext) != 0) { //check for multiple edges between two nodes (rare)
            if (nItNext == nIt2) {
                if (found == true)//ambiguous seq
                {
                    cerr << "Ambiguous branching from node: " << nodePath.at(n) << " to node: " << nodePath.at(n + 1) << endl;
                    path.clear();
                    return;
                }
                found = true;
                memIndex = nItNext.getAbsArrivalEdgeIndex();
            }
        }
        if (!found) {
            cerr << "Invalid path: " << nodePath.at(n) << " to node: " << nodePath.at(n + 1) << endl;
            path.clear();
            return;
        }

        path.push_back(memIndex);
    }

}

void OverlapsGraph::reverseEPath(vector<unsigned int> & path) {
    if (path.size() < 2)
        return;
    vector<unsigned int> reverse;

    unsigned int currentId = path[1];
    unsigned int reverseEdgeIndex;

    for (size_t i = 2; i < path.size(); i++) {
        reverseEdgeIndex = nodesTab[currentId].getReciprocal(path[i]);
        reverse.push_back(reverseEdgeIndex);
        currentId = nodesTab[currentId].getOverlapId(path[i]);
    }
    reverse.push_back(currentId);

    //path[0] stores the direction of the first node 1 direct, 0 reverse
    //this is used only when path is made up of a single node.
    if (path[0] == 1)
        reverse.push_back(0);
    else
        reverse.push_back(1);

    path.clear();
    for (size_t i = reverse.size(); i > 0; i--) {
        path.push_back(reverse.at(i - 1));
    }
}

void OverlapsGraph::ePathToCov(const vector<unsigned int> &path, vector<unsigned int> &cov) {
    cov.clear();
    vector<unsigned int> covTmp;

    if (path.size() < 2)
        return;

    if (path.size() == 2) {
        nodesTab[path.at(1)].getVCoverage(cov, true);
        return;
    }

    unsigned int id, size;
    bool dir, currentDir;
    unsigned int edgeIndex;
    size_t pos;

    unsigned int currentNode = path.at(1);

    if (path.at(2) < nodesTab[currentNode].getNRightOverlap())
        currentDir = true;
    else
        currentDir = false;

    nodesTab[currentNode].getVCoverage(cov, currentDir);

    for (size_t i = 2; i < path.size(); i++) {
        edgeIndex = path.at(i);
        nodesTab[currentNode].getOverlap(edgeIndex, id, size, dir);
        if (!dir)
            currentDir = !currentDir;
        currentNode = id;

        nodesTab[id].getVCoverage(covTmp, currentDir);

        //add covTmp to cov
        pos = cov.size() - size;
        for (size_t ii = 0; ii < size; ++ii) {
            cov[pos] += covTmp[ii];
            pos++;
        }
        cov.insert(cov.end(), covTmp.begin() + size, covTmp.end());
    }
}

void OverlapsGraph::ePathTabularInfo(const vector<unsigned int> &path, ostream &out) {

    unsigned int id, size;
    bool dir, currentDir;
    unsigned int currentNode = path.at(1);
    unsigned int start = 1;

    if (path.size() < 2)
        return;

    if (path.size() == 2) {
        //assume "direct" direction for the node
        currentNode = path.at(1);
        currentDir = true;
    } else if (path.at(2) < nodesTab[currentNode].getNRightOverlap())
        currentDir = true;
    else
        currentDir = false;

    if (!currentDir)
        out << '!';
    out << currentNode << '\t';
    out << start << '\t';
    out << start + nodesTab[currentNode].getSequenceLength() - 1 << '\t';
    out << "l=" << nodesTab[currentNode].getSequenceLength() << '\t';
    out << "cov=" << nodesTab[currentNode].getCoverage() << '\n';
    start += nodesTab[currentNode].getSequenceLength();

    for (unsigned int i = 2; i < path.size(); i++) {
        nodesTab[currentNode].getOverlap(path.at(i), id, size, dir);
        start -= size;
        if (!dir)
            currentDir = !currentDir;
        if (!currentDir)
            out << '!';
        currentNode = id;

        out << currentNode << '\t';
        out << start << '\t';
        out << start + nodesTab[currentNode].getSequenceLength() - 1 << '\t';
        out << "l=" << nodesTab[currentNode].getSequenceLength() << '\t';
        out << "cov=" << nodesTab[currentNode].getCoverage() << '\n';
        start += nodesTab[currentNode].getSequenceLength();
    }
}

void OverlapsGraph::estimateCoverage(double &minCoverage, double &targetSize) {

    LOG.oss << "Nodes coverage sampling:" << endl;

    double a = 0.0, previousA;
    double q = 0.0;
    double sLength;
    double cover;
    unsigned int nReads;
    double sumWeight = 0.0;
    double min = 1E30, max = 0.0;
    unsigned int nSample = 0;
    unsigned int rL = getReadLength();

    vector<orderCov> vCov;
    vCov.reserve(getNNodes());
    orderCov oc;

    //ofstream out("covStats");
    vector<int> samples;
    vector<float> weights;

    for (unsigned int i = 1; i <= getNNodes(); i++) {

        nReads = nodesTab[i].getNReads();
        cover = nodesTab[i].getCoverage();
        sLength = nodesTab[i].getSequenceLength() - rL + 1;

        if (cover == 0)
            continue;
        if (nReads < 10)
            continue;

        //out << cover << '\t' << sLength << '\n';
        samples.push_back(cover);
        weights.push_back(sLength);

        oc.cov = cover;
        oc.weight = sLength;
        vCov.push_back(oc);

        if (cover > max)
            max = cover;
        if (cover < min)
            min = cover;
        nSample++;
        sumWeight += sLength;
        previousA = a;
        a = previousA + (sLength / sumWeight)*(cover - previousA);
        q = q + sLength * (cover - previousA) * (cover - a);
    }
    // out.close();

    double sig = q / sumWeight;
    double sd;

    if (nSample < 2)
        sd = 0;
    else
        sd = sqrt(nSample / (nSample - 1) * sig);

    sort(vCov.begin(), vCov.end());
    double sum = 0.0;
    double median;
    double q01 = 0.0;
    double q10 = 0.0;

    for (size_t i = 0; i < vCov.size(); i++) {
        if (sum + vCov.at(i).weight < sumWeight * .01) {
            q01 = vCov.at(i).cov;
        }
        if (sum + vCov.at(i).weight < sumWeight * .10) {
            q10 = vCov.at(i).cov;
        }

        if (sum + vCov.at(i).weight > sumWeight * .5) {
            median = vCov.at(i).cov;
            break;
        }
        sum += vCov.at(i).weight;
    }

    LOG.oss << "   mean: " << a << '\n';
    LOG.oss << "   median: " << median << '\n';
    LOG.oss << "   1st percentile " << q01 << '\n';
    LOG.oss << "   1st decile " << q10 << '\n';
    LOG.oss << "   sd: " << sd << '\n';
    // cout << "   target size estimation: " << sumWeight << endl;

    LOG.oss << "   minimum required coverage";

    if (minCoverage == 0.0) {
        //maybe not conservative enough...
        minCoverage = a / 2.0;

        //far more conservative
        //        minCoverage = median/2.0;
        //        
        //        if (minCoverage > q10)
        //            minCoverage = q10;

        LOG.oss << " automatically set to: " << minCoverage << endl;
    } else {
        LOG.oss << ": " << minCoverage << " (user setting)" << endl;
    }

    unsigned int totLength = 0;
    unsigned int R;
    for (unsigned int i = 1; i <= getNNodes(); i++) {

        cover = nodesTab[i].getCoverage();
        if (cover < minCoverage)
            continue;

        sLength = nodesTab[i].getSequenceLength() - rL + 1;
        R = (int) cover / a;
        if (R == 0)
            R = 1;
        totLength += R*sLength;
    }

    if (targetSize == 0) {
        targetSize = totLength;
        LOG.oss << "   target size (roughly) estimated as " << totLength << endl;
    } else {
        LOG.oss << "   target size set by the user as " << targetSize;
    }

    LOG.flushStream(TOSTDOUT);
    targetSizeEstimation = targetSize;

}
