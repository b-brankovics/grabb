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
#include "Pairing.h"
#include "readsLayout.h"
#include "globalFunc.h"
#include "DevShell.h"
#include "stat.h"
#include "BeamSearchTree.h"
#include "BackwardsWalker.h"
#include "logWriter.h"
#include "Param.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <limits>

using namespace std;

bool DEV_INFO;

ofstream outGlob; //for dev purpose
unsigned long global_count;
logWriter LOG;

int main(int argc, char **argv) {

    //  outGlob.open("DEVLOG");
    //  outGlob << "node\tindex\tsize\tprob\n";

    //parse and check command line
    Param param;
    param.parseCommandLine(argc, argv);

    global_count = 0;

    float CCreliableCutoff = 1e-2;
    float CCsupectCutoff = 1e-4;
    float CCGincoherantCutoff = 1e-6;
    float nsd = 2.0; //used to determine min-max allowed distances for PE libraries

    Pairing *P = new Pairing();
    ReadsLayout *L = new ReadsLayout();
    ReadsStorage *R = new ReadsStorage();
    P->setRpointer(R);
    OverlapsGraph G;
    G.init(R, P, L);

    cout << setprecision(1) << fixed;

    //Overlapping mode
    if (param.overlappingMode()) {
        int maxReadLength = numeric_limits<short int>::max();

        string logFile = param.prefix + "_overlapping.log";
        if (!LOG.open(logFile, true)) {
            cerr << "[err] Failed opening log file" << endl;
        }

        LOG.oss << edenaVersion();
        LOG.flushStream(TOSTDOUT);
        for (int i = 1; i < argc; i++) {
            LOG.oss << "argv[" << i << "] " << argv[i] << '\n';
        }
        LOG.flushStream();

        //fast estimation of the total number of reads
        //for the Layout initiation
        cout << "Rapid file(s) examination... " << flush;

        unsigned int nReadsEstimation = 0; //upper bound
        unsigned int nL = 0;
        int rl;
        int minRl = numeric_limits<int>::max();
        unsigned int ret;

        vector<string> allReadsFile;
        allReadsFile.insert(allReadsFile.end(),
                param.drPairedEndFiles.begin(),
                param.drPairedEndFiles.end());
        allReadsFile.insert(allReadsFile.end(),
                param.rdPairedEndFiles.begin(),
                param.rdPairedEndFiles.end());
        allReadsFile.insert(allReadsFile.end(),
                param.singleEndFiles.begin(),
                param.singleEndFiles.end());

        for (size_t i = 0; i < allReadsFile.size(); i++) {
            ret = estimateNReads(allReadsFile[i], rl);
            if (ret == 0)
                exit(0);
            nReadsEstimation += ret;
            if (rl != minRl) {
                nL++;
                if (rl < minRl) {
                    minRl = rl;
                }
            }
        }

        cout << "done" << endl;
        LOG.oss << "Number of reads upper bound estimation: " << nReadsEstimation << '\n';
        LOG.flushStream(TOSTDOUT);

        if (minRl > maxReadLength) {
            LOG.oss << "   The maximum supported read length is " << maxReadLength << "bp\n";
            LOG.oss << "   All reads have been truncated accordingly" << endl;
            LOG.flushStream(TOSTDOUT);
            param.truncateTo = maxReadLength;
            minRl = maxReadLength;
        }

        if (param.truncateTo == 0) {
            if (nL > 1) {
                LOG.oss << "The shorter reads are " << minRl << " bases in length\n";
                LOG.oss << "Other reads will be truncated accordingly\n";
                LOG.flushStream(TOSTDOUT);
            }
            param.truncateTo = minRl;
        }

        if (param.truncateTo > minRl) {
            LOG.oss << "Adjusting reads truncation to " << minRl << endl;
            LOG.flushStream(TOSTDOUT);
            param.truncateTo = minRl;
        }

        cout << "Reads layout initialization... ";
        nReadsEstimation += nReadsEstimation * .1; //10% more
        L->init(R, nReadsEstimation);
        R->init(nReadsEstimation, param.truncateTo);
        cout << "done" << endl;

        //Load reads files
        int mateOrientation = 1; // DR

        for (size_t i = 0; i < param.singleEndFiles.size(); i++) {
            if (
                    R->loadReadsFiles(param.singleEndFiles[i], "", 0, P, L)
                    != 0)
                return EXIT_FAILURE;
        }

        for (size_t i = 0; i < param.drPairedEndFiles.size(); i += 2) {
            if (
                    R->loadReadsFiles(param.drPairedEndFiles[i],
                    param.drPairedEndFiles[i + 1],
                    mateOrientation,
                    P,
                    L) != 0)
                return EXIT_FAILURE;
        }

        mateOrientation = 2; //RD
        for (size_t i = 0; i < param.rdPairedEndFiles.size(); i += 2) {
            if (
                    R->loadReadsFiles(param.rdPairedEndFiles[i],
                    param.rdPairedEndFiles[i + 1],
                    mateOrientation,
                    P,
                    L) != 0)
                return EXIT_FAILURE;
        }

        R->adjustAllocation();

        cout << "done\n";
        LOG.oss << "   Number of reads: " << R->getEffectiveNReads() << endl;
        LOG.oss << "   Number of distinct sequences: " << R->getN_nrReads() << endl;
        LOG.oss << "   Average reads redundancy: " << (float) R->getEffectiveNReads() / R->getN_nrReads() << endl;
        LOG.flushStream(TOSTDOUT);


        if (param.drPairedEndFiles.size() >= 2 || param.rdPairedEndFiles.size() >= 2) {
            cout << "Building pairing index... " << flush;
            P->buildIndex(R->getEffectiveNReads());
            cout << "done" << endl;
        }

        cout << "Initializing prefix tables..." << flush;
        LOG.oss << "Initializing prefix tables";
        LOG.flushStream();
        R->adjustAllocation();
        R->initPrefixTables();
        cout << "done" << endl;

        if (param.minOverlap == 0)
            param.minOverlap = R->getReadsLength() / 2;
        R->setMinOvSize(param.minOverlap);

        LOG.oss << "Minimum overlaps size: " << param.minOverlap << endl;
        LOG.flushStream(TOSTDOUT);

        LOG.oss << "Start computing overlaps";
        LOG.flushStream();
        G.computeOverlaps(param.nThreads);
        G.removeTransitiveEdges();  
        LOG.oss << "Done computing overlaps";
        LOG.flushStream();
        G.countEdges();
        LOG.oss << "Number of nodes: " << G.getNNodes() << endl;
        LOG.flushStream();
        LOG.oss << "Number of edges: " << G.getNEdges() << endl;
        LOG.flushStream();

        if (param.checkGraph)
            G.checkConsistency();

        if (!G.saveData((param.prefix + ".ovl").c_str())) {
            cerr << "[err] failed writing the .ovl file\n";
            return EXIT_FAILURE;
        }
        LOG.oss << "OSG written to file " << param.prefix + ".ovl";
        LOG.flushStream();

        return EXIT_SUCCESS;
        //End overlapping mode
    }// Assembling mode
    else //ASSEMBLY MODE
    {

        string logFile = param.prefix + "_assembling.log";
        if (!LOG.open(logFile, true)) {
            cerr << "[err] Failed opening log file" << endl;
        }

        LOG.oss << edenaVersion();
        LOG.flushStream(TOSTDOUT);
        for (int i = 1; i < argc; i++) {
            LOG.oss << "argv[" << i << "] " << argv[i] << '\n';
        }
        LOG.flushStream();
        
        if (G.loadData(param.ovlFile.c_str()) == false) {
            LOG.oss << "[err] Cannot open \"" << param.ovlFile.c_str() << "\"" << endl;
            LOG.flushStream(TOSTDERR);
            exit(EXIT_FAILURE);
        }

        G.countEdges();
        Node::setEdgeSortedFlag(true);

        LOG.oss << "   reads length:             " << R->getReadsLength() << '\n';
        LOG.oss << "   number of reads:          " << R->getEffectiveNReads() << '\n';
        LOG.oss << "   number of nodes:          " << G.getNNodes() << '\n';
        LOG.oss << "   number of edges:          " << G.getNEdges() << '\n';
        LOG.oss << "   minimum overlap size:     " << G.getMinOverlap() << '\n';
        LOG.flushStream(TOSTDOUT);

        if (param.overlapCutoff != 0) {
            if (param.overlapCutoff > R->getReadsLength() - 1 ||
                    param.overlapCutoff < G.getMinOverlap()) {
                LOG.oss << "[err] improper overlaps cutoff value (-m)\n";
                LOG.oss << "[err] possible values are in the range ["
                        << G.getMinOverlap() << ";"
                        << R->getReadsLength() - 1 << "]" << endl;
                LOG.flushStream(TOSTDERR);
                exit(EXIT_FAILURE);
            }

            LOG.oss << "Overlaps size cutoff: " << param.overlapCutoff;
            LOG.flushStream(TOSTDOUT);
            cout << "Discarding overlaps shorter than " << param.overlapCutoff << "... " << flush;
            G.overlapSizeCutoff(param.overlapCutoff);
            cout << "done" << endl;
        }

        unsigned int nNodes = 0, nReads = 0;
        unsigned int sum = 0;

        if (param.cleanGraph) {

            if (param.checkGraph)
                G.checkConsistency();

            G.condense(true, false);
            LOG.oss << "   Graph has been condensed to "
                    << G.getNNodes()
                    << " nodes";
            LOG.flushStream(TOSTDOUT);


            if (param.discardNonUsable) {
                nNodes = G.discardShortOrphanNodes(nReads);
                LOG.oss << "Removed orphan non-usable nodes: "
                        << nNodes
                        << " ("
                        << (float) nReads / R->getEffectiveNReads()*100 << "% of the total number of reads)\n";
                LOG.flushStream(TOSTDOUT);

                //G.condense(true, true);
            }

            cout << "Removing dead-end paths... " << flush;
            sum = G.identifyDeadEnd(param.maxDeadEndLength, nReads);
            cout << "done\n";

            LOG.oss << "   Removed dead-ends: "
                    << sum
                    << " ("
                    << (float) nReads / R->getEffectiveNReads()*100
                    << "% of the total number of reads)" << endl;
            LOG.flushStream(TOSTDOUT);

            G.condense(true, true);
            LOG.oss << "   Graph has been condensed to "
                    << G.getNNodes()
                    << " nodes";
            LOG.flushStream(TOSTDOUT);

            if (param.checkGraph)
                G.checkConsistency();

            if (param.contextualCleaning) {
                //1) reliable cutoff
                //2) suspect cutoff (removed only if gas at least one reliable brother
                //3) G-incoherent cutoff (anyway removed)

                G.computeEdgesProb(CCreliableCutoff);
                sum = G.removeEdgesByValue(CCsupectCutoff, CCGincoherantCutoff);
                cout << "   " << sum << " edges have been cleaned out" << endl;

                LOG.oss << "Contextual cleaning cleared " << sum << " incoherent edges";
                LOG.flushStream();

                G.condense(true, true);
                LOG.oss << "   Graph has been condensed to "
                        << G.getNNodes()
                        << " nodes";
                LOG.flushStream(TOSTDOUT);


                cout << "Removing dead-end paths..." << flush;
                sum = G.identifyDeadEnd(param.maxDeadEndLength, nReads);
                cout << "done\n";

                LOG.oss << "   Removed dead-ends: "
                        << sum
                        << " ("
                        << (float) nReads / R->getEffectiveNReads()*100
                        << "% of the total number of reads)" << endl;
                LOG.flushStream(TOSTDOUT);

                G.condense(true, true);
                LOG.oss << "   Graph has been condensed to "
                        << G.getNNodes()
                        << " nodes";
                LOG.flushStream(TOSTDOUT);
            }

            if (param.discardNonUsable) {
                nNodes = G.discardShortOrphanNodes(nReads);
                LOG.oss << "Removed orphan non-usable nodes: " << nNodes << endl;
                LOG.flushStream(TOSTDOUT);
                G.condense(true, true);
                LOG.oss << "   Graph has been condensed to "
                        << G.getNNodes()
                        << " nodes";
                LOG.flushStream(TOSTDOUT);
            }

            G.estimateCoverage(param.minContigCoverage, param.targetSize);

            if (param.cleanBubbles) {
                nNodes = G.bubbles((double) param.minContigCoverage / 2);
                LOG.oss << "Resolved bubbles: " << nNodes << '\n';
                LOG.flushStream(TOSTDOUT);
                G.condense(true, true);
                LOG.oss << "   Graph has been condensed to "
                        << G.getNNodes()
                        << " nodes";
                LOG.flushStream(TOSTDOUT);
            }

            if (param.minNodeCoverage > 0.0) { //undocumented: NOT used by default: use with caution
                cout << "Node coverage cutoff..." << flush;
                nNodes = G.coverageCutoff(param.minNodeCoverage);
                cout << "done\n";
                LOG.oss << "Coverage cutoff: " << param.minNodeCoverage << " Nodes removed: " << nNodes << endl;
                G.condense(true, true);
                LOG.oss << "   Graph has been condensed to "
                        << G.getNNodes()
                        << " nodes";
                LOG.flushStream(TOSTDOUT);
            }

//            if (param.discardNonUsable) {
//                nNodes = G.discardShortOrphanNodes(nReads);
//                cout << "   " << nNodes << " nodes corresponding to " << nReads << " reads have been discarded ("
//                        << (float) nReads / R->getEffectiveNReads()*100 << "%)\n";
//                LOG.oss << "Removed orphan non-usable nodes: " << nNodes << endl;
//                LOG.flushStream(TOSTDOUT);
//            }

            if (param.discardSuspicousNodes) {
                cout << "Discarding suspicious nodes... " << flush;
                nNodes = G.discardSuspiciousNode(param.minContigCoverage);
                cout << "done\n";
                LOG.oss << "   Suspicious nodes removed: " << nNodes;
                LOG.flushStream(TOSTDOUT);
            }

            if (param.discardNonUsable || param.discardSuspicousNodes) {
                G.condense(true, true);
                LOG.oss << "   Graph has been condensed to "
                        << G.getNNodes()
                        << " nodes";
                LOG.flushStream(TOSTDOUT);
            }

            if (param.writeOSG) {
                //write overlaps string graph
                //smaller and can be used with the interactive mode
                cout << "Saving the overlaps string graph..." << flush;
                G.saveData((param.prefix + ".osg").c_str());
                cout << "done" << endl;
            }
        }

        if (param.cleanGraph == false)
            G.estimateCoverage(param.minContigCoverage, param.targetSize);

        G.estimatePairedDistance(param.shortPeHorizon, param.longPeHorizon, nsd, param.prefix);

        if (param.checkGraph)
            G.checkConsistency();

        for (unsigned int i = 1; i <= G.getNNodes(); i++)
            G.nodesTab[i].unsetVisited();

        if (!param.interactiveShell) {//OUTPUT CONTIGS
            if (param.minContigSize == 0)
                param.minContigSize = R->getReadsLength() + R->getReadsLength() / 2;

            if (param.trimRed == false)
                param.maxContigRedundancy = P->getGlobalMaxAllowedDistance();

            G.assemble(param.prefix,
                    param.minContigSize,
                    param.minContigCoverage,
                    param.trimContigEnds,
                    param.minNPair,
                    param.minRatio,
                    param.maxContigRedundancy);

        } else {//DEV SHELL

            //            for (unsigned int i = 1; i < G.getNNodes(); i++)
            //            {
            //                if (G.nodesTab[i].getNReads() > 10 &&
            //                        L->isUniDirectional(G.nodesTab[i].getLayout()))
            //                {
            //                    if (!done)
            //                    {
            //                        cout << "warning: unidirectional nodes have been detected\n";
            //                        done=true;
            //                    }
            //                    cout << "node " << i << " (" << L->getNReads(G.nodesTab[i].getLayout()) << " reads)" << endl;
            //                    // G.nodesTab[i].isolate();
            //                }
            //            }

            G.computeEdgesProb(CCreliableCutoff);
            DevShell DEV;
            DEV.init(&G, P, R, L);
            DEV.prompt();
        }
    }

    delete P;
    delete L;
    delete R;
    LOG.close();
    return 0;
}

