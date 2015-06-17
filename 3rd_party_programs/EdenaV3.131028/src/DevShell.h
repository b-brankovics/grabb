/* 
 * File:   DevShell.h
 * Author: david
 *
 * Created on November 26, 2012, 2:16 PM
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

#ifndef DEVSHELL_H
#define	DEVSHELL_H

#include "overlapGraph.h"
#include "Pairing.h"
#include "readsLayout.h"
#include "readsStorage.h"
#include "BackwardsWalker.h"

#include<vector>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

class DevShell
{
public:
    DevShell();
    DevShell(const DevShell& orig);
    virtual ~DevShell();
    
    
    void init(OverlapsGraph *g, Pairing *p, ReadsStorage *r, ReadsLayout *l);
    void prompt();
    
private:
    
    vector<unsigned int> ePath;
    size_t p1; //two pointer on ePath that indicates the boundaries for either
    size_t p2; //desired amplicon or PE matches counts
    vector<unsigned int> pathNode;
    vector<bool> pathDir;
    string nodeList;
    
    unsigned int nodeID;
    bool nodeDir;
    vector<string> args;
    istringstream iss;
    bool checkDistances;
    int minNPair;
    double minRatio;
    
    BackwardsWalker BKW;
    unsigned int chosen;
    
    OverlapsGraph *G;
    Pairing *P;
    ReadsStorage *R;
    ReadsLayout *L;
    
    void parsePath(string rawPath);
    void parseNodeId(string buffer);

    inline void whatDoYouMean(ostream &out)
    {
        out << "??" << endl;
    }

    void printHelp(ostream &out);
    
    void growSearchTree();
    void elongateSearchTree();
    void Tree2seq();
    void drawGraph();
    void path2Seq();
    void reversePath();
    void printReadsLayout();
    void locateRead();
    void setCheckDistances();
    void p3(); //get input for primer3 website
    bool set_p1p2_ePath(); //from args(1)
    void dev();
    

    

};

#endif	/* DEVSHELL_H */

