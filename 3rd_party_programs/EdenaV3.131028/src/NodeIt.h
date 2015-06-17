/* 
 * File:   NodeIt.h
 * Author: david
 *
 * Created on March 16, 2011, 11:40 AM
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

#ifndef NODEIT_H
#define	NODEIT_H
#include <cstdlib>
#include <string>
#include "node.h"
using namespace std;
class Node;

struct NodeIt
{
    friend class Node;

public:
    NodeIt();
    NodeIt(const NodeIt& orig);
    virtual ~NodeIt();

    void clear();
    void initNodeIt(unsigned int node, bool dir);
    void initIterator(); //set iterator to first son
    void setIteratorToEnd();
    unsigned int getNNeighbor();
    unsigned int getNBackNeighbor();

    inline bool getNext(NodeIt& myNodeIt)
    {
        if (edgeCounter == nEdge)
        {
            myNodeIt.clear();
            return false;
        }
        pNode->getNeighbor(dir, edgeCounter, myNodeIt);
        edgeCounter++;
        return true;
    }
    
    NodeIt& getNeighbor(size_t index, NodeIt&) const;
    
    //check weather index is consistent with current orientation
    NodeIt& getNeighborGivenAbsoluteIndex(size_t index, NodeIt&) const;
    
    void reverse();
    inline bool isTraversed() {return pNode->isVisited(dir);}
    inline bool isTraversedAnyDir() {return pNode->isVisited();}
    inline void setTraversed(){pNode->setVisited(dir);}
    inline bool isInplay() {return pNode->isVisited(dir);}
    inline void setInplay(){pNode->setVisited(dir);}
    
    inline bool isReduced(){return pNode->isReduced();}
    inline void setReduced(){pNode->setReduced();}
    inline void setDiscarded(){pNode->setDiscarded();}
    inline bool isDiscarded(){return pNode->isDiscarded();}
    inline void unsetDiscarded() {pNode->unsetDiscarded();}
    inline unsigned short int getValue(){return pNode->value;}
    inline void setValue(unsigned short int v){pNode->value=v;}

    inline unsigned int getArrivalEdgeSize() const {return arrivalEdgeSize;}
    inline unsigned int getAbsArrivalEdgeIndex() const {return absArrivalEdgeIndex;}
    unsigned int getNodeLength() const;
    void getSequence(string&);
    inline double getCoverage(){return pNode->getCoverage();}
    unsigned int getNodeId() const;
    unsigned int getLayout() const; //(layout = last read id in the node)
    inline unsigned int getLastReadInLayout() const {return getLayout();}
    unsigned int getFirstReadInLayout() const;
    bool getDirection() const;
    inline bool isNull() const {return pNode==0x0;}
    inline void discard(){pNode->isolate();pNode->setDiscarded();}
    inline void setMultipleEdges() {pNode->setMultipleEdges(dir);}
    inline bool hasMultipleEdges() {return pNode->hasMultipleEdges(dir);}
    inline void unsetMultiplesEdges(){pNode->unsetMultiplesEdges();}
    inline bool isEliminated() {return pNode->isEliminated(dir);}
    inline void setEliminated() {pNode->setEliminated(dir);}
    inline void unsetEliminated() {pNode->unsetEliminated(dir);}
    static Node *N;
    
    bool operator==(const NodeIt &a) const {
        return (dir==a.dir && pNode==a.pNode);
        }
    bool operator!=(const NodeIt &a) const {
        return !(dir!=a.dir || pNode!=a.pNode);
    }

private:

    Node *pNode;
    bool dir;
    unsigned int nEdge;
    unsigned int arrivalEdgeSize;
    unsigned int absArrivalEdgeIndex; //required for path
    size_t edgeCounter;

};

#endif	/* NODEIT_H */

