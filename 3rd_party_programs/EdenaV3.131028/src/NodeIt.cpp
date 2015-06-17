/* 
 * File:   NodeIt.cpp
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

#include "NodeIt.h"
//#include "node.h"

Node* NodeIt::N=0x0;

NodeIt::NodeIt()
{
    pNode=0x0;
}

NodeIt::NodeIt(const NodeIt& orig)
{
    pNode = orig.pNode;
    dir = orig.dir;
    nEdge = orig.nEdge;
    arrivalEdgeSize = orig.arrivalEdgeSize;
    absArrivalEdgeIndex = orig.absArrivalEdgeIndex; //required for path
    edgeCounter = orig.edgeCounter;
}

NodeIt::~NodeIt()
{
}

void NodeIt::clear()
{
    pNode=0x0;
    nEdge=0;
}

void NodeIt::initNodeIt(unsigned int node, bool direction)
{
    pNode=N+node;
    dir=direction;
    initIterator();
    absArrivalEdgeIndex=node; //initialized with current node ID
                               //required for paths construction
}

void NodeIt::initIterator()
{
    nEdge=pNode->getNEdges(dir);
    edgeCounter=0;
}

void NodeIt::setIteratorToEnd()
{
    edgeCounter = nEdge;
}

unsigned int NodeIt::getNNeighbor()
{
    return nEdge;
}

unsigned int NodeIt::getNBackNeighbor()
{
    return pNode->getNEdges(!dir);
}

//bool NodeIt::getNext(NodeIt &myNodeIt)
//{
//    if (edgeCounter == nEdge)
//    {
//        myNodeIt.clear();
//        //return myNodeIt;
//        return false;
//    }
//
//    pNode->getNeighbor(dir,edgeCounter,myNodeIt);
//    edgeCounter++;
//
//    return true;
//}

NodeIt& NodeIt::getNeighbor(size_t index, NodeIt &myNodeIt) const
{
    if (index >= nEdge)
    {
        myNodeIt.clear();
        return myNodeIt;
    }

    pNode->getNeighbor(dir,index,myNodeIt);
    return myNodeIt;
}

NodeIt& NodeIt::getNeighborGivenAbsoluteIndex(size_t absIndex, NodeIt &myNodeIt) const
{
    size_t index;
    
    if (!dir)
    {
        index=absIndex-pNode->nOvRight;
    }
    else
        index=absIndex;
        
  
    if (index >= nEdge)
    {
        myNodeIt.clear();
        return myNodeIt;
    }

    pNode->getNeighbor(dir,index,myNodeIt);
    return myNodeIt;
}
void NodeIt::reverse()
{
    dir=!dir;
    initIterator();
}

unsigned int NodeIt::getNodeLength() const
{
    return pNode->getSequenceLength();
}

void NodeIt::getSequence(string& s)
{
    s=pNode->getSequence(dir);
}

unsigned int NodeIt::getNodeId() const
{
    return pNode->getThisId();
}

unsigned int NodeIt::getLayout() const
{
    return pNode->getLayout();
}

unsigned int NodeIt::getFirstReadInLayout() const
{
    return  pNode->L->getBegin(pNode->getLayout());
}

bool NodeIt::getDirection() const
{
    return dir;
}
