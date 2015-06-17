/* 
 * File:   PEMatches.cpp
 * Author: david
 * 
 * Created on April 8, 2012, 11:34 PM
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

#include "PEMatches.h"
#include "globalFunc.h"


#include <cstdlib>
#include <limits>
using namespace std;


Pairing* PEMatches::P = 0x0;
unsigned int PEMatches::nLibrary = 0;

PEMatchesStorage::PEMatchesStorage(Pairing *p)
{
    P=p;
    PEMatches::P=p;
    PEMatches::nLibrary=p->getNLibrary();
    nElement=0;  
    maxSize=10000;
    
    allocate();
}

PEMatchesStorage::~PEMatchesStorage()
{
    freeMemory();
}

void PEMatchesStorage::allocate()
{
    V = (PEMatches**)malloc(maxSize*sizeof(PEMatches*));
    for (size_t i=0; i<maxSize; i++)
        V[i]=0x0;
}

void PEMatchesStorage::doubleAllocatedSize()
{
    void *p=(PEMatches**)realloc((void*)V,(maxSize*2)*sizeof(PEMatches*));
    
    if (p==0x0)
    {
        cout << "PEMatchesStorage::doubleAllocatedSize()\n";
        cout << "Out of memory\n";
        exit(0);        
    }
    
    V=(PEMatches**)p;
    
    for (size_t i=maxSize; i<(maxSize*2); i++)
        V[i]=0x0;
    
    maxSize*=2;
}

void PEMatchesStorage::freeMemory()
{
    for (unsigned int i=0; i<maxSize; i++)
    {
        if (V[i]!=0x0)
            delete V[i]; //allocated by "new"    
    }
    free(V);
}

void PEMatchesStorage::initStorage()
{
    nElement=0;
}

PEMatches* PEMatchesStorage::getPPEMatches(unsigned int nChoice)
{
    nElement++;
    
    if (nElement >= maxSize)
    {
        doubleAllocatedSize();
    }
   
     if (V[nElement-1]==0x0)
         V[nElement-1]=new PEMatches;
         
     
    V[nElement-1]->init(nChoice);
   
    
    return V[nElement-1];
}


PEMatches::PEMatches()
{
    counts=0x0;
    sumDist=0x0;
    sumDistSquare=0x0;
    maxBranch=10;
    nBranch=0;

    allocate();
}

PEMatches::PEMatches(const PEMatches& orig)
{
}

PEMatches::~PEMatches()
{
    freeMemory();
}

void PEMatches::allocate()
{
    if (nLibrary==0)
    {
        cout << "PEMatches::allocate(...)\n";
        cout << "This should mot happen\n";
        sendBugReportPlease(cout);
        exit(0);
    }
    
    counts = (unsigned int**)malloc(nLibrary*sizeof(unsigned int*));
    sumDist = (double**)malloc(nLibrary*sizeof(double*));
    sumDistSquare = (double**)malloc(nLibrary*sizeof(double*));
    
    for (unsigned int i=0; i<nLibrary; i++)
    {
        counts[i] = (unsigned int*)malloc(maxBranch*sizeof(unsigned int));
        sumDist[i] = (double*)malloc(maxBranch*sizeof(double));
        sumDistSquare[i] = (double*)malloc(maxBranch*sizeof(double));
    } 
}

void PEMatches::doubleAllocationSize()
{
    
    void *p;
    for (unsigned int i=0; i<nLibrary; i++)
    {
        p=realloc( (void*)counts[i],(maxBranch*2)*sizeof(unsigned int));
        if (p==0x0)
            break;
        counts[i]=(unsigned int*)p;
        p=realloc( (void*)sumDist[i],(maxBranch*2)*sizeof(double));
        if (p==0x0)
            break;
        sumDist[i]=(double*)p;
        
        p=realloc( (void*)sumDistSquare[i],(maxBranch*2)*sizeof(double));
        if (p==0x0)
            break;
        sumDistSquare[i]=(double*)p;
    }
    
    if (p==0x0)
    {
        cout << "PEMatches::doubleAllocationSize()\n";
        cout << "Out of memory\n";
        exit(0);
    }
    maxBranch*=2;
}

void PEMatches::freeMemory()
{
     for (unsigned int i=0; i<nLibrary; i++)
     {
         free(counts[i]);
         free(sumDist[i]);
         free(sumDistSquare[i]);  
     }
     
     free(counts);
     free(sumDist);
     free(sumDistSquare);
     
     counts=0x0;
     sumDist=0x0;
     sumDistSquare=0x0;
}

void PEMatches::init(unsigned int n)
{
    if (n > maxBranch)
    {
        while (n>maxBranch)
        {
            doubleAllocationSize();
        }
    }

    nBranch = n;

    for (size_t i = 0; i < nLibrary; i++)
    {
        for (size_t j = 0; j < nBranch; j++)
        {
            counts[i][j] = 0;
            sumDist[i][j] = 0.0;
            sumDistSquare[i][j] = 0.0;
        }
    }
}

double PEMatches::getVar(unsigned int lib, unsigned int i)
{
    lib--;
    
    if (counts[lib][i]<2)
        return 0.0;
    
    double mean=sumDist[lib][i]/counts[lib][i];

    return (1.0 / (counts[lib][i] - 1))*(sumDistSquare[lib][i] - counts[lib][i] * mean * mean);
}

void PEMatches::getMaxCount(unsigned lib, unsigned int &sumCounts, unsigned int &max, size_t&maxIndex)
{
    sumCounts=0;
    max=counts[lib-1][0];
    sumCounts+=max;
    maxIndex=0;
    
    for (unsigned int i=1; i<nBranch; i++)
    {
        if (counts[lib-1][i] > max)
        {
            max=counts[lib-1][i];
            maxIndex=i;
        }
        sumCounts+=counts[lib-1][i];
    }
}