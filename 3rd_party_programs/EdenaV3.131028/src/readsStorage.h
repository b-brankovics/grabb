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

#ifndef STREAD_H
#define STREAD_H
#include <cstdlib>
#include <string>
#include <cstring>
#include <cstring>
#include <iostream>
#include <fstream>
#include <set>
#include "Pairing.h"
#include "readsLayout.h"
#include "customString.h"
using namespace std;

struct _OV {
    unsigned int id;
    short size;
};

struct orderSet
{

    bool operator ()(const _OV& a, const _OV & b) const
    {
        return a.id > b.id;
    }
};

struct mateInsert
{
    unsigned int r1;
    char dir1;
    unsigned int r2;
    char dir2;
    
    inline bool operator< (const mateInsert& b) const
    {      
        if (dir1 != b.dir1)
            return dir1 > b.dir1;
        else if (dir2 != b.dir2)
            return dir2 > b.dir2;
        else if (r1 != b.r1)
            return r1 > b.r1;
        else
            return r2 > b.r2;
    }
};

class ReadsStorage
{
    friend class OverlapsGraph;
    friend class Node;

    //Store the reads sequences
    //redundant reads are stored once
    //Assume all reads with the same length

public:
    ReadsStorage();
    virtual ~ReadsStorage();

protected:
private:

    unsigned int nNR_reads; //in the non-redundant set
    unsigned int nReadsUpperBound; //total number
    unsigned int effectiveNReads;
    char* pp;
    static char *READS; //non redundant

    unsigned int *nrReadIds; //position in the NR reads set
    char *nrReadDirection;
    char *flags;

    unsigned int nDiscardedReads;
    int readsLength;

    char** prefixTableD;
    char** prefixTableR;

    int minOvSize;
    static char REVERSE[256];
    static char TOUPPER[256];

    inline void compReverse(const char* a, char *b, int l)
    {
        b[l]='\0';
        l--;
        
        while(l>=0)
        {
            b[l]=REVERSE[(size_t)*a];
            a++;
            l--;
        }
    }
    
    inline bool toUpperAndCheckDNA(char *a)
    {
        while (*a!='\0')
        {
            *a=TOUPPER[(size_t)*a];
            if (*a == 'X')
                return false;
            a++;
        }
        return true;
    }

   
    inline void r_strcpy(char* a, const char*b)
    {
        while (*b!='\0')
        {
            *a=REVERSE[(size_t)*b];
            a++;
            b--;
        }
        *a='\0';
    }

   
   
    inline unsigned int getId(const char*p) const
    {
        return 1+(p-READS-1) / (readsLength+1);
    }
    inline unsigned int getId(size_t p) const
    {
        return 1+ (p-1) / (readsLength+1);
    }
    
    struct Read_pointer
    {
        char* p;
    };

    struct d_orderSet
    {

        bool operator ()(const char* a, const char* b) const
        {

            if (strcmp(a, b) < 0)
                return true;
            else
                return false;

        }
    };

    struct r_orderSet
    {

        bool operator ()(char* a, char* b) const
        {

            if (r_strcmp(a, b) < 0)
                return true;
            else
                return false;

        }
    };

    //ordered non-redundant reads
    set<char*, d_orderSet> d_set;

    //ordered reverse but NOT complement non.redundant reads
    //Read_pointer.p refer to the LAST char of the read 
    set<char*, r_orderSet> r_set;

   
    set<mateInsert> mateInsertBTree;



public:

    void init(unsigned int, unsigned int);
    void allocate();
    void adjustAllocation();
    void initPrefixTables();
    void cleanPrefixTables();
    void freeMemory();
    
    bool load(istream&);
    void save(ostream&);

    unsigned int insertRead(const char* p);
    inline void initDuplicateCheck(){mateInsertBTree.clear();}
    bool isDuplicateMate(unsigned int r1, unsigned int r2);
   
    //use STL::set
    unsigned int determineOverlaps(unsigned int nrId,
            vector<unsigned int>& ID,
            vector<short>& OV,
            multiset<_OV, orderSet>& ovSet);
    
    //use prefix tables
    unsigned int determineOverlaps2(unsigned int nrId,
            vector<unsigned int>& ID,
            vector<short>& OV,
            multiset<_OV, orderSet>& ovSet);

    int loadReadsFiles(
            string infile1,
            string infile2,
            int mateOrientation,
            Pairing*,
            ReadsLayout*);
    
    char** d_prefix_lowerBound(const char*);
    char** r_prefix_lowerBound(const char*);

    inline void getDirectNR(char* dest, size_t nrId)
    {
        strcpy(dest,(1+READS) + (nrId - 1)*(readsLength + 1));
    }

    inline void getReverseNR(char* dest, size_t nrId)
    {
        r_strcpy(dest,(1+READS) + (nrId - 1)*(readsLength + 1) +readsLength-1);
    }
    
    inline void getDirectSequence(char* dest, size_t id) {
        if (nrReadDirection[id])
            getDirectNR(dest, nrReadIds[id]);
        else
            getReverseNR(dest, nrReadIds[id]);
    }

    inline void getReverseSequence(char* dest, size_t id) {
        if (!(nrReadDirection[id]))
            getDirectNR(dest, nrReadIds[id]);
        else
            getReverseNR(dest, nrReadIds[id]);
    }
    
    inline int getReadsLength() const
    {
        return readsLength;
    }

    inline unsigned int getN_nrReads() const
    {
        return nNR_reads;
    };

    inline unsigned int getEffectiveNReads() const
    {
        return effectiveNReads;
    };

    inline unsigned int getNDiscaredReads() const
    {
        return nDiscardedReads;
    }
    inline void setMinOvSize(unsigned int m){minOvSize=m;}
    
    inline void setFlag(unsigned int rId, char state){flags[rId]=state;}
    inline char getFlag(unsigned int rId){return flags[rId];}
};


#endif // STREAD_H
