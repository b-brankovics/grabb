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

#include "readsStorage.h"
#include "globalFunc.h"
#include "logWriter.h"
#include "crc.h"
#include <cstdlib>
using namespace std;

extern logWriter LOG;

char *ReadsStorage::READS = 0x0;
char ReadsStorage::REVERSE[256];
char ReadsStorage::TOUPPER[256];

ReadsStorage::ReadsStorage() {
    READS = 0x0;
    prefixTableD = 0x0;
    prefixTableR = 0x0;
    nrReadIds = 0x0;
    nrReadDirection = 0x0;
    flags = 0x0;
    readsLength = 0;

    memset(REVERSE, (int) 'X', 256 * sizeof (char));
    memset(TOUPPER, (int) 'X', 256 * sizeof (char));

    REVERSE[(int) 'A'] = 'T';
    REVERSE[(int) 'T'] = 'A';
    REVERSE[(int) 'C'] = 'G';
    REVERSE[(int) 'G'] = 'C';
    TOUPPER[(int) 'a'] = 'A';
    TOUPPER[(int) 't'] = 'T';
    TOUPPER[(int) 'c'] = 'C';
    TOUPPER[(int) 'g'] = 'G';
    TOUPPER[(int) 'A'] = 'A';
    TOUPPER[(int) 'T'] = 'T';
    TOUPPER[(int) 'C'] = 'C';
    TOUPPER[(int) 'G'] = 'G';
}

ReadsStorage::~ReadsStorage() {
    freeMemory();
}

void ReadsStorage::freeMemory() {
    cleanPrefixTables();
    if (READS != 0x0) {
        free(READS);
        READS = 0x0;
    }
    if (nrReadIds != 0x0) {
        free(nrReadIds);
        nrReadIds = 0x0;
    }
    if (nrReadDirection != 0x0) {
        free(nrReadDirection);
        nrReadDirection = 0x0;
    }
    if (flags != 0x0) {
        free(flags);
        flags = 0x0;
    }
}

//Mandatory to call

void ReadsStorage::init(unsigned int upperBoundSize, unsigned int rl) {
    if (upperBoundSize == 0 || rl == 0)
        return;
    readsLength = rl;
    nReadsUpperBound = upperBoundSize;
    nNR_reads = upperBoundSize;
    effectiveNReads = upperBoundSize;
    nDiscardedReads = 0;

    allocate();
    pp = READS + 1;

    nNR_reads = 0;
    effectiveNReads = 0;
}

void ReadsStorage::allocate() {
    freeMemory();
    READS = (char*) malloc(1 + (size_t) (nNR_reads)*(readsLength + 1) * sizeof (char));
    READS[0] = '\0';
    nrReadIds = (unsigned int*) malloc((effectiveNReads + 1) * sizeof (unsigned int));
    nrReadDirection = (char*) malloc((effectiveNReads + 1) * sizeof (char));
    flags = (char*) malloc((effectiveNReads + 1) * sizeof (char));
    memset(flags, 0, (effectiveNReads + 1) * sizeof (char));
}

void ReadsStorage::adjustAllocation() {
    void* p;
    p = realloc((void*) READS, 1 + (size_t) (nNR_reads) * (readsLength + 1));
    READS = (char*) p;
    p = realloc((void*) nrReadIds, (effectiveNReads + 1) * sizeof (unsigned int));
    nrReadIds = (unsigned int*) p;
    p = realloc((void*) nrReadDirection, (effectiveNReads + 1) * sizeof (char));
    nrReadDirection = (char*) p;
    p = realloc((void*) flags, (effectiveNReads + 1) * sizeof (char));
    flags = (char*) p;
    memset(flags, 0, (effectiveNReads + 1) * sizeof (char));
}

void ReadsStorage::initPrefixTables() {
    cleanPrefixTables();

    try {
        //prefix tables are 1-indexed
        prefixTableD = new char* [nNR_reads + 1];
        prefixTableR = new char* [nNR_reads + 1];
    }    catch (bad_alloc ex) {
        LOG.oss << ex.what() << "\n[err] Not enough memory to allocate suffix table" << '\n';
        LOG.flushStream(TOSTDERR);
        exit(EXIT_FAILURE);
    }

    set<char*, d_orderSet>::iterator dsIt;
    set<char*, r_orderSet>::iterator rsIt;
    char** p;

    if (d_set.size() != r_set.size()) {
        cerr << "[bug] ReadsStorage::initPrefixTables() problem\n";
        sendBugReportPlease(cerr);
    }

    prefixTableD[0] = 0x0;
    prefixTableR[0] = 0x0;

    p = prefixTableD + 1;
    for (dsIt = d_set.begin(); dsIt != d_set.end(); dsIt++) {
        *p = *dsIt;
        p++;
    }

    p = prefixTableR + 1;
    for (rsIt = r_set.begin(); rsIt != r_set.end(); rsIt++) {
        *p = *rsIt;
        p++;
    }

    d_set.clear();
    r_set.clear();

}

void ReadsStorage::cleanPrefixTables() {
    if (prefixTableD != 0x0) {
        delete [] prefixTableD;
        prefixTableD = 0x0;
    }
    if (prefixTableR != 0x0) {
        delete [] prefixTableR;
        prefixTableR = 0x0;
    }
}

bool ReadsStorage::load(istream &in_bin) {

    unsigned int crcCheck;
    in_bin.read((char*) &crcCheck, sizeof (unsigned int));


    in_bin.read((char*) &nNR_reads, sizeof (unsigned int));
    in_bin.read((char*) &effectiveNReads, sizeof (unsigned int));
    in_bin.read((char*) &readsLength, sizeof (unsigned int));


    if (readsLength > 10000) {
        LOG.oss << "[err] The file is not compatible with current version.\n";
        LOG.oss << "[err] Please, generate a new ovl file\n";
        LOG.flushStream(TOSTDERR);
        exit(EXIT_FAILURE);
    }
    allocate();

    in_bin.read((char*) READS, sizeof (char) *(size_t) nNR_reads * (readsLength + 1));
    in_bin.read((char*) nrReadIds, sizeof (unsigned int)* (effectiveNReads + 1));
    in_bin.read((char*) nrReadDirection, sizeof (char)* (effectiveNReads + 1));

    Crc32 CRC;

    CRC.AddData((uint8_t*) & nNR_reads, sizeof (unsigned int));
    CRC.AddData((uint8_t*) & effectiveNReads, sizeof (unsigned int));
    CRC.AddData((uint8_t*) & readsLength, sizeof (unsigned int));
    CRC.AddData((uint8_t*) READS, sizeof (char)*(size_t) nNR_reads * (readsLength + 1));
    CRC.AddData((uint8_t*) nrReadIds, sizeof (unsigned int)* (effectiveNReads + 1));
    CRC.AddData((uint8_t*) nrReadDirection, sizeof (char)* (effectiveNReads + 1));

    unsigned int crc = CRC.GetCrc32();
    if (crcCheck != crc)
        return false;

    return true;
}

void ReadsStorage::save(ostream &out_bin) {
    Crc32 CRC;

    CRC.AddData((uint8_t*) & nNR_reads, sizeof (unsigned int));
    CRC.AddData((uint8_t*) & effectiveNReads, sizeof (unsigned int));
    CRC.AddData((uint8_t*) & readsLength, sizeof (unsigned int));
    CRC.AddData((uint8_t*) READS, sizeof (char) *(size_t) nNR_reads * (readsLength + 1));
    CRC.AddData((uint8_t*) nrReadIds, sizeof (unsigned int)* (effectiveNReads + 1));
    CRC.AddData((uint8_t*) nrReadDirection, sizeof (char)* (effectiveNReads + 1));

    unsigned int crcCheck = CRC.GetCrc32();
    out_bin.write((char*) &crcCheck, sizeof (unsigned int));

    out_bin.write((char*) &nNR_reads, sizeof (unsigned int));
    out_bin.write((char*) &effectiveNReads, sizeof (unsigned int));
    out_bin.write((char*) &readsLength, sizeof (unsigned int));
    out_bin.write((char*) READS, sizeof (char)*(size_t) nNR_reads * (readsLength + 1));
    out_bin.write((char*) nrReadIds, sizeof (unsigned int)* (effectiveNReads + 1));
    out_bin.write((char*) nrReadDirection, sizeof (char)* (effectiveNReads + 1));
}

unsigned int ReadsStorage::insertRead(const char* p) {
    set<char*, d_orderSet>::iterator dsIt;
    pair < set<char*, d_orderSet>::iterator, bool> ret;
    bool dir = true;
    unsigned int id = 0;
  

    strcpy(pp, p);
    compReverse(pp, pp + readsLength + 1, readsLength);
    *(pp + readsLength + 1 + readsLength) = '\0';

    //search for the reverse comp seq
    dsIt = d_set.find(pp + readsLength + 1);

    if (dsIt == d_set.end()) {
        dir = true;
        ret = d_set.insert(pp);
        r_set.insert(pp + readsLength - 1); //last nucl ("reverse" c-string)

        if (ret.second == true)//new element inserted
        {
            pp += readsLength + 1;
            nNR_reads++;
            id = nNR_reads;
        } else {
            id = getId(*(ret.first));
        }
    } else {//already exist reverse complemented
        id = getId(*dsIt);
        dir = false;
    }

    if (nNR_reads != (pp - READS) / (readsLength + 1)) {
        cerr << "[bug] ReadsStorage::insertRead(const char* p)" << endl;
        sendBugReportPlease(cerr);
    }

    effectiveNReads++;
    nrReadIds[effectiveNReads] = id;
    nrReadDirection[effectiveNReads] = dir;

    return id;
}

bool ReadsStorage::isDuplicateMate(unsigned int r1, unsigned int r2) {
    static mateInsert M;
    pair < set<mateInsert>::iterator, bool> ret;

    M.dir1 = nrReadDirection[r1];
    M.dir2 = nrReadDirection[r2];
    M.r1 = nrReadIds[r1];
    M.r2 = nrReadIds[r2];

    ret = mateInsertBTree.insert(M);

    //ret.second is set to true if a new element has been insert
    return !ret.second;
}

unsigned int ReadsStorage::determineOverlaps(unsigned int nrId,
        vector<unsigned int>& ID,
        vector<short>& OV,
        multiset<_OV, orderSet>& ovSet) {
    //BINARY TREE VERSION

    //determine all overlaps for reads nrId and store
    //them in ID and OV. ovSet is provided as a parameter
    //to be thread safe

    set<char*, d_orderSet>::iterator dsIt;
    set<char*, r_orderSet>::iterator rsIt;

    char d_suffix[1000];
    char rc_suffix[1000];

    getDirectNR(d_suffix, nrId);
    getReverseNR(rc_suffix, nrId);

    unsigned int nOv = 0;
    unsigned int nOvRight = 0;
    _OV myOv;
    //    Read_pointer d_index; //for searching in d_set
    //    Read_pointer r_index; //for searching in r_set
    char* d_index; //for searching in d_set
    char* r_index; //for searching in r_set

    nOv = 0;
    nOvRight = 0;

    ID.clear();
    OV.clear();

    d_index = d_suffix;
    r_index = rc_suffix + readsLength - 1;

    if (minOvSize > readsLength)
        minOvSize = readsLength;

    for (short j = 1; j <= readsLength - minOvSize; j++) {
        //   >>>>>----
        //        >>>>>>>>>
        ovSet.clear();

        d_index++;
        r_index--;

        dsIt = d_set.lower_bound(d_index);

        while (dsIt != d_set.end()) {
            if (d_prefixOf(d_index, *dsIt)) {
                myOv.id = getId(*dsIt);
                myOv.size = (readsLength - j);
                ovSet.insert(myOv);
                dsIt++;
                nOv++;
                nOvRight++;
            } else
                break;
        }

        //   >>>>>----
        //        <<<<<<<<<

        rsIt = r_set.lower_bound(r_index);

        while (rsIt != r_set.end()) {
            if (r_prefixOf(r_index, *rsIt)) {
                myOv.id = getId(*rsIt);
                myOv.size = -(readsLength - j);
                ovSet.insert(myOv);
                rsIt++;
                nOv++;
                nOvRight++;
            } else
                break;
        }

        for (set<_OV>::iterator sIt = ovSet.begin(); sIt != ovSet.end(); sIt++) {
            ID.push_back(sIt->id);
            OV.push_back(sIt->size);
        }
    }

    d_index = rc_suffix;
    r_index = d_suffix + readsLength - 1;

    for (short j = 1; j <= readsLength - minOvSize; j++) {
        ovSet.clear();
        d_index++;
        r_index--;

        //      ---->>>>>
        //  <<<<<<<<

        rsIt = d_set.lower_bound(d_index);

        while (rsIt != d_set.end()) {
            if (d_prefixOf(d_index, *rsIt)) {
                myOv.id = getId(*rsIt);
                myOv.size = -(readsLength - j);
                ovSet.insert(myOv);
                rsIt++;
                nOv++;
            } else
                break;
        }

        //      ---->>>>>
        //    >>>>>>>

        dsIt = r_set.lower_bound(r_index);

        while (dsIt != r_set.end()) {
            if (r_prefixOf(r_index, *dsIt)) {
                myOv.id = getId(*dsIt);
                myOv.size = (readsLength - j);
                ovSet.insert(myOv);
                dsIt++;
                nOv++;
            } else
                break;
        }

        for (set<_OV>::iterator sIt = ovSet.begin(); sIt != ovSet.end(); sIt++) {
            ID.push_back(sIt->id);
            OV.push_back(sIt->size);
        }
    }

    return nOvRight;
}

unsigned int ReadsStorage::determineOverlaps2(unsigned int nrId,
        vector<unsigned int>& ID,
        vector<short>& OV,
        multiset<_OV, orderSet>& ovSet) {
    
    //PREFIX TABLE VERSION (slightly more efficient)

    //determine all overlaps for reads nrId and store
    //them in ID and OV. ovSet is provided as a parameter
    //to be thread safe
    //ovSet keeps IDs sorted
    
    //final overlaps (edges) vectors ID and OV are sorted:
    // firstly on the overlap size
    // secondly on the overlapped reads ID.
    
    char d_suffix[1000];
    char rc_suffix[1000];

    // "reverse" c-string
    d_suffix[0] = '\0';
    rc_suffix[0] = '\0';
    getDirectNR(d_suffix + 1, nrId);
    getReverseNR(rc_suffix + 1, nrId);

    unsigned int nOv = 0;
    unsigned int nOvRight = 0;
    _OV myOv;

    char *d_index; //for searching in d_set
    char *r_index; //for searching in r_set
    char** pTableIndex;
    nOv = 0;
    nOvRight = 0;

    ID.clear();
    OV.clear();

    d_index = d_suffix + 1;
    r_index = rc_suffix + readsLength;

    if (minOvSize > readsLength)
        minOvSize = readsLength;

    for (short j = 1; j <= readsLength - minOvSize; j++) {
        //   >>>>>----
        //        I>>>>>>>>>I

        ovSet.clear();

        d_index++;
        r_index--;

        pTableIndex = d_prefix_lowerBound(d_index);

        while (pTableIndex <= prefixTableD + nNR_reads) {
            if (d_prefixOf(d_index, *pTableIndex)) {
                myOv.id = getId(*pTableIndex);
                myOv.size = (readsLength - j);
                ovSet.insert(myOv);
                pTableIndex++;
                nOv++;
                nOvRight++;
            } else
                break;
        }

        //   >>>>>----
        //        I<<<<<<<<<I

        pTableIndex = r_prefix_lowerBound(r_index);

        while (pTableIndex <= prefixTableR + nNR_reads) {
            if (r_prefixOf(r_index, *pTableIndex)) {
                myOv.id = getId(*pTableIndex);
                myOv.size = -(readsLength - j);
                ovSet.insert(myOv);
                pTableIndex++;
                nOv++;
                nOvRight++;
            } else
                break;
        }

        for (set<_OV>::iterator sIt = ovSet.begin(); sIt != ovSet.end(); sIt++) {
            ID.push_back(sIt->id);
            OV.push_back(sIt->size);
        }
    }

    d_index = rc_suffix + 1;
    r_index = d_suffix + readsLength;

    for (short j = 1; j <= readsLength - minOvSize; j++) {
        ovSet.clear();
        d_index++;
        r_index--;

        //   <<<<<<-----
        //        I>>>>>>>>>>I

        pTableIndex = d_prefix_lowerBound(d_index);

        while (pTableIndex <= prefixTableD + nNR_reads) {
            if (d_prefixOf(d_index, *pTableIndex)) {
                myOv.id = getId(*pTableIndex);
                myOv.size = -(readsLength - j);
                ovSet.insert(myOv);
                pTableIndex++;
                nOv++;
            } else
                break;
        }

        //   <<<<<<-----
        //        I<<<<<<<<<<<I

        pTableIndex = r_prefix_lowerBound(r_index);

        while (pTableIndex <= prefixTableR + nNR_reads) {
            if (r_prefixOf(r_index, *pTableIndex)) {
                myOv.id = getId(*pTableIndex);
                myOv.size = (readsLength - j);
                ovSet.insert(myOv);
                pTableIndex++;
                nOv++;
            } else
                break;
        }

        for (set<_OV>::iterator sIt = ovSet.begin(); sIt != ovSet.end(); sIt++) {
            ID.push_back(sIt->id);
            OV.push_back(sIt->size);
        }
    }

    return nOvRight;
}

int ReadsStorage::loadReadsFiles(
        string infile1,
        string infile2,
        int mateOrientation,
        Pairing* P,
        ReadsLayout* L) {

    if (readsLength == 0) {
        LOG.oss << "[err] ReadsStorage problem\n";
        LOG.flushStream(TOSTDERR);
        sendBugReportPlease(cerr);
    }
    ifstream in1, in2;
    string buffer;
    bool paired;
    unsigned int nFastaEntries = 0;

    initDuplicateCheck();
    unsigned int nDuplicate = 0;
    bool dupInsert;

    if (infile2 == "")
        paired = false;
    else
        paired = true;

    in1.open(infile1.c_str());

    if (paired)
        in2.open(infile2.c_str());

    if (!in1) {
        LOG.oss << "[err] Cannot open the file \"" << infile1 << "\"" << endl;
        LOG.flushStream(TOSTDERR);
        return 1;
    }

    if (paired && !in2) {
        LOG.oss << "[err] Cannot open the file \"" << infile2 << "\"" << endl;
        LOG.flushStream(TOSTDERR);
        return 1;
    }

    if (paired)
        LOG.oss << "Opening paired files:\n   " << infile1 << "\n   " << infile2 << endl;
    else
        LOG.oss << "Opening file:\n   " << infile1 << endl;
    LOG.flushStream(TOSTDOUT);

    char *s1 = new char[1000];
    char *s2 = new char[1000];
    char *sComp = new char[1000];

    bool fastQ = false;

    if (paired) {
        if (in1.peek() == '>' && in2.peek() == '>')
            fastQ = false;
        else if (in1.peek() == '@' && in2.peek() == '@')
            fastQ = true;
        else {
            LOG.oss << "[err] Read file unrecognized format or paired files inconsistency\n";
            LOG.oss << "[err] Paired files must be either in FASTA or FASTQ format\n";
            LOG.flushStream(TOSTDERR);
            return 1;
        }
    } else {
        if (in1.peek() == '>')
            fastQ = false;
        else if (in1.peek() == '@')
            fastQ = true;
        else {
            LOG.oss << "[err] Read file unrecognized format\n";
            LOG.oss << "[err] Files must be either in FASTA or FASTQ format\n";
            LOG.flushStream(TOSTDERR);
            return 1;
        }
    }

    if (paired)
        P->startNewLibrary(mateOrientation);

    getline(in1, buffer); //header
    getline(in1, buffer); //read

    int originalRL = (int)buffer.length();
    if (originalRL > 0) {
        if (buffer.at(originalRL - 1) == '\r')
            originalRL--;
    }
    in1.seekg(0, ios::beg);

    bool clean1 = false, clean2 = false;
    bool dir1 = true, dir2 = true;
    unsigned int id1 = 0, id2 = 0;
    unsigned int nrId1 = 0, nrId2 = 0;
    static unsigned int nUnambiguousReads;
    unsigned int nAmbiguous = 0;
    while (!in1.eof()) {
        getline(in1, buffer); //header
        in1.read(s1, readsLength);
        getline(in1, buffer); //'\r' and/or '\n'
        if (in1.eof())
            break;

        s1[readsLength] = '\0';


        if (paired) {
            getline(in2, buffer); //header
            in2.read(s2, readsLength);
            getline(in2, buffer); //'\r' and/or '\n'
            s2[readsLength] = '\0';

        }

        clean1 = toUpperAndCheckDNA(s1);
        if (paired)
            clean2 = toUpperAndCheckDNA(s2);

        nrId1 = nrId2 = 0;

        if (clean1) {
            id1 = insertRead(s1);
            nrId1 = effectiveNReads;
            nUnambiguousReads++;
        }

        if (clean2) {
            id2 = insertRead(s2);
            nrId2 = effectiveNReads;
            nUnambiguousReads++;
        }

        dupInsert = false;

        bool removeDuplicate = false;

        if (paired && clean1 && clean2) {

            if (removeDuplicate) {

                if (isDuplicateMate(nrId1, nrId2)) {
                    nDuplicate++;
                    effectiveNReads -= 2;
                    dupInsert = true;
                }
            }

            if (!dupInsert)
                P->addPair(effectiveNReads - 1, effectiveNReads);
        }

        if (!dupInsert) {
            if (clean1) {
                //id1 = insertRead(s1);

                // nUnambiguousReads++;
                L->initLayout(nrId1, 1, true, id1);

                if (L->getLastIdentical(id1) != 0) {
                    dir1 = nrReadDirection[nrId1];
                    L->merge(L->getLastIdentical(id1), nrId1, dir1, 0);
                }
                L->setLastIdentical(id1, nrId1);
            } else
                nAmbiguous++;

            if (paired) {
                if (clean2) {
                    // id2 = insertRead(s2);
                    // nUnambiguousReads++;
                    L->initLayout(nrId2, 1, true, id2);

                    if (L->getLastIdentical(id2) != 0) {
                        dir2 = nrReadDirection[nrId2];
                        L->merge(L->getLastIdentical(id2), nrId2, dir2, 0);
                    }

                    L->setLastIdentical(id2, nrId2);
                } else
                    nAmbiguous++;
            }

            if (paired)
                nFastaEntries += 2;
            else
                nFastaEntries++;

        }

        if ((nFastaEntries % 1000) == 0) {
            if (fastQ)
                cerr << "reading fastQ entries: ";
            else
                cerr << "reading fasta entries: ";
            cerr << nFastaEntries << '\r' << flush;
        }

        if (fastQ) {
            getline(in1, buffer);
            getline(in1, buffer);
            if (paired) {
                getline(in2, buffer);
                getline(in2, buffer);
            }
        }

        if (in1.eof()) {

            if (paired && !in2.eof()) {
                cout << "   Error: paired files must exactly have the same number of reads\n";
                exit(0);
            }

            break;
        }

        if (in2.eof()) {
            cout << "   Error: paired files must exactly have the same number of reads\n";
            exit(0);
        }
    }

    if (paired)
        P->endLibrary();

    nDiscardedReads += nAmbiguous;
    LOG.oss << "   Loaded: " << nFastaEntries << ", discarded: " << nAmbiguous << endl;
    // cout << " nDuplicate: " << nDuplicate << endl;

    if (readsLength == originalRL)
        LOG.oss << "   Reads length: " << readsLength << endl;
    else
        LOG.oss << "   Reads length: " << originalRL << " truncated to " << readsLength << endl;

    LOG.flushStream(TOSTDOUT);
    
    in1.close();
    in1.clear();
    in2.close();
    in2.clear();
    delete[] s1;
    delete[] s2;
    delete[] sComp;

    return 0;
}

char** ReadsStorage::d_prefix_lowerBound(const char *r) {
    //Returns the index to the first element in the table which does not compare
    //less than r, i.e. it is either equal or greater.

    char** high, **i, **low;

    for (low = prefixTableD, high = prefixTableD + nNR_reads; high - low > 1;) {
        i = low + (high - low) / 2;

        if (strcmp(r, *i) <= 0) //r <= *i
            high = i;
        else
            low = i;
    }
    return high;
}

char** ReadsStorage::r_prefix_lowerBound(const char *r) {
    //Returns the index to the first element in the table which does not compare
    //less than r, i.e. it is either equal or greater.
    //Both r and return value point to the END of a read
    // (prefixTableR order the reads REVERSED but NOT COMPLEMENTED

    char** high, **i, **low;

    for (low = prefixTableR, high = prefixTableR + nNR_reads; high - low > 1;) {
        i = low + (high - low) / 2;

        if (r_strcmp(r, *i) <= 0) //r <= *i
            high = i;
        else
            low = i;
    }
    return high;
}



