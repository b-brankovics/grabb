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


#include "readsLayout.h"
#include "globalFunc.h"
#include "readsStorage.h"
#include "Pairing.h"
#include "crc.h"

ReadsLayout::ReadsLayout() {
    
    R = 0x0;
    lastIdentical = 0x0;
    next = 0x0;
    previous = 0x0;
    nodeId = 0x0;
    position = 0x0;
    flags = 0x0;
    tabSize = 0;
}

ReadsLayout::~ReadsLayout() {
    cleanMemory();
    //dtor
}

void ReadsLayout::init(ReadsStorage *r, unsigned int n) {
    R = r;
    tabSize = n;
    try {

        lastIdentical = (unsigned int*) malloc((tabSize) * sizeof (unsigned int));
        next = (unsigned int*) malloc((tabSize) * sizeof (unsigned int));
        previous = (unsigned int*) malloc((tabSize) * sizeof (unsigned int));
        nodeId = (unsigned int*) malloc((tabSize) * sizeof (unsigned int));
        position = (int*) malloc((tabSize) * sizeof (int));
        flags = (unsigned char*) malloc((tabSize) * sizeof (unsigned char));
        //        lastIdentical = new unsigned int[tabSize];
        //        next = new unsigned int[tabSize];
        //        previous = new unsigned int[tabSize];
        //        nodeId = new unsigned int[tabSize];
        //        position = new int[tabSize];
        //        flags = new unsigned char[tabSize];
    } catch (bad_alloc ex) {
        cout << ex.what() << "\nnot enough memory. ReadsLayout::init(STReads *r)" << '\n';
        exit(0);
    }
    memset(lastIdentical, 0, tabSize * sizeof (unsigned int));
    memset(next, 0, tabSize * sizeof (unsigned int));
    memset(previous, 0, tabSize * sizeof (unsigned int));
    memset(nodeId, 0, tabSize * sizeof (unsigned int));
    memset(position, 0, tabSize * sizeof (int));
    memset(flags,0,tabSize * sizeof (unsigned char));
}

void ReadsLayout::cleanMemory() {
    
    if (next != 0x0) {
     //   delete [] next;
        free(next);
        next = 0x0;
    }
    if (previous != 0x0) {
     //   delete [] previous;
        free(previous);
        previous = 0x0;
    }
    if (nodeId != 0x0) {
     //   delete [] nodeId;
        free(nodeId);
        nodeId = 0x0;
    }
    if (position != 0x0) {
     //   delete [] position;
        free(position);
        position = 0x0;
    }
    if (lastIdentical != 0x0) {
  //      delete [] lastIdentical;
        free(lastIdentical);
        lastIdentical = 0x0;
    }
    if (flags != 0x0) {
    //    free(flags);
        delete [] flags;
        flags = 0x0;
    }
}

void ReadsLayout::save(ostream& out_bin)
{
    unsigned int crcCheck;
    Crc32 CRC;
  
    tabSize = R->getEffectiveNReads() + 1;
    
    CRC.AddData((uint8_t*) & tabSize, sizeof (unsigned int));
    CRC.AddData((uint8_t*) lastIdentical, sizeof (unsigned int) *tabSize);
    CRC.AddData((uint8_t*) next, sizeof (unsigned int) *tabSize);
    CRC.AddData((uint8_t*) previous, sizeof (unsigned int) *tabSize);
    CRC.AddData((uint8_t*) nodeId, sizeof (unsigned int) *tabSize);
    CRC.AddData((uint8_t*) position, sizeof (int) *tabSize);
    CRC.AddData((uint8_t*) flags, sizeof (unsigned char) *tabSize);

    crcCheck=CRC.GetCrc32();
    
    out_bin.write((char*) &crcCheck, sizeof (unsigned int));
      
    out_bin.write((char*) &tabSize, sizeof (unsigned int));
    out_bin.write((char*) lastIdentical, sizeof (unsigned int) *tabSize);
    out_bin.write((char*) next, sizeof (unsigned int) *tabSize);
    out_bin.write((char*) previous, sizeof (unsigned int) *tabSize);
    out_bin.write((char*) nodeId, sizeof (unsigned int) *tabSize);
    out_bin.write((char*) position, sizeof (int) *tabSize);
    out_bin.write((char*) flags, sizeof (unsigned char) *tabSize);
}

bool ReadsLayout::load(istream& in_bin, ReadsStorage*r) {
    
  
    unsigned int crcCheck;
    Crc32 CRC;
    
    
    cleanMemory();
    unsigned int tabsize;

    in_bin.read((char*) &crcCheck, sizeof (unsigned int));
    
    in_bin.read((char*) &tabsize, sizeof (unsigned int));

    init(r, tabsize);

    in_bin.read((char*) lastIdentical, sizeof (unsigned int) *tabSize);
    in_bin.read((char*) next, sizeof (unsigned int) *tabSize);
    in_bin.read((char*) previous, sizeof (unsigned int) *tabSize);
    in_bin.read((char*) nodeId, sizeof (unsigned int) *tabSize);
    in_bin.read((char*) position, sizeof (int) *tabSize);
    in_bin.read((char*) flags, sizeof (unsigned char) *tabSize);
    
    CRC.AddData((uint8_t*) & tabSize, sizeof (unsigned int));
    CRC.AddData((uint8_t*) lastIdentical, sizeof (unsigned int) *tabSize);
    CRC.AddData((uint8_t*) next, sizeof (unsigned int) *tabSize);
    CRC.AddData((uint8_t*) previous, sizeof (unsigned int) *tabSize);
    CRC.AddData((uint8_t*) nodeId, sizeof (unsigned int) *tabSize);
    CRC.AddData((uint8_t*) position, sizeof (int) *tabSize);
    CRC.AddData((uint8_t*) flags, sizeof (unsigned char) *tabSize);
    
    unsigned int crc2=CRC.GetCrc32();
    
    if (crc2 != crcCheck)
        return false;
    
    return true;
}

void ReadsLayout::initLayout(size_t layout, int pos, bool dir, unsigned int nodeId) {
    
    if (getNext(layout) != 0 || getPrevious(layout) != 0) {
        cout << "void ReadsLayout::initLayout(size_t layout, int pos, bool dir, unsigned int nodeId) problem\n";
        cout << "layout=" << layout << " next=" << getNext(layout) << " previous=" << getPrevious(layout) << endl;
        cout << "layout-1=" << layout-1 << " next=" << getNext(layout-1) << " previous=" << getPrevious(layout-1) << endl;
        cout << "layout+1=" << layout+1 << " next=" << getNext(layout+1) << " previous=" << getPrevious(layout+1) << endl;
        sendBugReportPlease(cerr);
    }

    setPosition(layout, pos, dir);
    setNodeId(layout, nodeId);
}

void ReadsLayout::setLayoutNodeId(size_t layout, unsigned int nodeId) {
    
    if (getNext(layout) != 0) {
        cout << "void ReadsLayout::setLayoutNodeId(size_t layout, unsigned int nodeId) problem\n";
        sendBugReportPlease(cerr);
    }

    do {
        setNodeId(layout, nodeId);
        layout = getPrevious(layout);

    } while (layout != 0);
}

string ReadsLayout::getDirectRead(size_t i) const {
    
    char s[1000];
    R->getDirectSequence(s,i);
    return s;
}

string ReadsLayout::getReverseRead(size_t i) const {
    
    char s[1000];
    R->getReverseSequence(s,i);
    return s;
}

void ReadsLayout::getPCharRead(char *dest ,size_t i, bool dir) const {

    if (dir)
        R->getDirectSequence(dest,i);
    else
        R->getReverseSequence(dest,i);
}

size_t ReadsLayout::getBegin(size_t index) const {
    
    while (getPrevious(index) != 0) {
        index = getPrevious(index);
    }

    return index;
}

size_t ReadsLayout::getEnd(size_t index) const {
    
    while (getNext(index) != 0) {
        index = getNext(index);
    }

    return index;
}

size_t ReadsLayout::reverseComplement(size_t index) {
    
    if (getNext(index) != 0) {
        cout << "size_t ReadsLayout::reverseComplement(size_t index) problem\n";
        sendBugReportPlease(cerr);
    }

    unsigned int reversePos = getPosition(index) + 1;
    size_t pTmp;
    size_t previous;

    do {
        setDirection(index, !getDirection(index));
        setPosition(index, reversePos - getPosition(index));
        //swap links
        pTmp = getNext(index);
        setNext(index, getPrevious(index));
        setPrevious(index, pTmp);
        previous = index;
        index = getNext(index);

    } while (index != 0);

    return previous;
}

//merge relative to l1:
//example: direct=false, shift=-3
// <<<<<<<<<<<<<
//    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// (shift is zero based)

size_t ReadsLayout::merge(size_t l1, size_t l2, bool direct, int shift) {
    
    if (getNext(l1) != 0 || getNext(l2) != 0 || l1 == 0 || l2 == 0) {
        cout << "size_t ReadsLayout::merge(size_t l1, size_t l2, bool direct, int shift) problem\n";
        sendBugReportPlease(cerr);
    }
    size_t pTmp;

    if (direct == false)
        l2 = reverseComplement(l2);

    if (shift < 0) {
        shift = -shift;
        pTmp = l1;
        l1 = l2;
        l2 = pTmp;
    }

    pTmp = l2;

    while (pTmp != 0) {
        setPosition(pTmp, getPosition(pTmp) + shift);
        pTmp = getPrevious(pTmp);
    }

    size_t p1 = l1;
    size_t p2 = l2;


    while (getPosition(getPrevious(p2)) >= getPosition(p1)) {
        p2 = getPrevious(p2);
    }

    if (getPosition(p2) >= getPosition(p1)) //add the overhanging part
    {
        pTmp = getPrevious(p2);
        setPrevious(p2, p1);
        setNext(p1, p2);

        p2 = pTmp;
        l1 = l2;
    }

    while (p2 != 0) {
        if (getPosition(p2) >= getPosition(p1)) {
            //insert p2 after p1
            pTmp = getPrevious(p2);

            setPrevious(p2, p1);
            setNext(p2, getNext(p1));
            setNext(p1, p2);
            setPrevious(getNext(p2), p2);

            p1 = p2;
            p2 = pTmp;
        } else
            p1 = getPrevious(p1);

    }

    return l1;
}

bool ReadsLayout::isUniDirectional(size_t index) const {
    size_t p = getBegin(index);
    bool direction = getDirection(p);

    do {
        if (getDirection(p) != direction)
            return false;

        p = getNext(p);

    } while (p != 0);

    return true;
}

unsigned int ReadsLayout::getSequenceLength(size_t listIndex) const {
    
    return getPosition(listIndex) + R->getReadsLength() - 1;
}

unsigned int ReadsLayout::getNReads(size_t listIndex) const {
    //return the total number of bases in the layout
    if (getNext(listIndex) != 0) {
        cout << "size_t ReadsLayout::getCoverage(size_t listIndex) problem\n";
        sendBugReportPlease(cerr);
    }

    unsigned int count = 1;

    while (getPrevious(listIndex) != 0) {
        count++;
        listIndex = getPrevious(listIndex);
    }

    return count;

}

unsigned int ReadsLayout::getNext(size_t i) const {
    return next[i];
}

int ReadsLayout::getPosition(size_t i) const {
    return position[i] > 0 ? position[i] : -position[i];
}

bool ReadsLayout::getDirection(size_t i) const {
    return position[i] > 0 ? true : false;
}

unsigned int ReadsLayout::getNodeId(size_t i) const {
    return nodeId[i];
}

void ReadsLayout::setNext(size_t i, unsigned int id) {
    next[i] = id;
}

void ReadsLayout::setPrevious(size_t i, unsigned int id) {
    previous[i] = id;
}

void ReadsLayout::setPosition(size_t i, int pos) {
    position[i] > 0 ? position[i] = pos : position[i] = -pos;
}

void ReadsLayout::setDirection(size_t i, bool dir) {
    (dir == (position[i] > 0)) ? position[i] = position[i] : position[i] = -position[i];
}

void ReadsLayout::setPosition(size_t i, int pos, bool dir) {
    dir ? position[i] = pos : position[i] = -pos;
}

void ReadsLayout::setNodeId(size_t i, unsigned int id) {
    nodeId[i] = id;
}

void ReadsLayout::unChain(size_t i) {
    next[i] = 0;
    previous[i] = 0;
}

unsigned int ReadsLayout::getLastIdentical(size_t i) const {
    return lastIdentical[i];
}

void ReadsLayout::setLastIdentical(size_t i, unsigned int lastId) {
    lastIdentical[i] = lastId;
}

string ReadsLayout::getSequence(size_t listIndex) {
    //basic sequence construction, assume ordered position

    if (getNext(listIndex) != 0) {
        cout << "string ReadsLayout::getSequence(size_t listIndex) problem\n";
        sendBugReportPlease(cerr);
    }
    string s;
    static char s2[1000];
    unsigned int pos = 0;
    unsigned int ov;
    size_t index = getBegin(listIndex);

    if (getDirection(index))
        s = getDirectRead(index);
    else
        s = getReverseRead(index);

    index = getNext(index);
    pos = R->getReadsLength();
    while (index != 0) {
        ov = pos - getPosition(index) + 1;
        
      //  getSequence(s2,index,getDirectionIndex,ov);
        
        getPCharRead(s2,index,getDirection(index));
        s += s2+ov;//(p+ov);
        
//        if (getDirection(index))
//            s += getDirectRead(index).substr(ov);
//        else
//            s += getReverseRead(index).substr(ov);

        pos += R->getReadsLength() - ov;
        index = getNext(index);
    }

    return s;
}

string ReadsLayout::getReverseSequence(size_t listIndex) {
    string s;
    getCompReverse(getSequence(listIndex), s);
    return s;
}

void ReadsLayout::getVcoverage(size_t listIndex, bool direction, vector<unsigned int> &cov)
{
    unsigned int pos;
    unsigned int sequenceLength=getSequenceLength(listIndex);
    int rl=R->getReadsLength();
    size_t index;
    
    cov.clear();
    cov.assign(sequenceLength,0);
    vector<unsigned int> tmp;
    tmp.assign(sequenceLength,0);

    if (direction)
    {
        index = getBegin(listIndex);

        while (index != 0)
        {
            pos = getPosition(index);
            cov[pos - 1]++;
            index = getNext(index);
        }
    }
    else
    {
        index = listIndex;
        sequenceLength-= (rl-1);

        while (index != 0)
        {
            pos = getPosition(index);
            pos=sequenceLength-pos;
            cov[pos]++;
            index = getPrevious(index);
        }
    }
    
    tmp[0]=cov[0];
   
    for (size_t i=1; i<cov.size(); i++)
    {
        cov[i]=cov[i]+cov[i-1];
        tmp[i]=cov[i];
    }
    
    for (size_t i=rl; i<cov.size(); i++)
        cov[i]-=tmp[i-rl];
}

void ReadsLayout::print(size_t index, ostream &out, bool dir, unsigned int start, unsigned int maxD, Pairing *P) {

    if (getNext(index) != 0) {
        cerr << "void ReadsLayout::print(size_t index) problem\n";
        sendBugReportPlease(cerr);
    }
    if (!dir)
        index = reverseComplement(index);

    size_t p = getBegin(index);
    size_t tmp;

    do {
        
        unsigned int position=getPosition(p);
        
        if (position > maxD)
            break;

        if (position < start)
        {
            tmp = p;
            p = getNext(p);
            continue;
        }

        unsigned int pairedRead=0;
        unsigned int pairedNode=0;
        int lib=0;

        if (P->getNLibrary() != 0)
        {
            pairedRead = P->getPairing(p);
            pairedNode = getNodeId(pairedRead);
            lib = P->getPeLibraryID(p);
        }
        
        if (getDirection(p))
            out << '>';
        else
            out << '<';

        for (int i = 0; i < getPosition(p) % 120; i++)
            out << " ";
        if (getDirection(p))
            out << getDirectRead(p) << " " << p << ' ' << lib << ' ' << pairedNode << '\n';
        else
            out << getReverseRead(p) << " " << p << ' ' << lib << ' ' << pairedNode << '\n';
        
        tmp = p;
        p = getNext(p);

    } while (tmp != index);

    out << flush;

    if (!dir) //back to initial direction
        index = reverseComplement(index);
}

void ReadsLayout::statOverlaps(size_t listIndex, double &s, double &ss, unsigned int * distr) {

    if (getNext(listIndex) != 0) {
        cout << "string ReadsLayout::statOverlaps(size_t i, double &s,double &ss) problem\n";
        sendBugReportPlease(cerr);
    }

    s = ss = 0.0;

    unsigned int pos = 0;
    unsigned int ov;
    size_t index = getBegin(listIndex);

    index = getNext(index);
    pos = R->getReadsLength();

    while (index != 0) {
        ov = pos - getPosition(index) + 1;
        distr[ov]++;
        s += ov;
        ss += (ov * ov);
        pos += R->getReadsLength() - ov;
        index = getNext(index);
    }
}

void ReadsLayout::statOverlaps2(size_t listIndex,
                            unsigned int &nOverlap,
                            unsigned int &nSample,
                            double &mean,
                            unsigned int* distr)
{
    unsigned int pos = 0;
    unsigned int oh; //over hanging
    unsigned int mult;
    size_t index = getBegin(listIndex);
    
    mean=0.0;
    nOverlap=nSample=0;

    index = getNext(index);
    pos = R->getReadsLength();
    unsigned int rl=R->getReadsLength();

    while (index != 0) {
        
        mult=1;
        oh = rl - (pos - getPosition(index) + 1);
        
        while(oh==0)
        {
            mult++;
            index = getNext(index);
            if (index == 0)
                break;
            pos += oh;
            oh = rl - (pos - getPosition(index) + 1);
        }

        if (mult > 1)
        {
            nOverlap+=mult-1;
            mult = (mult * mult) - mult;
            distr[0] += mult;
            nSample += mult;
            if (index == 0)
                break;
        }
        
        distr[oh]++;
        mean+=oh;
        nSample++;
        nOverlap++;
        pos += oh;
        index = getNext(index);
    }
    if (nSample==0)
        mean=0.0;
    else
        mean/=nSample;
}

unsigned int ReadsLayout::sampleOH(size_t listIndex,
                               bool dir,
                               unsigned int maxD,
                               unsigned int maxN,
                               unsigned int *distr)
{
    unsigned int index;
    unsigned int pos;
    int nodeLength=getSequenceLength(listIndex);
    unsigned int rl=R->getReadsLength();
    unsigned int nOverlap=0;
    unsigned int mult;
    int oh; //overhang

    if (dir)
    {
        index = getBegin(listIndex);
        pos = 1;
        index = getNext(index);

        while (index != 0)
        {
            if (pos > maxD && maxD!=0)
                break;
            if (nOverlap > maxN && nOverlap!=0)
                break;
            
            mult = 1;
            oh = getPosition(index) - pos;

            while (oh == 0)
            {
                mult++;
                index = getNext(index);
                if (index == 0)
                    break;

                oh = getPosition(index) - pos;
            }

            if (mult > 1)
            {
                nOverlap += mult - 1;
                mult = (mult * mult) - mult;
                distr[0] += mult;

                if (index == 0)
                    break;
            }

            distr[oh]++;
            nOverlap++;
            pos += oh;
            index = getNext(index);
        }
    }
    else
    {
        pos=nodeLength-rl+1;
        index = getPrevious(listIndex);

        while (index != 0)
        {
            if (nodeLength - (pos+rl-1) > maxD && maxD!=0)
                break;
            if (nOverlap > maxN && nOverlap!=0)
                break;
            
            mult = 1;
            oh = pos - getPosition(index);

            while (oh == 0)
            {
                mult++;
                index = getPrevious(index);
                if (index == 0)
                    break;

                oh = pos - getPosition(index);
            }

            if (mult > 1)
            {
                nOverlap += mult - 1;
                mult = (mult * mult) - mult;
                distr[0] += mult;

                if (index == 0)
                    break;
            }

            distr[oh]++;
            nOverlap++;
            pos -= oh;
            index = getPrevious(index);
        }
    }
    
    return nOverlap;
}

double ReadsLayout::getMeanOverlapSize(size_t listIndex)
{
    double sum=0.0;
    size_t index = getBegin(listIndex);
    size_t pos=R->getReadsLength();
    unsigned int ov;
    unsigned int nEchant=0;
    
    index = getNext(index);
    
    while (index != 0) {
        ov=pos - getPosition(index) + 1;
        sum+=ov;
        pos += R->getReadsLength() - ov;
        index = getNext(index);
        nEchant++;
    }
    
    if (sum == 0)
        return 0.0;
    else
        return sum/nEchant;
}


bool ReadsLayout::checkLayout(size_t index) {
    if (getNext(index) != 0) {
        cerr << "next!=0" << endl;
        return false;
    }

    int lastPos = getPosition(index);
    lastIdentical[index] = 1000000000; //used as a flag

    while (getPrevious(index) != 0) {
        index = getPrevious(index);

        if (lastIdentical[index] == 1000000000) {
            cerr << "cycle!" << endl;
            return false;
        } else
            lastIdentical[index] = 1000000000;

        if (getPosition(index) > lastPos) {
            cerr << "posSorting!" << endl;
            return false;
        } else
            lastPos = getPosition(index);
    }

    return index;
}

void ReadsLayout::writeReadsAsFasta(size_t index, ostream &out)
{
    while (index!=0)
    {
        out << '@' << index << '\n';
        out << getDirectRead(index) << endl;
        out << '+' << index << '\n';
        for (int i=0; i<54; i++) //fastq
            out << 'C';
        out << '\n';
        index=getPrevious(index);
    }
}

void ReadsLayout::flagReads(size_t index, char state)
{
    while (index!=0)
    {
        R->setFlag(index,state);
        index=getPrevious(index);
    }
}