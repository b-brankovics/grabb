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

#include "Pairing.h"
#include "globalFunc.h"
#include "logWriter.h"
#include "crc.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
using namespace std;

extern logWriter LOG;
int PELibrary::readLength;

PELibrary::PELibrary() {
    // mateOrientation=0; //direct-reverse by default
    startIndex = 0;
    endIndex = 0;
    meanCloneSize = 0.0;
    sdCloneSize = 0.0;
    nUsableMates = 0;
    expectedMateCoverage = 0.0;
}

PELibrary::~PELibrary() {

}

void PELibrary::save(ostream& out_bin) {
    Crc32 CRC;
    unsigned int crcCheck;

    CRC.AddData((uint8_t*) & mateOrientation, sizeof (int));
    CRC.AddData((uint8_t*) & startIndex, sizeof (size_t));
    CRC.AddData((uint8_t*) & endIndex, sizeof (size_t));
    CRC.AddData((uint8_t*) & meanCloneSize, sizeof (double));
    CRC.AddData((uint8_t*) & sdCloneSize, sizeof (double));
    CRC.AddData((uint8_t*) & distanceMin, sizeof (unsigned int));
    CRC.AddData((uint8_t*) & distanceMax, sizeof (unsigned int));

    crcCheck = CRC.GetCrc32();
    out_bin.write((char*) &crcCheck, sizeof (unsigned int));

    out_bin.write((char*) &mateOrientation, sizeof (int));
    out_bin.write((char*) &startIndex, sizeof (size_t));
    out_bin.write((char*) &endIndex, sizeof (size_t));
    out_bin.write((char*) &meanCloneSize, sizeof (double));
    out_bin.write((char*) &sdCloneSize, sizeof (double));
    out_bin.write((char*) &distanceMin, sizeof (unsigned int));
    out_bin.write((char*) &distanceMax, sizeof (unsigned int));
}

bool PELibrary::load(istream& in_bin) {
    Crc32 CRC;
    unsigned int crcCheck;

    in_bin.read((char*) &crcCheck, sizeof (unsigned int));

    in_bin.read((char*) &mateOrientation, sizeof (int));
    in_bin.read((char*) &startIndex, sizeof (size_t));
    in_bin.read((char*) &endIndex, sizeof (size_t));
    in_bin.read((char*) &meanCloneSize, sizeof (double));
    in_bin.read((char*) &sdCloneSize, sizeof (double));
    in_bin.read((char*) &distanceMin, sizeof (unsigned int));
    in_bin.read((char*) &distanceMax, sizeof (unsigned int));

    CRC.AddData((uint8_t*) & mateOrientation, sizeof (int));
    CRC.AddData((uint8_t*) & startIndex, sizeof (size_t));
    CRC.AddData((uint8_t*) & endIndex, sizeof (size_t));
    CRC.AddData((uint8_t*) & meanCloneSize, sizeof (double));
    CRC.AddData((uint8_t*) & sdCloneSize, sizeof (double));
    CRC.AddData((uint8_t*) & distanceMin, sizeof (unsigned int));
    CRC.AddData((uint8_t*) & distanceMax, sizeof (unsigned int));

    unsigned int crc = CRC.GetCrc32();
    if (crc != crcCheck)
        return false;

    return true;
}

Pairing::Pairing() {
    pairingIndex = 0x0;
    fastPairing = 0x0;
    libraryId = 0x0;
    nDRLib = 0;
}

Pairing::~Pairing() {
    cleanMemory();
}

void Pairing::init() {
    cleanMemory();
    R1.clear();
    R2.clear();
    vPeLibrary.clear();
    pairedEndsMaxD = 0;
    matePairsMaxD = 0;
    nDRLib = 0;
}

void Pairing::cleanMemory() {
    if (pairingIndex != 0x0) {
        delete [] pairingIndex;
        pairingIndex = 0x0;
    }

    if (fastPairing != 0x0) {
        delete [] fastPairing;
        fastPairing = 0x0;
    }

    if (libraryId != 0x0) {
        delete [] libraryId;
        libraryId = 0x0;
    }
}

unsigned int Pairing::getNPairing() const {
    return R1.size();
}

unsigned int Pairing::getR1(size_t index) const {
    return R1.at(index);
}

unsigned int Pairing::getR2(size_t index) const {
    return R2.at(index);
}

unsigned int Pairing::getPairing(unsigned int readId, bool &pair2) const {
    if (readId < 1 || readId > numberOfReads)
        return 0;

    unsigned int index;
    index = pairingIndex[readId];

    if (index == 0)
        return 0;

    if (R1[index] == readId) {
        pair2 = true;
        return R2[index];
    }

    if (R2[index] == readId) {
        pair2 = false;
        return R1[index];
    }

    cout << "unsigned int Pairing::getPairing(...) problem\n";
    sendBugReportPlease(cout);

    return 0;
};

void Pairing::getDistanceRange(unsigned int readId, unsigned int &min, unsigned int &max, int &mateOrientation) {
    size_t peLib = libraryId[readId];

    if (peLib == 0) {
        min = max = 0;
        return;
    }
    peLib--;
    min = vPeLibrary.at(peLib).getDistanceMin();
    max = vPeLibrary.at(peLib).getDistanceMax();
    mateOrientation = vPeLibrary.at(peLib).getMateOrientation();

    return;
}

void Pairing::updatePERange(double nsd) {
    vector<PELibrary>::iterator it;
    int min, max;
    double tolerance;
    vector<unsigned int>::iterator distrIt;

    setGlobalMaxAllowedDistance(0);
    LOG.oss << "   [allowed distance range] (mean,sd) (tot,valid,sampled,usable)\n";
    for (size_t i = 0; i < vPeLibrary.size(); i++) {
        it = getLibraryIt(i);
        tolerance = nsd * it->getSdCloneSize();
        min = (int) (it->getMeanCloneSize() - tolerance);
        if (min < 0)
            min = 0;
        max = (int) (it->getMeanCloneSize() + tolerance);
        if (max > (int) getGlobalMaxAllowedDistance())
            setGlobalMaxAllowedDistance(max);
        it->setDistanceMin(min);
        it->setDistanceMax(max);

        LOG.oss << "   lib" << i + 1 << " ";
        if (it->getMateOrientation() == 1)
            LOG.oss << "><";
        else if (it->getMateOrientation() == 2)
            LOG.oss << "<>";

        LOG.oss << " [" << min << "," << max << "] ("
                << setprecision(3)
                << (it->getMeanCloneSize())
                << "," << it->getSdCloneSize() << ") ("
                << it->getSize()
                << "," << countNValidPair(i)
                << "," << it->getNEchant()
                << "," << it->getNUsableMates()
                << ")" << endl;

        LOG.flushStream(TOSTDOUT);

        if (it->getMateOrientation() == 1) {
            if (max > (int)pairedEndsMaxD)
                pairedEndsMaxD = max;
        }
        else if (it->getMateOrientation() == 2) {
            if (max > (int)matePairsMaxD)
                matePairsMaxD = max;
        }
    }
}

void Pairing::addPair(unsigned int r1, unsigned int r2) {
    R1.push_back(r1);
    R2.push_back(r2);
}

void Pairing::startNewLibrary(int mateOrientation) {
    if (vPeLibrary.size() == 255) //should never happen
    {
        cerr << "void Pairing::startNewLibrary() problem" << endl;
        exit(0);
    }

    PELibrary myLib;
    myLib.setStartIndex(R1.size());
    myLib.setMateOrientation(mateOrientation);
    vPeLibrary.push_back(myLib);

    if (mateOrientation == 1)
        nDRLib++;
}

void Pairing::endLibrary() {
    (vPeLibrary.end() - 1)->setEndIndex(R1.size() - 1);
}

void Pairing::buildIndex(unsigned int n) {
    if (R1.size() == 0)
        return;

    numberOfReads = n;
    cleanMemory();
    try {
        pairingIndex = new unsigned int [numberOfReads + 1];
        fastPairing = new unsigned int [numberOfReads + 1];
        libraryId = new unsigned char [numberOfReads + 1];
    }    catch (bad_alloc ex) {
        cout << ex.what() << "\not enough memory. Pairing::buildIndex(...)" << '\n';
        exit(0);
    }
    memset(pairingIndex, 0, (numberOfReads + 1) * sizeof (unsigned int));
    memset(fastPairing, 0, (numberOfReads + 1) * sizeof (unsigned int));
    memset(libraryId, 0, (numberOfReads + 1) * sizeof (unsigned char));

    for (unsigned int i = 0; i < R1.size(); i++) {
        pairingIndex[R1.at(i)] = i;
        pairingIndex[R2.at(i)] = i;

        fastPairing[R1.at(i)] = R2.at(i);
        fastPairing[R2.at(i)] = R1.at(i);
    }

    if (vPeLibrary.at(0).getStartIndex() != 0) {
        cerr << "void Pairing::buildIndex(unsigned int n) problem" << endl;
        exit(0);
    }

    updatePeLibraryIndex();

}

void Pairing::updatePeLibraryIndex() {
    size_t start;
    size_t end;

    nDRLib = 0;

    for (unsigned int i = 0; i < vPeLibrary.size(); i++) {

        if (vPeLibrary.at(i).mateOrientation == 1)
            nDRLib++;

        start = vPeLibrary.at(i).getStartIndex();
        end = vPeLibrary.at(i).getEndIndex();

        for (size_t j = start; j <= end; j++) {
            libraryId[R1.at(j)] = i + 1;
            libraryId[R2.at(j)] = i + 1;
        }
    }

}

void Pairing::unPair(size_t lib) {

    size_t start;
    size_t end;

    start = vPeLibrary.at(lib).getStartIndex();
    end = vPeLibrary.at(lib).getEndIndex();

    for (size_t j = start; j <= end; j++) {
        libraryId[R1.at(j)] = 0;
        libraryId[R2.at(j)] = 0;
    }

}

unsigned int Pairing::countNValidPair(size_t lib) {
    size_t start;
    size_t end;

    start = vPeLibrary.at(lib).getStartIndex();
    end = vPeLibrary.at(lib).getEndIndex();
    unsigned int nValid = end - start + 1;

    for (size_t j = start; j <= end; j++) {

        if (R->getFlag(R1.at(j)) == 1) {
            nValid--;
            continue;
        }

        if (R->getFlag(R2.at(j)) == 1)
            nValid--;
    }
    return nValid;
}

void Pairing::save(ostream& out_bin) {
    unsigned int crcCheck;

    unsigned int uint = getNPairing();

    Crc32 CRC;

    CRC.AddData((uint8_t*) & numberOfReads, sizeof (unsigned int));
    CRC.AddData((uint8_t*) & uint, sizeof (unsigned int));
    CRC.AddData((uint8_t*)&(*R1.begin()), R1.size() * sizeof (unsigned int));
    CRC.AddData((uint8_t*)&(*R2.begin()), R1.size() * sizeof (unsigned int));

    crcCheck = CRC.GetCrc32();
    out_bin.write((char*) &crcCheck, sizeof (unsigned int));

    out_bin.write((char*) &numberOfReads, sizeof (unsigned int));
    out_bin.write((char*) &uint, sizeof (unsigned int));

    //this suppose the element in the vector to be stored contiguously in memory
    //this should be actually the case but not in the spec. of the STL
    out_bin.write((char*) &(*R1.begin()), R1.size() * sizeof (unsigned int));
    out_bin.write((char*) &(*R2.begin()), R2.size() * sizeof (unsigned int));

    uint = vPeLibrary.size();
    out_bin.write((char*) &uint, sizeof (unsigned int));

    for (unsigned int i = 0; i < vPeLibrary.size(); i++)
        vPeLibrary.at(i).save(out_bin);
}

bool Pairing::load(istream& in_bin) {
    unsigned int uint;
    unsigned int crcCheck;

    init();

    in_bin.read((char*) &crcCheck, sizeof (unsigned int));
    in_bin.read((char*) &numberOfReads, sizeof (unsigned int));
    in_bin.read((char*) &uint, sizeof (unsigned int));

    R1.assign(uint, 0);
    R2.assign(uint, 0);
    //this suppose the element in the vector to be stored contiguously in memory
    //this should be actually the case but not in the spec. of the STL
    in_bin.read((char*) &(*R1.begin()), R1.size() * sizeof (unsigned int));
    in_bin.read((char*) &(*R2.begin()), R2.size() * sizeof (unsigned int));

    Crc32 CRC;

    CRC.AddData((uint8_t*) & numberOfReads, sizeof (unsigned int));
    CRC.AddData((uint8_t*) & uint, sizeof (unsigned int));
    CRC.AddData((uint8_t*)&(*R1.begin()), R1.size() * sizeof (unsigned int));
    CRC.AddData((uint8_t*)&(*R2.begin()), R1.size() * sizeof (unsigned int));

    unsigned int crc = CRC.GetCrc32();

    if (crc != crcCheck)
        return false;

    in_bin.read((char*) &uint, sizeof (unsigned int));
    PELibrary peLib;
    vPeLibrary.clear();

    for (unsigned int i = 0; i < uint; i++) {
        if (peLib.load(in_bin) == false)
            return false;
        vPeLibrary.push_back(peLib);
    }

    buildIndex(numberOfReads);

    return true;
}

