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

#include <algorithm>
#include <iomanip>

#include "stat.h"
#include "globalFunc.h"

void stats(vector<int>::iterator start, vector<int>::iterator end, ostream &out) {
    int nEchant = end - start;
    if (nEchant == 0)
        return;

    int minEchant = *start, maxEchant = *start;
    double mean = 0;
    unsigned int sum = 0, N50 = 0;
    unsigned int cumul = 0;

    sort(start, end);

    for (vector<int>::iterator iIt = start; iIt != end; iIt++) {

        if (*iIt < minEchant)
            minEchant = *iIt;
        if (*iIt > maxEchant)
            maxEchant = *iIt;

        sum += *iIt;
    }

    mean = (double) sum / nEchant;

    for (vector<int>::iterator iIt = start; iIt != end; iIt++) {
        cumul += *iIt;
        if (cumul > sum / 2) {
            N50 = *iIt;
            break;
        }
    }

    out << "   sum:  " << smartDNALength(sum) << endl;
    out << "   N50:  " << smartDNALength(N50) << endl;
    out << "   mean: " << smartDNALength(mean) << endl;
    out << "   max:  " << smartDNALength(maxEchant) << endl;
    out << "   min:  " << smartDNALength(minEchant) << endl;
}

//geometric distribution of overHanging sizes
//p  probability (=1/E(X))
//k  overhanging value
//c  correction factor

double cdfOH(double p, unsigned int k, double c) {
    if (k == 0)
        return 1.0;
    else
        return c * pow(1 - p, (int) k);
}

double pdfOH(double p, unsigned int k, double c) {
    if (k == 0)
        return p * c - c + 1;
    else
        return p * pow((1 - p), (int) k) * c;
}