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

#ifndef STAT_H_INCLUDED
#define STAT_H_INCLUDED

#include <vector>
#include <ostream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstring>
using namespace std;

void stats(vector<int>::iterator start, vector<int>::iterator end, ostream &out);
double cdfOH(double p, unsigned int k, double c); 
double pdfOH(double p, unsigned int k, double c);

#endif // STAT_H_INCLUDED
