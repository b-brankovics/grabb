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

#ifndef GLOBALFUNC_H_INCLUDED
#define GLOBALFUNC_H_INCLUDED

#define MAXLOG 10000

#include <iostream>
#include <string>
using namespace std;

string edenaVersion();
void edenaVersion(ostream &out);
void edenaAuthors(ostream &out);
void edenaUsage(int exitValue);
void sendBugReportPlease(ostream &out);

void getCompReverse(string source, string &comp);

char compReverse(char);
void compReverse(const char* source, char* dest, size_t length);
void compReverse(string source, string &comp);

void lineWrap(ostream&, const string&, int);
double fastLog(unsigned int v);
unsigned int hamming(const char* a, const char* b);
bool isNotAmbiguous(char*s, unsigned int length);
unsigned int estimateNReads(string fileName, int &rl);
string smartDNALength(double l);
void getline_c(FILE *p, string &line);

#endif // GLOBALFUNC_H_INCLUDED
