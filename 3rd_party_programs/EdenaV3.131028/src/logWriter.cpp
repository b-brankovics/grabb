/* 
 * File:   logWriter.cpp
 * Author: david
 * 
 * Created on February 21, 2013, 3:55 PM
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

#include "logWriter.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
using namespace std;

logWriter::logWriter() {

}

logWriter::logWriter(const logWriter& orig) {
}

logWriter::~logWriter() {
    close();
}

void logWriter::write(string message, int flag) {

    if (message == "")
        return;

    if (*(message.end() - 1) != '\n')
        message += '\n';

    if (flag == TOSTDOUT) {
        cout << message;
        cout << flush;
    } else if (flag == TOSTDERR) {
        cerr << message;
        cerr << flush;
    }

    vector<string> lines;
    string buffer;
    istringstream iss;
    iss.str(message);
    while (!iss.fail()) {
        getline(iss, buffer);

        while (buffer.size() > 0) {
            if (buffer[0] == ' ')
                buffer = buffer.substr(1);
            else
                break;
        }
        if (buffer.size() > 0)
            lines.push_back(buffer);
    }
    
    if (out) {
        for (size_t i=0; i<lines.size(); i++)
        {
            out << localTime() << " ";
            out << lines[i] << '\n';
        }
        out << flush;
    }

//
//    if (out) {
//        out << localTime() << " --- ";
//        out << message;
//        out << flush;
//    }
}

void logWriter::writeRaw(string message)
{
    if (out)
    {
        out << message;
    }
}

void logWriter::flushStream(int flag) {
    string mes = oss.str();
    write(mes, flag);
    oss.clear();
    oss.str("");
}

void logWriter::flushRawStream()
{
    if (out)
    {
        out << oss.str();
    }
    oss.clear();
    oss.str("");
}

bool logWriter::open(string fileName, bool overwrite) {

    if (overwrite) {
        out.open(fileName.c_str(), fstream::out);
    } else
        out.open(fileName.c_str(), fstream::out | fstream::app);
    
    oss << setprecision(1) << fixed;

    return out; //T/F
}

void logWriter::close() {
    if (out.is_open()) {
        out.close();
        out.clear();
    }
}

const char* logWriter::localTime() {

    time_t rawtime;
    struct tm * timeinfo;
    static char buffer[80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer, 80, "[%m/%d/%Y:%H:%M:%S]", timeinfo);
    return buffer;
}