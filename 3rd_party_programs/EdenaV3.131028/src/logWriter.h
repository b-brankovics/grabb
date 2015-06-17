/* 
 * File:   logWriter.h
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

#ifndef LOGWRITER_H
#define	LOGWRITER_H

#define TOSTDOUT 1
#define TOSTDERR 2

#include <fstream>
#include <string>
#include <sstream>

using namespace std;

class logWriter
{
public:
    logWriter();
    logWriter(const logWriter& orig);
    virtual ~logWriter();
    
    void write(string message, int flag);
    inline void write(string message){write(message,0);};
    void writeRaw(string message);
    void flushStream(int flag);
    void flushRawStream();
    inline void flushStream(){flushStream(0);};
    bool open(string fileName, bool overwrite);
    void close();
   
    ostringstream oss;
    
private:
    
    const char* localTime();
    fstream out;
  
};

#endif	/* LOGWRITER_H */

