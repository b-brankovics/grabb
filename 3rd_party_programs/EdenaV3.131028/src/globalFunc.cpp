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

#define _FILE_OFFSET_BITS 64

#include "globalFunc.h"
#include "logWriter.h"
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <sys/types.h>
#include <unistd.h>

extern logWriter LOG;

string edenaVersion() {
    return "Edena v3.131028";
}

void edenaVersion(ostream &out) {
    out << edenaVersion() << endl;
}

void edenaAuthors(ostream &out) {
    out << "Copyright (C) 2008,2011,2012,2013\nDavid Hernandez, Patrice Francois, Jacques Schrenzel\n";
    out << "Genomic Research Laboratory, Geneva University Hospitals, Switzerland\n";
    out << "All rights reserved.\n";
}

void edenaUsage(int exitValue) {
    edenaVersion(cout);
    edenaAuthors(cout);

    cout << "\n";
    cout << "PROGRAM OPTIONS:\n"
            << "  1) Overlapping mode:\n"
            << "    -nThreads <int>         Number of threads to use. Default 2\n"
            << "    -r\n"
            << "    -singleEnd <file1> <file2> ...\n"
            << "                            Reads file(s) in FASTA or FASTQ format.\n"
            << "                            Several files can be specified\n"
            << "    -DRpairs\n"
            << "    -paired <file1_1> <file1_2> <file2_1> <file_2_2> ...\n"
            << "                            Direct-Reverse paired reads files. Several\n"
            << "                            pairs of files can be specified.\n"
            << "    -RDpairs\n"
            << "    -matePairs <file1_1> <file1_2> <file2_1> <file_2_2> ...\n"
            << "                            Reverse-Direct paired reads files. Several\n"
            << "                            pairs of files can be specified.\n"
            << "    -p\n"
            << "    -prefix <name>          Prefix for the output files. Default is \"out\".\n"
            << "    -M\n"
            << "    -minOverlap <int>       Minimum size of the overlaps to compute.\n"
            << "                            Default is half of the reads length.\n"
            << "    -t\n"
            << "    -truncate <int>         Truncate the 3' end of the reads TO the specified\n"
            << "                            length\n"
            << "  2) Assembler mode:\n"
            << "    -e\n"
            << "    -edenaFile <file.ovl>   Edena overlap (.ovl) file. Required.\n"
            << "    -p\n"
            << "    -prefix <name>          Prefix for the output files.\n"
            << "    -m\n"
            << "    -overlapCutoff <int>    Only consider overlaps >= than the specified size.\n"
            << "                            The optimal setting of this parameter depends on the\n"
            << "                            coverage that was achieved by the sequencing run.\n"
            << "                            You should therefore try different values in order\n"
            << "                            to get the optimal one.\n"
            << "    -c\n"
            << "    -minContigSize <int>    Minimum size of the contigs to output.\n"
            << "                            Default is 1.5*readLength.\n"
            << "    -cc\n"
            << "    -contextualCleaning\n"
            << "                 <yes/no>   Contextual cleaning of spurious edges.\n"
            << "    -minCoverage <int>      Minimum required coverage for the contigs. This\n"
            << "                            value is automatically determined if not specified.\n"
            << "    -sph\n"
            << "    -shortPeHorizon <int>   Maximum search distance for short >< paired-end\n"
            << "                            sampling. Default: 1000\n"
            << "    -lph\n"
            << "    -longPeHorizon <int>    Maximum search distance for long <> paired-end\n"
            << "                            sampling. Default: 15000\n"
            << "    -peHorizon <int>        obsolete: Maximum search distance for both short\n"
            << "                            and long paired-end reads sampling.\n"
            << "    -trimRed <yes/no>       By default, possible redundant sequences, caused by\n"
            << "                            unresolved repeats, are trimmed from contigs ends.\n"
            << "                            Setting this flag to 'no' allows keeping such \n"
            << "                            redundancies up to the length of the largest insert\n"
            << "                            size (maxJump). !! setting this setting to 'no'\n"
            << "                            produces an artificially increased assembly size !!\n"
            << "    -maxRed <int>           Max ending redundancy length. Default: 0 (equivalent\n"
            << "                            to '-trimRed yes'. Overrides -trimRed.\n"          
            << "    -d\n"
            << "    -deadEnds <int>         Maximum length for dead-end paths removal.\n"
            << "                            Default value is set to 2*readLength-1.\n"
            << "    -discardNonUsable\n"
            << "                 <yes/no>   Reads that are involved in orphan nodes smaller than\n"
            << "                            1.5*readLength are considered as \"non-usable\".\n"
            << "                            This filter discards such nodes. Default: enabled\n"
            << "    -trim <int>             Coverage cutoff for contigs ends:\n"
            << "                            Contig ends ending in a gap may contain errors due\n"
            << "                            to low coverage. This option trim a few bases from\n"
            << "                            these ends until a minimum coverage is reached.\n"
            << "                            Default is 4. Ends are not trimmed if set to 1.\n"
            << "    -shell                  Interactive shell. Allows querying the assembly\n"
            << "                            graph.\n"          
            << "REPORT BUGS:\n"
            << "   david.hernandez@genomic.ch\n";
    exit(exitValue);
};

void sendBugReportPlease(ostream &out) {
    out << "Please, report the problem to\ndavid.hernandez@genomic.ch\n";
    exit(EXIT_FAILURE);
}

void getCompReverse(string source, string &comp) {


    static int done;
    static char REVERSE[256];

    if (done == 0) // first call
    {
        done = 1;

        memset(REVERSE, (int) 'X', 256 * sizeof (char));

        REVERSE[(int) 'A'] = 'T';
        REVERSE[(int) 'T'] = 'A';
        REVERSE[(int) 'C'] = 'G';
        REVERSE[(int) 'G'] = 'C';
        REVERSE[(int) 'R'] = 'Y';
        REVERSE[(int) 'Y'] = 'R';
        REVERSE[(int) 'M'] = 'K';
        REVERSE[(int) 'K'] = 'M';
        REVERSE[(int) 'W'] = 'W';
        REVERSE[(int) 'S'] = 'S';
        REVERSE[(int) 'B'] = 'V';
        REVERSE[(int) 'V'] = 'B';
        REVERSE[(int) 'D'] = 'H';
        REVERSE[(int) 'H'] = 'D';
        REVERSE[(int) 'N'] = 'N';

        REVERSE[(int) 'a'] = 'T';
        REVERSE[(int) 't'] = 'A';
        REVERSE[(int) 'c'] = 'G';
        REVERSE[(int) 'g'] = 'C';
        REVERSE[(int) 'r'] = 'Y';
        REVERSE[(int) 'y'] = 'R';
        REVERSE[(int) 'm'] = 'K';
        REVERSE[(int) 'k'] = 'M';
        REVERSE[(int) 'w'] = 'W';
        REVERSE[(int) 's'] = 'S';
        REVERSE[(int) 'b'] = 'V';
        REVERSE[(int) 'v'] = 'B';
        REVERSE[(int) 'd'] = 'H';
        REVERSE[(int) 'h'] = 'D';
        REVERSE[(int) 'n'] = 'N';

    }

    comp = "";

    for (int i = source.length() - 1; i >= 0; i--) {
        comp += REVERSE[(size_t) source[i]];
    }
}

char compReverse(const char source) {
    static int done;
    static char REVERSE[256];

    if (done == 0) // first call
    {
        done = 1;

        memset(REVERSE, (int) 'X', 256 * sizeof (char));

        REVERSE[(int) 'A'] = 'T';
        REVERSE[(int) 'T'] = 'A';
        REVERSE[(int) 'C'] = 'G';
        REVERSE[(int) 'G'] = 'C';
        REVERSE[(int) 'R'] = 'Y';
        REVERSE[(int) 'Y'] = 'R';
        REVERSE[(int) 'M'] = 'K';
        REVERSE[(int) 'K'] = 'M';
        REVERSE[(int) 'W'] = 'W';
        REVERSE[(int) 'S'] = 'S';
        REVERSE[(int) 'B'] = 'V';
        REVERSE[(int) 'V'] = 'B';
        REVERSE[(int) 'D'] = 'H';
        REVERSE[(int) 'H'] = 'D';
        REVERSE[(int) 'N'] = 'N';

        REVERSE[(int) 'a'] = 'T';
        REVERSE[(int) 't'] = 'A';
        REVERSE[(int) 'c'] = 'G';
        REVERSE[(int) 'g'] = 'C';
        REVERSE[(int) 'r'] = 'Y';
        REVERSE[(int) 'y'] = 'R';
        REVERSE[(int) 'm'] = 'K';
        REVERSE[(int) 'k'] = 'M';
        REVERSE[(int) 'w'] = 'W';
        REVERSE[(int) 's'] = 'S';
        REVERSE[(int) 'b'] = 'V';
        REVERSE[(int) 'v'] = 'B';
        REVERSE[(int) 'd'] = 'H';
        REVERSE[(int) 'h'] = 'D';
        REVERSE[(int) 'n'] = 'N';
    }

    return REVERSE[(size_t)source];
}

void compReverse(const char* source, char* dest, size_t length) {
    for (size_t i = 0; i < length; i++)
        *(dest + i) = compReverse(*(source + length - i - 1));
}

void compReverse(string source, string &comp) {
    comp.clear();
    for (int i = source.length() - 1; i >= 0; i--)
        comp += compReverse(source.at(i));
}

void lineWrap(ostream &out, const string &s, int v) {
    unsigned int p = 0;

    while (p < s.length()) {
        out << s.substr(p, v);
        out << '\n';
        p += v;
    }
}

double fastLog(unsigned int v) {
    static int done;
    static double fLog[MAXLOG + 1];

    if (done == 0) {
        done = 1;
        for (size_t i = 0; i <= MAXLOG; i++)
            fLog[i] = log(i);
        //fLog[0]=0;
    }

    if (v <= MAXLOG)
        return fLog[v];
    else
        //        if (v!=0)
        //                return log(v);
        //        else
        return 0.0;
}

unsigned int hamming(const char* a, const char* b) {
    unsigned int h = 0;
    while (*a != '\0' && *b != '\0') {
        if (*a != *b)
            h++;
        a++;
        b++;
    }
    return h;
}

bool isNotAmbiguous(char*s, unsigned int length) {
    static int done;
    static char TOUPPER[256];

    if (done == 0) // first call
    {
        done = 1;
        memset(TOUPPER, (int) 'X', 256 * sizeof (char));

        TOUPPER[(int) 'A'] = 'A';
        TOUPPER[(int) 'T'] = 'T';
        TOUPPER[(int) 'C'] = 'C';
        TOUPPER[(int) 'G'] = 'G';
        TOUPPER[(int) 'a'] = 'A';
        TOUPPER[(int) 't'] = 'T';
        TOUPPER[(int) 'c'] = 'C';
        TOUPPER[(int) 'g'] = 'G';
    }

    for (size_t i = 0; i < length; i++) {
        s[i] = TOUPPER[(size_t)s[i]];

        if (s[i] == 'X')
            return false;
    }
    return true;
}

//estimate an upper bound value for the number of reads in the file

unsigned int estimateNReads(string fileName, int &rl) {
  
    bool fastQ;

    FILE* in;
    if ((in = fopen(fileName.c_str(), "r")) == NULL) {
        cerr << "[err] Could not open file " << fileName << endl;
        LOG.oss << "[err] Could not open file " << fileName << endl;
        LOG.flushStream();
        return 0;
    }
    
    char c=getc(in);
    fseeko ( in, 0, SEEK_SET );
    if (c== '>')
        fastQ=false;
    else if (c == '@')
        fastQ=true;
    else {
        cout << "Reads file unrecognized format\n";
        cout << "Files must be either FASTA or FASTQ format\n";
        return 0;
    }

    //readLength
    string line;
    
    getline_c(in,line);
    getline_c(in,line);
    
    rl = line.size();
    if (rl > 0) {
        if (line.at(rl - 1) == '\r')
            rl--;
    }

    off_t sp1;

    size_t fileSize;
    size_t chunkSize = 1e6;
    size_t nOctets;
    size_t nChunk = 10;
    unsigned int nLines, maxLines = 0;

    fseeko(in,0,SEEK_END);
    sp1=ftello(in);
    
    if ((int) sp1 == -1) {
        cerr << "File: " << fileName << " cannot be handled\n";
        cerr << "This may arise on 32 bits OS.\n";
        cerr << "if so, you must switch to a 64 bits OS\n";
        return 0;
    } else
        fileSize = sp1;

    if (nChunk * chunkSize > fileSize) {
        nChunk = 1;
        chunkSize = fileSize;
    }

    char *buffer = new char[chunkSize];
    size_t last_p = 0;
    int checkRl;

    char idChar = '>';
    if (fastQ)
        idChar = '@';

    for (size_t i = 0; i < nChunk; i++) {
        sp1 = i * (fileSize / nChunk);
        fseeko(in,0,SEEK_SET);
        nOctets=fread(buffer,1,chunkSize,in);
        nLines = 0;
        last_p = 0;
        bool checkNextLine = false;

        for (size_t p = 0; p < nOctets; p++) {
            if (*(buffer + p) == '\n') {
                nLines++;

                if (checkNextLine && *(buffer + last_p) != idChar) {//last '@' could have occurred in quality line
                    checkRl = p - last_p;

                    if (*(buffer + (p - 1)) == '\r')
                        checkRl--;
                    if (checkRl != rl) {
                        cerr << checkRl << " " << rl << endl;
                        cerr << "[err] All reads within a file must be the same length\n";
                        LOG.write("[err] All reads within a file must be the same length");
                        fclose(in);
                        delete [] buffer;
                        return 0;
                    }
                }

                if (*(buffer + last_p) == idChar)
                    checkNextLine = true;
                else
                    checkNextLine = false;

                last_p = p + 1;
            }
        }
        if (nLines > maxLines)
            maxLines = nLines;
    }

    maxLines /= 2;
    if (fastQ)
        maxLines /= 2;

    fclose(in);
    delete[]buffer;

    float r = (float) fileSize / chunkSize;

    nLines = (unsigned int) (maxLines * r);
    nLines = (unsigned int) (nLines + r);

    return nLines;
}

string smartDNALength(double l) {

    static int done;
    static vector <string> units;

    if (done == 0) {
        units.push_back(" bp ");
        units.push_back(" Kbp");
        units.push_back(" Mbp");
        units.push_back(" Gbp");
        units.push_back(" Tbp");
        done = 1;
    }

    ostringstream oss;
    float d = 1000;
    size_t index = 0;

    while (l >= 1000) {
        l /= d;
        index++;
        if (index == 4)
            break;
    }

    if (index == 0)
        oss << setprecision(0) << fixed << l << units.at(index);
    else
        oss << setprecision(2) << fixed << l << units.at(index);

    return oss.str();
}

void getline_c(FILE *p, string &line)
{
    char c;
    line="";
    c=getc(p); 
    while (!feof(p) && !ferror(p) && c != '\n')
    {
        line+=c;  
        c=getc(p);
    }
}

