/* 
 * File:   customString.h
 * Author: david
 *
 * Created on June 10, 2012, 2:54 PM
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

#ifndef CUSTOMSTRING_H
#define	CUSTOMSTRING_H

inline int r_strcmp(const char* a, const char* b)
{

    register const unsigned char *s1 = (const unsigned char *) a;
    register const unsigned char *s2 = (const unsigned char *) b;
    unsigned char c1, c2;

    do
    {
        c1 = (unsigned char) *s1--;
        c2 = (unsigned char) *s2--;
        if (c1 == '\0')
            return c1 - c2;
    }
    while (c1 == c2);

    return c1 - c2;

}

//test
inline int my_strcmp(const char* a, const char* b)
{

    register const unsigned char *s1 = (const unsigned char *) a;
    register const unsigned char *s2 = (const unsigned char *) b;
    unsigned char c1, c2;

    do
    {
        c1 = (unsigned char) *s1++;
        c2 = (unsigned char) *s2++;
        if (c1 == '\0')
            return c1 - c2;
    }
    while (c1 == c2);

    return c1 - c2;
}

inline bool d_prefixOf(const char *p1, const char *p2)
{
    while (*p1 != '\0')
    {
        if (*p1 != *p2)
            return false;
        p1++;
        p2++;
    }
    return true;
}

inline bool r_prefixOf(const char *p1, const char *p2)
{
    while (*p1 != '\0')
    {
        if (*p1 != *p2)
            return false;
        p1--;
        p2--;
    }
    return true;
}

#endif	/* CUSTOMSTRING_H */

