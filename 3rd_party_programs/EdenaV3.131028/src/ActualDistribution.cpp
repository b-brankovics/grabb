/* 
 * File:   ActualDistribution.cpp
 * Author: david
 * 
 * Created on June 28, 2012, 10:37 AM
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


#include "ActualDistribution.h"
#include "globalFunc.h"
#include <iostream>
using namespace std;

ActualDistribution::ActualDistribution()
{
}

ActualDistribution::ActualDistribution(const ActualDistribution& orig)
{
}

ActualDistribution::~ActualDistribution()
{
}

void ActualDistribution::setDistribution(vector<unsigned int> &d, int min, int max)
{
    minV = min;
    maxV = max;
    unsigned int range=maxV-minV+1;

    if (maxV - minV + 1 != (int)d.size())
    {
        cout << "ActualDistribution::setDistribution(...) problem\n";
        sendBugReportPlease(cout);
    }
    
    counts.clear();
    CDF.clear();
    sumiPi.clear();
    CDF.assign(range, 0.0);
    sumiPi.assign(range, 0.0);
    nSample=0;
    
    counts.assign(d.begin(), d.end());

    CDF.assign(range, 0.0);
    sumiPi.assign(range, 0.0);

    if (range > 0)
    {
        CDF[0] = counts[0];
        sumiPi[0] = minV * counts[0];

        for (unsigned int i = 1; i < range; i++)
        {
            CDF[i] = counts[i] + CDF[i - 1];
            sumiPi[i] = (minV+i) * counts[i] + sumiPi[i - 1];
        }
        
        nSample=CDF[range-1];
    }
}

double ActualDistribution::cdf(int a)
{
    if (a < minV)
        return 0.0;
    if (a > maxV)
        return 1.0;
    
    return (double) CDF[a-minV]/nSample;

}

double ActualDistribution::count_cdf(int a)
{
    if (a < minV)
        return 0.0;
    if (a > maxV)
        return nSample;

    return (double) CDF[a - minV];
}

double ActualDistribution::cdf(int a, int b)
{
    if (a<=b)
        return cdf(b)-cdf(a);
    return 0.0;
}

double ActualDistribution::count_cdf(int a, int b)
{
    if (a<=b)
        return count_cdf(b)-count_cdf(a);
    return 0.0;
}

double ActualDistribution::sumCdf(int a, int b, int shift)
{
    double sump=0.0;
    
    for (int i=0; i<shift; i++)
    {
        if (a+i > maxV)
            break;
        sump+=count_cdf(a+i, b+i);
    }
        
    return sump/nSample;
}

double ActualDistribution::count_sumiPi(int a) // sum(x) for x<=a
{
    if (a < minV)
        return 0.0;
    if (a > maxV)
        return sumiPi[maxV-minV];

    return (double) sumiPi[a - minV];
}

double ActualDistribution::mu()
{
    return (double)sumiPi[maxV-minV]/nSample;
}

double ActualDistribution::mu(int a, int b)
{
    if (b < minV)
        return 0.0;
    if (a > maxV)
        return 0.0;
    if (b<a)
        return 0.0;
    
    return (count_sumiPi(b) - count_sumiPi(a)) / count_cdf(a,b);
   //return (double)(sumiPi[a-minV]-sumiPi[b-minV]) / (CDF[a-minV]-CDF[b-minV]);
}

double ActualDistribution::mu(int a, int b, int shift)
{
    double sum=0.0,nS=0.0;
    
    for (int i=0; i<shift; i++)
    {
        if (b+i < minV)
            continue;
        if (a+i > maxV)
            break;
        sum += count_sumiPi(b+i) - count_sumiPi(a+i);
        nS += count_cdf(a+i,b+i);
    }
    
    return sum/nS;
}
 

