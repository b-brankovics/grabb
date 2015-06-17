/* 
 * File:   ActualDistribution.h
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

#ifndef ACTUALDISTRIBUTION_H
#define	ACTUALDISTRIBUTION_H

#include <vector>

using namespace std;

class ActualDistribution
{
public:
    ActualDistribution();
    ActualDistribution(const ActualDistribution& orig);
    virtual ~ActualDistribution();
    
    void setDistribution(vector<unsigned int> &d, int min, int max);
    double cdf(int a); //p(x<=a)
    double count_cdf(int a); // number of element x such that x<=a 
    
    double cdf(int a, int b); //p(a > x >= b)
    double count_cdf(int a, int b); //number of element x such that a > x >= b
    //used to compute the expected number of PE hits in a given node
    double sumCdf(int a, int b, int shift);
    
    double count_sumiPi(int a); // sum(x) for x<=a
    double mu();
    double mu(int a, int b); //mu for sample \in ]a;b]
    //used to compute the expected mean distance of PE hits in a given node
    double mu(int a, int b, int shift);
//    double sigma();
//    double sigma(int a, int b); //sigma for sample \in ]a;b]
   
    
private:


    vector<unsigned int> counts;
    vector<double> CDF; //cumulative distribution function * nSample -> CDF[maxV-1] = nSampe
    vector<double> sumiPi;
    int minV;
    int maxV;
    unsigned int nSample;

    

};

#endif	/* ACTUALDISTRIBUTION_H */

