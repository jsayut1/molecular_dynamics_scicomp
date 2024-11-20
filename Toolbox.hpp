#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>

#include "MersenneTwiser.h"


using namespace std;

#ifndef TOOLBOX_HPP
#define TOOLBOX_HPP

const double  nav     = 6.02*pow(10.0,23);        // Avagadro's number
const double  cm2m    = 1.0/100.0;                // x cm * cm2m = m
const double  A2m     = 1.0/pow(10.0,10);         // x Angstrom * A to m = m    
const double  kB      = 1.380649*pow(10,-23.0);   // Units: J⋅K^−1 

struct xyz
{
    double x;
    double y;
    double z;
};

double get_dist(const xyz & a1, const xyz & a2, const xyz & boxdim, xyz & rij_vec)
{
    /* Write code to calculate the distance vector and scalar distance between two atoms, using the minimum image convention.
    The distance vector is rij_vec - it is passed by reference so it can be modified directly.
    The convention is that the vector should point from atom 1 to atom 2
    
    The function should return the scalar distance.
    */
    rij_vec.x = a2.x - a1.x;
    rij_vec.x -= round(rij_vec.x/boxdim.x)*boxdim.x;
    rij_vec.y = a2.y - a1.y;
    rij_vec.y -= round(rij_vec.y/boxdim.y)*boxdim.y;
    rij_vec.z = a2.z - a1.z;
    rij_vec.z -= round(rij_vec.z/boxdim.z)*boxdim.z;
    return sqrt((rij_vec.x) * (rij_vec.x) + (rij_vec.y) * (rij_vec.y) + (rij_vec.z) * (rij_vec.z));
}

double update_max_displacement(double fraction_accepted, double boxdim, double max_displacement)
{
    /* Write code to update the maximum displacement based on current acceptance crtieria. 
    The function should return the the new maximum displacement.
    */
    double scale_factor = 0.0;
    if(fraction_accepted<0.5)
    {
        scale_factor = 0.9;
    }
    else
    {
        scale_factor = 1.1;
    }
    return std::max(0.01*boxdim, std::min(0.5*boxdim, max_displacement*scale_factor));
}
#endif