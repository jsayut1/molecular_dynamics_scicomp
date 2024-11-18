#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>

#include "MersenneTwiser.h"
#include "Toolbox.hpp"


using namespace std;

class LJ_model
{
public:
    
    ///units for epsilon: e/Kb?
    
    double  sigma;
    double  epsilon;
    double  rcut;
    
    double  term1;   // Placeholder for (sigma/rij)^12
    double  term2;   // Placeholder for (sigma/rij)^6
    
    double  get_eij(double rij);
    void    get_fij(double rij, const xyz & rij_vec, xyz & fij);
    void    get_single_particle_contributions(const vector<xyz> & coords, int selected_atom, const xyz & selected_atom_coords, const xyz & boxdim, double & energy_selected, xyz & stress_selected);
    
    LJ_model(double in_sigma, double in_epsilon, double in_rcut);
    ~LJ_model();
};

LJ_model::LJ_model(double in_sigma, double in_epsilon, double in_rcut)
{
    sigma   = in_sigma;
    epsilon = in_epsilon;
    rcut    = in_rcut;
}
LJ_model::~LJ_model(){}

double LJ_model::get_eij(double rij)
{
    /* Write code to compute the Lennard Jones energy between the two atoms.
    It should account for the user-specified cutoff distance.
    The function should return the computed energy.
    */
    return (4*epsilon*((pow((sigma/rij),12)) - (pow((sigma/rij),6))));
}

void LJ_model::get_fij(double rij, const xyz & rij_vec, xyz & fij)
{
    /* Write code to compute the Lennard Jones force between the two atoms.
    It should account for the user-specified cutoff distance.
    The function should update the the force and distance vectors directly
    The function should not return anything.
    */
    fij.x = -(4*epsilon*((-12*pow((sigma),12)*pow((rij),-13)) - (-6*pow((sigma),6)*pow((rij),-7))))*rij_vec.x/abs(rij);
    fij.y = -(4*epsilon*((-12*pow((sigma),12)*pow((rij),-13)) - (-6*pow((sigma),6)*pow((rij),-7))))*rij_vec.y/abs(rij);
    fij.z = -(4*epsilon*((-12*pow((sigma),12)*pow((rij),-13)) - (-6*pow((sigma),6)*pow((rij),-7))))*rij_vec.z/abs(rij);
    return;
}

void LJ_model::get_single_particle_contributions(const vector<xyz> & coords, int selected_atom, const xyz & selected_atom_coords, const xyz & boxdim, double & energy_selected, xyz & stress_selected)
{
    energy_selected   = 0;
    stress_selected.x = 0;
    stress_selected.y = 0;
    stress_selected.z = 0;
    xyz force;
    force.x = 0;
    force.y = 0;
    force.z = 0;
    
    static double rij;
    static xyz    rij_vec;

    /* Write code to determine the contributions to the total system energy, forces, and stresses due to the selected atom.
    Self interactions should not be included.

    // Loop over all atoms. Within the loop:
        
        // Get the scalar distance and distance vector between atoms, using MIC


        // Determine pair energy, but only if interaction is within cutoff idstance
    
    
        // Determine the atom pair's contribution to the total system pressure - again, only perform 
        // if within the model cutoff - we'll use this if the move is accepted
    
    */
    for(int j=0; j<coords.size(); j++)
    {
        // if(coords[j].x==selected_atom_coords.x && coords[j].y==selected_atom_coords.y && coords[j].z==selected_atom_coords.z)
        //     continue;
        if(j!=selected_atom)
        {
            rij = get_dist(selected_atom_coords, coords[j],boxdim,rij_vec);

            if(rij<rcut)
            {
                energy_selected += get_eij(rij);
                get_fij(rij,rij_vec,force);
                stress_selected.x += force.x * rij_vec.x;
                stress_selected.y += force.y * rij_vec.y;
                stress_selected.z += force.z * rij_vec.z;
            }
        }
        //still calculating its distance with itself, too.
    }
}