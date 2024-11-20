#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>

#include "MersenneTwiser.h"
#include "Toolbox.hpp"

using namespace std;

class system_coordinates
{
    
public: 
    // General system coordinate definitions
    
    int         natoms;     // Number of atoms in the system
    string      atmtyp;     // Atom type (atomic symbol)
    double      density;    // System density in g/cm^3
    double      numden;     // System number density (atoms/Ang^3)
    double      molmass;    // Molar mass of an atom of atmtyp
    xyz         boxdim;     // Simulation box x, y, and z lengths
    vector<xyz> coords;     // Coordinates for all atoms in our system
    
    // Variables used for file i/o
    
    ofstream    trajstream; // Trajectory file - uses the LAMMPS lammpstrj format

    void generate_coords();
    void write_frame(int mcstep);
    
    // Constructor and deconstructor
    
    system_coordinates(int in_natoms, string in_atmtyp, double in_density, double in_molmass);
    ~system_coordinates();
};

system_coordinates::system_coordinates(int in_natoms, string in_atmtyp, double in_density, double in_molmass)
{
    natoms  = in_natoms;
    atmtyp  = in_atmtyp;
    density = in_density;
    molmass = in_molmass;
    
    trajstream.open("MC_traj.lammpstrj");
}
system_coordinates::~system_coordinates()
{
    trajstream.close();
}

void system_coordinates::generate_coords()
 ///////////////////////////////////////////////////////////////////////////////////////////////
 // Generate the system coordinates
 // Generates an initial configuation of atom on a 3D lattice  with a user-specified number of 
 // atoms, at the user-specified density. 
 ///////////////////////////////////////////////////////////////////////////////////////////////
{
    // Determine the box length that will yield the user-specified density. Start by computing the target number density.

    numden = density / molmass * nav / pow(cm2m,3.0) * pow(A2m,3.0); // Density in atoms per Ang^3
    

    cout << "# Num. den. (atoms/Ang): " << numden << endl;
    
    boxdim.x = pow(natoms/numden, 1.0/3.0); 
    boxdim.y = boxdim.x;
    boxdim.z = boxdim.x;
    
    cout << "# Box length (x):        " << boxdim.x << endl;
    cout << "# Box length (y):        " << boxdim.y << endl;
    cout << "# Box length (z):        " << boxdim.z << endl;

    // Compute the number of gridpoints to use in each direction (ceiling of the cube root of number of atoms).
    // Set the spacing in each dimension based on the box length and the number of gridpoints+1 (to prevent overlap
    // across the periodic boundary. 
    
    int     ngridpoints = ceil(pow(natoms,1.0/3.0));
    double  spacing     = boxdim.x / (ngridpoints+1);
    
    cout << "# Init. spacing (Ang):   " << spacing << endl;

    xyz     tmp_atom;
    int     added_atoms = 0;
    
    for (int x=0; x<ngridpoints; x++)
    {
        for(int y=0; y<ngridpoints; y++)
        {
            for(int z=0; z<ngridpoints; z++)  
            {
                if(added_atoms >= natoms)
                    continue;
                
                tmp_atom.x = x * spacing;
                tmp_atom.y = y * spacing;
                tmp_atom.z = z * spacing;

                coords.push_back(tmp_atom);
                
                added_atoms++;
            }
        }
    }
                
    
    // Print this initial configuration to a file
    
    write_frame(-1);   
}

void system_coordinates::write_frame(int mcstep)
{
    trajstream << "ITEM: TIMESTEP" << endl;
    trajstream << mcstep << endl;
    trajstream << "ITEM: NUMBER OF ATOMS" << endl;
    trajstream << natoms << endl;
    trajstream << "ITEM: BOX BOUNDS pp pp pp" << endl;
    trajstream << "0 " << boxdim.x << endl;
    trajstream << "0 " << boxdim.y << endl;
    trajstream << "0 " << boxdim.z << endl;
    trajstream << "ITEM: ATOMS id type element xu yu zu" << endl;
    
    for (int i=0; i<natoms; i++)
        trajstream << i << " 1 " << atmtyp << " " << coords[i].x << " "<< coords[i].y << " "  << coords[i].z << endl;
}