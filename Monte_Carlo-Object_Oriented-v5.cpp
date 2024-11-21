/* 
This code runs Monte Carlo simulations for Lennard Jones particles of a single type, in the canonical (NVT) ensemble.

The simulation cell is assumed to be cubic.

R. K. Lindsey (2023)
*/


#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include <chrono>  // For high_resolution_clock, duration, etc.

#include "MersenneTwiser.h"
#include "Toolbox.hpp"
#include "system_coordintaes.hpp"
#include "LJ_model.hpp"
#include "mpi.h"

using namespace std;




int main(int argc, char* argv[])
{
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Set user-defined variables (Read in from input file at some point)
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    int     seed    = stoi(argv[1]);// Seed for random number generator - read from commandline 
    double  redden  = stod(argv[2]);// Reduced density - read from commandline
    
    MTRand  mtrand(seed); // Seed the (psuedo) random number generator

    //////// MPI Initialization
    int rank, nprocs, my_natoms, mpierr, rank_IC, rank_owner;
    MPI_Init(&argc, &argv);                          
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);           
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rank_IC=0;

    //////// Static variables 
    
    double  sigma   = 3.4;          // Units: Angstrom
    double  epsilon = 120;          // Units: K/k_B
    double  rcut    = 4*sigma;      // Model outer cutoff
    
    int     natoms  = 500;          // Number of atoms in the system
    string  atmtyp  = "Ar";         // Atom type (atomic symbol)
    double  molmass = 39.948;       // Mass of an particle of atmtyp (Atomic mass in the case of an argon atom)
    double  density = redden / pow(A2m,3.0) * pow(cm2m,3.0) / nav * molmass / pow(sigma,3.0);     // System density; Units: g/cm^3
    double  numden;                 // System number density; Units: atoms/Ang^3 - will be calculated later on

    double  temp    = 1.2*epsilon;  // Temperature, in K
    double  nsteps  = 5e6;          // Number of MC steps
    int     iofrq   = 2e3;          // Frequency to output statistics and trajectory
    int     nequil  = 1e6;          // Equilibration period (chemical potential and heat capacity only collected after this many steps)
    
    //////// Print information for user

    cout << "# Number of atoms:       " << natoms     << endl;
    cout << "# Atom type:             " << atmtyp     << endl;
    cout << "# Molar Mass (g/mol):    " << molmass    << endl;
    cout << "# Density (g/cm^3):      " << density    << endl;
    cout << "# LJ sigma (Angstrom):   " << sigma      << endl;
    cout << "# LJ epsilon/kB (K):     " << epsilon    << endl;
    cout << "# LJ cutoff (Angstrom):  " << rcut       << endl;
    cout << "# LJ cutoff (sigma):     " << rcut/sigma << endl;
    cout << "# Temperature (K):       " << temp       << endl;
    cout << "# Number MC steps:       " << nsteps     << endl;
    cout << "# Output frequency:      " << iofrq      << endl;
    
    //////// Set up MPI and print

    
    my_natoms=natoms/nprocs;
    if ( natoms%nprocs >rank) my_natoms=my_natoms+1;

    cout << "---------" << endl;
    cout << "MPI !" << endl;
    cout << "Rank: " << rank     << endl;
    cout << "My natoms " << my_natoms    << endl;
    cout << "---------" << endl;





    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Initialize the system - arguments are natoms, atom type, density (g/cc), molar mass, 
    // and whether to print a simulation trajectory for object
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    system_coordinates system(natoms, atmtyp, density, molmass);
    system.generate_coords();
    
    if(rank==0){
    system_coordinates ICsystem(natoms, atmtyp, density, molmass);
    ICsystem.generate_coords();
    }
    system_coordinates my_system(my_natoms, atmtyp, density, molmass);


    cout << "# reduc. density (rho*): " << redden /*write this*/ << endl;
    cout << "# reduc. temp (T*):      " << temp / epsilon /*write this*/ << endl;

    // Intialize the model - arguments are sigma (Angstroms), epsilon (K), and outer cutoff (Angstroms)
    
    LJ_model LJ(sigma, epsilon, rcut);
    LJ_model my_LJ(sigma, epsilon, rcut);


    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Determine initial system energy and pressure
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    double  rij;                // Scalar distance between atoms
    xyz     rij_vec;            // Distance vector between atoms
    
    double  energy = 0;         // Energy for the configuration
   
    double widom_factor = 0;    // The exponential term in the Widom chemical potential expression
    double widom_trials = 0;    // Number of Widom insertion attempts
    double widom_avg    = 0;    // Average value of the wWidom factor
        
    double Eavg = 0;            // Average energy
    double Esqavg = 0;          // Square of average energy

   
    xyz     force;              // Force on a given atom
    xyz     stensor;            // System stress tensor (diagonal components only, i.e., xx, yy, zz)
    stensor.x = 0;
    stensor.y = 0;
    stensor.z = 0;



    if(rank==0){

    for (int i=0; i<ICsystem.natoms; i++)                                                                         
    {
        for (int j=i+1; j<ICsystem.natoms; j++)
        {
            // Get the scalar distance and distance vector between atoms, using MIC

            /*write this*/
            rij = get_dist(system.coords[i],system.coords[j],system.boxdim,rij_vec);
                      
            // Determine atom pair's contirbution to total system energy - remember to only perform the 
            // calculation if the pair distance is within the model cutoff
            
            /*write this*/
            if(rij<rcut)
            {
                energy+=my_LJ.get_eij(rij);
            }
            
            // Determine the atom pair's contribution to the total system pressure - again, only perform 
            // if within the model cutoff
                        
            /*write this*/
            if(rij<rcut)
            {
                my_LJ.get_fij(rij, rij_vec, force);
                stensor.x += force.x * rij_vec.x;
                stensor.y += force.y * rij_vec.y;
                stensor.z += force.z * rij_vec.z;
            

            }
        }
    }

    }




    MPI_Finalize();
    return rank_IC;



    //cout << energy/system.natoms/LJ.epsilon << endl;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Begin the simulation
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    double max_displacement = 0.5*system.boxdim.x;  // Start trial displacements at one angstrom

    int     selected_atom;
    
    double  eold_selected;  // Pre-trial-move energy of the selected atom
    double  enew_selected;  // Post-trial-move energy of the selected atom
    double  delta_energy;   // Difference between old and new (trial) energy
    xyz     sold_selected;  // Pre-trial-move stress tensor diagonal of the selected atom
    xyz     snew_selected;  // Post-trial-move stress tensor diagonal of the selected atom

    xyz     trial_displacement;    
    xyz     trial_position;
    
    int     naccepted_moves = 0;    // Number of accepted moves
    double  fraction_accepted;      // Fraction of attempted moves that have been accepted
    int     nrunningav_moves = 0;   // Move index for running averages (doesn't start until equilibration period has ended)

    double pressure;
    double Cv;

    double stat_avgE   = 0;
    double stat_avgEsq = 0;
    double stat_avgP   = 0;
    







    for (int i=0; i<nsteps; i++)
    {
        // MPI: rank_IC selects which rank owns the particle that gets pushed (randomly)
        // if (rank==rank_IC) {
        //     rank_owner=int(mtrand()*nprocs)

        // MPI: rank_IC sends the owner rank to everyone else (MPI_Bcast)
        // int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, 
        //    MPI_Comm comm )
        // } 
        // MPI: everyone receives the owner rank from rank_IC (MPI_Receive)

        // MPI: everyone calls broadcast 
        //     mpierr=MPI_Bcast(&rank_owner,1,MPI_INT,rank_IC,MPI_COMM_WORLD);

        // MPI: owner chooses a particle randomly
        // if (rank==rank_owner) {
            // Select a random particle. The syntax below shows how to use the random number generator. This generate a random integer between 0 and natoms-1
            // selected_atom = int(mtrand()*system.natoms);

            // MPI: owner calculates the trial displacement and position

            // 1. Generate the trial **displacement** in x, y, and z - the particle should be able to move in positive
            // and negative directions, i.e., +/- the maximum displacement

            trial_displacement.x = (mtrand()*2.0-1.0)*max_displacement;
            trial_displacement.y = (mtrand()*2.0-1.0)*max_displacement;
            trial_displacement.z = (mtrand()*2.0-1.0)*max_displacement;
            
            // // 2. Generate the trial **position** in x, y, and z based on the displacement
            
            trial_position.x = trial_displacement.x + system.coords[selected_atom].x;
            trial_position.y = trial_displacement.y + system.coords[selected_atom].y;
            trial_position.z = trial_displacement.z + system.coords[selected_atom].z;
            
            // // 3. Apply PBC if the particle has moved outside the box
            
            trial_position.x -= floor(trial_position.x/system.boxdim.x)*system.boxdim.x;
            trial_position.y -= floor(trial_position.y/system.boxdim.y)*system.boxdim.y;
            trial_position.z -= floor(trial_position.z/system.boxdim.z)*system.boxdim.z;

            // MPI: Convert relevant quantities to C / MPI supported data types
            // nope - selectedAtomCoords=&system.coords[selected_atom][0];
            // selectedAtomCoords[0]=system.coords[selected_atom].x;
            // selectedAtomCoords[1]=system.coords[selected_atom].y;
            // selectedAtomCoords[2]=system.coords[selected_atom].z;
            // trialPosition[0]=trial_position.x;
            // trialPosition[1]=trial_position.y;
            // trialPosition[2]=trial_position.z;

        // }

        // declarations for MPI communication variables
        // double selectedAtomCoords[3];
        // double trialPosition[3];
        // xyz selectedAtomCoordsXYZ


        // MPI: owner sends trial displacement and position to everyone else and everyone else receives (MPI_Bcast)
        // int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, 
        //    MPI_Comm comm )        
        // mpierr=MPI_Bcast(selectedAtomCoords,3,MPI_DOUBLE,rank_owner,MPI_COMM_WORLD);
        // mpierr=MPI_Bcast(trialPosition,3,MPI_DOUBLE,rank_owner,MPI_COMM_WORLD);

        // convert back to C data types
        // system.coords[selected_atom].x=selectedAtomCoords[0];
        // system.coords[selected_atom].y=selectedAtomCoords[1];
        // system.coords[selected_atom].z=selectedAtomCoords[2];
        // trial_position.x=trialPosition[0];
        // trial_position.y=trialPosition[1];
        // trial_position.z=trialPosition[2];

        // MPI: all ranks calculate the energy contribution between the owner's selected atom and all of their atoms
        // MPI: owner calls the normal function, everyone else calls the not_owner version that doesn't need the index of the selected atom

        // if (rank==rank_owner) {
            // Determine contributions to the system's energy, force, and stress due to the selected atom
            LJ.get_single_particle_contributions(system.coords, selected_atom, system.coords[selected_atom], system.boxdim, eold_selected, sold_selected);
            
            // 4. Determine the energy contribution of that particle with the system **in it's trial position**
            LJ.get_single_particle_contributions(system.coords, selected_atom, trial_position, system.boxdim, enew_selected, snew_selected);
        // } else {
            // // Determine contributions to the system's energy, force, and stress due to the selected atom
            // LJ.get_single_particle_contributions_not_owner(system.coords, selectedAtomCoordsXYZ, system.boxdim, eold_selected, sold_selected);
            
            // // 4. Determine the energy contribution of that particle with the system **in it's trial position**
            // LJ.get_single_particle_contributions_not_owner(system.coords, trial_position, system.boxdim, enew_selected, snew_selected);
        // }

        // MPI: all ranks other than owner send their energy to the owner (MPI_Reduce(MPI_SUM)) - synchronization point 2
        // MPI: note - for now, we're not sending the stress. that will require converting to and from xyz again
        // int MPI_Reduce(const void *sendbuf, void *recvbuf, int count,
        //        MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
        // mpierr = MPI_Reduce(enew_selected,enew_selected_total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        if (i >= nequil) // Only do Widom tests for the equilibrated portion of the simulation
        {        
            // 5. Do a widom insertion test
            
            double ewidom;
            xyz    swidom;
            xyz    widom_position;
            
            // 5.a Generate another position for a new ghost atom
            
            /*write this*/
            
            widom_position.x = mtrand()*system.boxdim.x;
            widom_position.y = mtrand()*system.boxdim.y;
            widom_position.z = mtrand()*system.boxdim.z;
            
            // 5.b Calculate change in system energy due to insertion of new ghost atom 
            
            LJ.get_single_particle_contributions(system.coords, -1, widom_position, system.boxdim, ewidom, swidom);
            
            // 5.c Update the Widom factor
            
            widom_factor = exp((-1.0 / (temp)) * ewidom);/*write this*/
            widom_avg   += widom_factor;
            widom_trials++;
        }
        else
        {
            // Needed to avoid a divide-by-zero during equilibration phase when these values aren't collected
            widom_avg    = 1;
            widom_trials = 1;
        }
        
        // MPI: owner decides whether to update the position of the selected atom
        // if so, owner updates the atom's coordinates in its own data (no communication needed)
        // if(rank==rank_owner) {

        // 6. Accept or reject the move
        // If E_old is the energy of the original system and E_new is the system energy when the 
        // particle is at it's trial position, E_old - eold_selected + enew_selected = E_new
        // Therefore delta E, which = E_new - E_old is just enew_selected - eold_selected
        
        delta_energy = enew_selected - eold_selected;

        // MPI: owner decides whether to update the position of the selected atom
        // MPI: if so, owner updates the atom's coordinates in its own data (no communication needed)

        if ( mtrand() < exp(-(1.0/temp)*delta_energy)/*write the acceptance criteria*/ ) // Then accept
        {
            // Then the system energy has decreased **or** our random number is less than our probability to accept
            // Update the system position, energy, and stress tensor, and number of accepted moves
            
            /*write this*/

            system.coords[selected_atom] = trial_position;

            energy += delta_energy;

            stensor.x = stensor.x + snew_selected.x - sold_selected.x;
            stensor.y = stensor.y + snew_selected.y - sold_selected.y;
            stensor.z = stensor.z + snew_selected.z - sold_selected.z;

            naccepted_moves++;
        }

    // }

        // owner boadcasts naccepted_moves to everyone, then everyone updates max_displacement
        // int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, 
        //    MPI_Comm comm )        
        // mpierr=MPI_Bcast(naccepted_moves,1,MPI_INT,rank_owner,MPI_COMM_WORLD);

        // Update maximum diplacement to target a 50% acceptance rate        
        fraction_accepted = float(naccepted_moves)/float(i+1);
        max_displacement = update_max_displacement(fraction_accepted, system.boxdim.x, max_displacement);

        // print statistics if ineeded - don't forget to conver the stress tensor to pressure 
        // Compute instantaneous properties
        
        pressure = numden*temp + 1.0/3.0/pow(system.boxdim.x,3.0)*(stensor.x + stensor.y + stensor.z);
        Cv       = 0;

        if (i >= nequil) // Compute values for running averages, only using the equilibrated portion of the trajectory
        {
            stat_avgE   += energy;
            stat_avgEsq += energy*energy;        
            stat_avgP   += pressure/LJ.epsilon*pow(LJ.sigma,3.0); // Convert to reduced units! ********
            nrunningav_moves++;
            
            double avgE   = stat_avgE /float(nrunningav_moves);
            double avgEsq = stat_avgEsq / float(nrunningav_moves);
            
            Cv = (avgEsq - pow(avgE,2)) / (pow(temp,2)*kB);/*write this - this should only be the dE/dT portion*/
        }
        
   
        if ( (i+1) % iofrq == 0)
        {
            // MPI: everyone communicates things for writing to rank_IC
            // system.write_frame(i);

	    cout << "Debug: " << endl;
	    cout << sigma << endl;
	    cout << pressure << endl;
	    cout << temp << endl; 
	    cout << redden << endl;


            cout << "Step:  " << setw(10) << left << i;
            cout << " NAcc:  " << setw(10) << left << setprecision(3) <<  naccepted_moves;
            cout << " fAcc:  " << setw(10) << left << fixed << setprecision(3) << fraction_accepted;
            cout << " Maxd:  " << setw(10) << left << fixed << setprecision(5) << max_displacement;
            cout << " E*/N:  " << setw(10) << left << fixed << setprecision(5) << energy/natoms/LJ.epsilon;
            cout << " P*:     " << setw(10) << left << fixed << setprecision(5) << (pressure + temp * redden / pow(LJ.sigma,3.0))/LJ.epsilon*pow(LJ.sigma,3.0)/*write this - this is the reduced pressure*/;
            cout << " P*cold: " << setw(10) << left << fixed << setprecision(5) << (pressure)/LJ.epsilon*pow(LJ.sigma,3.0) /*write this - this is the reduced Virial component of pressure*/;
            cout << " Mu*_xs: " << setw(10) << left << fixed << setprecision(5) << (-temp) * log(widom_avg / widom_trials)/*write this - this is excess chemical potential*/; 
            cout << " Cv*/N_xs:  " << setw(15) << left << fixed << setprecision(5) << Cv/system.natoms*kB/*write this - this is excess heat capacity per atom*/;
            cout << " E(kJ/mol): " << setw(10) << left << fixed << setprecision(3) << energy * 0.008314; // KJ/mol per K 
            cout << " P(bar):    " << setw(10) << left << fixed << setprecision(3) << pressure * 0.008314 * 10.0e30 * 1000/(6.02*10.0e23)*1.0e-5; // KJ/mol/A^3 to bar
            cout << endl;
        }

    }
    
    stat_avgE    /= float(nrunningav_moves);
    stat_avgEsq  /= float(nrunningav_moves);
    Cv            =  (stat_avgEsq - pow(stat_avgE,2)) / (kB * pow(temp,2));/*write this, based on average energies - this should only be the dE/dT portion*/
    stat_avgE    *= 1.0/natoms/LJ.epsilon;

    // Finalize MPI
    // MPI_finalize();

    cout << "# Computed average properties: " << endl;
    cout << " # E*/N:  "      << setw(10) << left << fixed << setprecision(5) << stat_avgE << endl;
    cout << " # P*:     "     << setw(10) << left << fixed << setprecision(5) << (stat_avgP / float(nrunningav_moves))+(temp * redden / LJ.epsilon) << endl;
    cout << " # Cv*/N_xs:   " << setw(15) << left << fixed << setprecision(5) << Cv/system.natoms*kB/*write this - this is excess heat capacity per atom based on average energies*/ << endl;
    cout << " # Mu_xs:  "     << setw(10) << left << fixed << setprecision(5) << (-temp) * log(widom_avg / widom_trials)/*write this - this is excess chemical potential*/ << endl;       
    
}
















