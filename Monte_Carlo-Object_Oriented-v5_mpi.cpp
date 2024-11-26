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




MPI_Datatype system_coordinates_type;

void create_system_coordinates_mpi(MPI_Datatype *system_coordinates_type) {
    int block_lengths[5] = {1, 1, 1, 1, 3};  // Exclude string and vector
    MPI_Aint displacements[5];
    MPI_Datatype types[5] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; 

    displacements[0] = offsetof(system_coordinates, natoms);
    displacements[1] = offsetof(system_coordinates, density);
    displacements[2] = offsetof(system_coordinates, numden);
    displacements[3] = offsetof(system_coordinates, molmass);
    displacements[4] = offsetof(system_coordinates, boxdim);

    MPI_Type_create_struct(5, block_lengths, displacements, types, system_coordinates_type);
    MPI_Type_commit(system_coordinates_type);
}

int main(int argc, char* argv[])
{
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Set user-defined variables (Read in from input file at some point)
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    int     seed    = stoi(argv[1]);// Seed for random number generator - read from commandline 
    double  redden  = stod(argv[2]);// Reduced density - read from commandline
    
    MTRand  mtrand(seed); // Seed the (psuedo) random number generator

    //////// MPI Initialization
    int rank, nprocs, my_natoms, mpierr, rank_IC, rank_owner,my_atom_start;
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
    double stensor_init_x,stensor_init_y,stensor_init_z;
    double stensor_x,stensor_y,stensor_z;
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
    my_atom_start=rank*(natoms/nprocs) + min(rank,natoms%nprocs);


    cout << "---------" << endl;
    cout << "MPI !" << endl;
    cout << "Rank: " << rank     << endl;
    cout << "My natoms " << my_natoms    << endl;
    cout << "My starting atom " << my_atom_start    << endl;
    cout << "---------" << endl;





    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Initialize the system - arguments are natoms, atom type, density (g/cc), molar mass, 
    // and whether to print a simulation trajectory for object
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    system_coordinates system(natoms, atmtyp, density, molmass);
    //system.generate_coords();
    

    create_system_coordinates_mpi(&system_coordinates_type);
    system_coordinates ICsystem(natoms, atmtyp, density, molmass);
    

    if(rank==0){
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
            rij = get_dist(ICsystem.coords[i],ICsystem.coords[j],ICsystem.boxdim,rij_vec);
                      
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

    stensor_init_x=stensor.x;
    stensor_init_y=stensor.y;
    stensor_init_z=stensor.z;

    }


    // Broadcast ICsystem
    MPI_Bcast(&ICsystem, 1, system_coordinates_type, 0, MPI_COMM_WORLD);

    // Broadcast coords (std::vector<xyz>) manually
    int coords_size = ICsystem.coords.size();
    MPI_Bcast(&coords_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        ICsystem.coords.resize(coords_size);
    }
    MPI_Bcast(ICsystem.coords.data(), coords_size * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if (rank != 0) {
        std::cout << "Process " << rank << " received data: " << ICsystem.natoms << ", " 
                  << ICsystem.density << ", " << ICsystem.numden << ", " << ICsystem.molmass << ", " 
                  << ICsystem.boxdim.x << ", " << ICsystem.boxdim.y << ", " << ICsystem.boxdim.z << std::endl;
    }
    // TODO: fix the local system coordinates below, should not generate coordinates more than once


    // cout << "---------" << endl;
    // cout << "Local system coords!" << endl;
    // cout << "Rank: " << rank     << endl;
    // //my_system.push_back();
    // cout << "   " << endl;
    // cout << "IC system coords " << endl;
    // cout << ICsystem.coords[1].z << endl;
    // cout << ICsystem.coords[2].z << endl;
    // cout << ICsystem.coords[497].x << endl;
    // cout << ICsystem.coords[498].x << endl;
    // cout << ICsystem.coords[499].z << endl;
    // cout << "   " << endl;

    int icnt;
    xyz blank;
    cout << "Hello" << endl;
    cout << my_system.coords.size() << endl;

    for (icnt=0; icnt<my_natoms; icnt++){
        blank.x=ICsystem.coords[icnt+my_atom_start].x;
        blank.y=ICsystem.coords[icnt+my_atom_start].y;
        blank.z=ICsystem.coords[icnt+my_atom_start].z;
        my_system.coords.push_back(blank);

    }

    // set boxdim of my_system based on ICsystem
    my_system.boxdim.x=ICsystem.boxdim.x;
    my_system.boxdim.y=ICsystem.boxdim.y;
    my_system.boxdim.z=ICsystem.boxdim.z;

    // cout << "Local system coords " << endl;
    // cout << my_system.coords[1].z << endl;
    // cout << my_system.coords[2].z << endl;
    // cout << my_system.coords[247].x << endl;
    // cout << my_system.coords[248].x << endl;
    // cout << my_system.coords[249].z << endl;
    // cout << "   " << endl;

    // broadcast stress tensor values
    mpierr=MPI_Bcast(&stensor_init_x,1,MPI_DOUBLE,rank_IC,MPI_COMM_WORLD);
    mpierr=MPI_Bcast(&stensor_init_y,1,MPI_DOUBLE,rank_IC,MPI_COMM_WORLD);
    mpierr=MPI_Bcast(&stensor_init_z,1,MPI_DOUBLE,rank_IC,MPI_COMM_WORLD);


    // set stensor of the local stress tensor to be equal to what it is for the IC_system
    if(rank!=0){
        stensor.x=stensor_init_x;
        stensor.y=stensor_init_y;
        stensor.z=stensor_init_z;

    }


    //MPI_Finalize();
    //return rank_IC;



    //cout << energy/system.natoms/LJ.epsilon << endl;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Begin the simulation
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    double max_displacement = 0.5*ICsystem.boxdim.x;  // Start trial displacements at one angstrom

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
    int entered_loop;

    double pressure;
    double Cv;

    double stat_avgE   = 0;
    double stat_avgEsq = 0;
    double stat_avgP   = 0;
    




    // declarations for MPI communication variables
    double selectedAtomCoords[3];
    double trialPosition[3];
    xyz selectedAtomCoordsXYZ;
    double enew_selected_total;
    double eold_selected_total;
    double snew_total_x;
    double snew_total_y;
    double snew_total_z;
    double sold_total_x;
    double sold_total_y;
    double sold_total_z;
    double snew_selected_x;
    double snew_selected_y;
    double snew_selected_z;
    double sold_selected_x;
    double sold_selected_y;
    double sold_selected_z;



    for (int i=0; i<nsteps; i++)
    {
        // Select a random particle. The syntax below shows how to use the random number generator. This generate a random integer between 0 and natoms-1
        // MPI: rank_IC selects which rank owns the particle that gets pushed (randomly)
        if (rank==rank_IC) {
            rank_owner=int(mtrand()*nprocs);

        // MPI: rank_IC sends the owner rank to everyone else (MPI_Bcast)
        } 

        // MPI: everyone calls broadcast 
        mpierr=MPI_Bcast(&rank_owner,1,MPI_INT,rank_IC,MPI_COMM_WORLD);

        // MPI: owner chooses a particle randomly

        


        // // MPI: owner chooses a particle randomly
        if (rank==rank_owner) {
            // Select a random particle. The syntax below shows how to use the random number generator. This generate a random integer between 0 and natoms-1
            selected_atom = int(mtrand()*my_natoms);

            // MPI: owner calculates the trial displacement and position
            // cout << "selected atom local index " << endl;
            // cout << selected_atom << endl;
            // cout << my_system.coords[selected_atom].x << endl;
            // cout << my_system.coords[selected_atom].y << endl;
            // cout << my_system.coords[selected_atom].z << endl;

            // 1. Generate the trial **displacement** in x, y, and z - the particle should be able to move in positive
            // and negative directions, i.e., +/- the maximum displacement

            trial_displacement.x = (mtrand()*2.0-1.0)*max_displacement;
            trial_displacement.y = (mtrand()*2.0-1.0)*max_displacement;
            trial_displacement.z = (mtrand()*2.0-1.0)*max_displacement;
            
            // // 2. Generate the trial **position** in x, y, and z based on the displacement
            
            trial_position.x = trial_displacement.x + my_system.coords[selected_atom].x;
            trial_position.y = trial_displacement.y + my_system.coords[selected_atom].y;
            trial_position.z = trial_displacement.z + my_system.coords[selected_atom].z;
            
            // // 3. Apply PBC if the particle has moved outside the box

            trial_position.x -= floor(trial_position.x/my_system.boxdim.x)*my_system.boxdim.x;
            trial_position.y -= floor(trial_position.y/my_system.boxdim.y)*my_system.boxdim.y;
            trial_position.z -= floor(trial_position.z/my_system.boxdim.z)*my_system.boxdim.z;

            // MPI: Convert relevant quantities to C / MPI supported data types
            // nope - selectedAtomCoords=&system.coords[selected_atom][0];
            selectedAtomCoords[0]=my_system.coords[selected_atom].x;
            selectedAtomCoords[1]=my_system.coords[selected_atom].y;
            selectedAtomCoords[2]=my_system.coords[selected_atom].z;
            trialPosition[0]=trial_position.x;
            trialPosition[1]=trial_position.y;
            trialPosition[2]=trial_position.z;

            // cout << "selected atom coords " << endl;
            // cout << selectedAtomCoords[0] << endl;
            // cout << selectedAtomCoords[1] << endl;
            // cout << selectedAtomCoords[2] << endl;
        }


        // MPI: owner sends trial displacement and position to everyone else and everyone else receives (MPI_Bcast)
        // int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, 
        //    MPI_Comm comm )        
        mpierr=MPI_Bcast(selectedAtomCoords,3,MPI_DOUBLE,rank_owner,MPI_COMM_WORLD);
        mpierr=MPI_Bcast(trialPosition,3,MPI_DOUBLE,rank_owner,MPI_COMM_WORLD);

        // convert back to C data types
        selectedAtomCoordsXYZ.x=selectedAtomCoords[0];
        selectedAtomCoordsXYZ.y=selectedAtomCoords[1];
        selectedAtomCoordsXYZ.z=selectedAtomCoords[2];
        trial_position.x=trialPosition[0];
        trial_position.y=trialPosition[1];
        trial_position.z=trialPosition[2];


        // cout << "selected atom coords " << endl;
        // cout << selectedAtomCoordsXYZ.x << endl;
        // cout << selectedAtomCoordsXYZ.y << endl;
        // cout << selectedAtomCoordsXYZ.z << endl;
        // cout << my_system.coords[248].x << endl;
        // cout << my_system.coords[249].z << endl;
        


        // MPI: all ranks calculate the energy contribution between the owner's selected atom and all of their atoms
        // MPI: owner calls the normal function, everyone else calls the not_owner version that doesn't need the index of the selected atom

        if (rank==rank_owner) {
            // Determine contributions to the system's energy, force, and stress due to the selected atom
            LJ.get_single_particle_contributions(my_system.coords, selected_atom, my_system.coords[selected_atom], my_system.boxdim, eold_selected, sold_selected);
            
            // 4. Determine the energy contribution of that particle with the system **in it's trial position**
            LJ.get_single_particle_contributions(my_system.coords, selected_atom, trial_position, my_system.boxdim, enew_selected, snew_selected);
        
            // make stress into doubles
            sold_selected_x=sold_selected.x;
            sold_selected_y=sold_selected.y;
            sold_selected_z=sold_selected.z;
            snew_selected_x=snew_selected.x;
            snew_selected_y=snew_selected.y;
            snew_selected_z=snew_selected.z;

        
        } else {
            // // Determine contributions to the system's energy, force, and stress due to the selected atom
            LJ.get_single_particle_contributions_not_owner(my_system.coords, selectedAtomCoordsXYZ, my_system.boxdim, eold_selected, sold_selected);
            
            // // 4. Determine the energy contribution of that particle with the system **in it's trial position**
            LJ.get_single_particle_contributions_not_owner(my_system.coords, trial_position, my_system.boxdim, enew_selected, snew_selected);

            // make stress into doubles
            sold_selected_x=sold_selected.x;
            sold_selected_y=sold_selected.y;
            sold_selected_z=sold_selected.z;
            snew_selected_x=snew_selected.x;
            snew_selected_y=snew_selected.y;
            snew_selected_z=snew_selected.z;
        }


        // if (rank==rank_owner) {
        //     cout << "enew_selected_total" <<  enew_selected_total << endl;
        // }
        

        // MPI: all ranks other than owner send their energy to the owner (MPI_Reduce(MPI_SUM)) - synchronization point 2
        // MPI: note - for now, we're not sending the stress. that will require converting to and from xyz again
        // int MPI_Reduce(const void *sendbuf, void *recvbuf, int count,
        //        MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
        mpierr = MPI_Reduce(&enew_selected,&enew_selected_total,1,MPI_DOUBLE,MPI_SUM,rank_owner,MPI_COMM_WORLD);
        mpierr = MPI_Reduce(&eold_selected,&eold_selected_total,1,MPI_DOUBLE,MPI_SUM,rank_owner,MPI_COMM_WORLD);

        // MPI reduces for each component of the stress
        mpierr = MPI_Reduce(&sold_selected_x,&sold_total_x,1,MPI_DOUBLE,MPI_SUM,rank_owner,MPI_COMM_WORLD);
        mpierr = MPI_Reduce(&sold_selected_y,&sold_total_y,1,MPI_DOUBLE,MPI_SUM,rank_owner,MPI_COMM_WORLD);
        mpierr = MPI_Reduce(&sold_selected_z,&sold_total_z,1,MPI_DOUBLE,MPI_SUM,rank_owner,MPI_COMM_WORLD);
        mpierr = MPI_Reduce(&snew_selected_x,&snew_total_x,1,MPI_DOUBLE,MPI_SUM,rank_owner,MPI_COMM_WORLD);
        mpierr = MPI_Reduce(&snew_selected_y,&snew_total_y,1,MPI_DOUBLE,MPI_SUM,rank_owner,MPI_COMM_WORLD);
        mpierr = MPI_Reduce(&snew_selected_z,&snew_total_z,1,MPI_DOUBLE,MPI_SUM,rank_owner,MPI_COMM_WORLD);


        // MPI_Finalize();
        // return rank_IC;

        //if (i >= nequil) // Only do Widom tests for the equilibrated portion of the simulation
        if (i>=10000000000)
        {        
            // 5. Do a widom insertion test
            
            double ewidom;
            xyz    swidom;
            xyz    widom_position;
            
            // 5.a Generate another position for a new ghost atom
            
            /*write this*/
            
            widom_position.x = mtrand()*my_system.boxdim.x;
            widom_position.y = mtrand()*my_system.boxdim.y;
            widom_position.z = mtrand()*my_system.boxdim.z;
            
            // 5.b Calculate change in system energy due to insertion of new ghost atom 
            
            LJ.get_single_particle_contributions(my_system.coords, -1, widom_position, my_system.boxdim, ewidom, swidom);
            
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
        if (rank==rank_owner){
        delta_energy = enew_selected - eold_selected;
        
        entered_loop=0;
        // MPI: owner decides whether to update the position of the selected atom
        // MPI: if so, owner updates the atom's coordinates in its own data (no communication needed)

        if ( mtrand() < exp(-(1.0/temp)*delta_energy)/*write the acceptance criteria*/ ) // Then accept
        {
            // Then the system energy has decreased **or** our random number is less than our probability to accept
            // Update the system position, energy, and stress tensor, and number of accepted moves
            
            /*write this*/

            my_system.coords[selected_atom] = trial_position;

            energy += delta_energy;



            stensor.x = stensor.x + snew_total_x - sold_total_x;
            stensor.y = stensor.y + snew_total_y - sold_total_y;
            stensor.z = stensor.z + snew_total_z - sold_total_z;

            naccepted_moves++;
            entered_loop=1;

            stensor_x=stensor.x;
            stensor_y=stensor.y;
            stensor_z=stensor.z;




        }
        
        }


    // }

        // owner boadcasts naccepted_moves to everyone, then everyone updates max_displacement
        // int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, 
        //    MPI_Comm comm )       

        MPI_Bcast(&stensor_x,1,MPI_DOUBLE,rank_owner,MPI_COMM_WORLD);
        MPI_Bcast(&stensor_y,1,MPI_DOUBLE,rank_owner,MPI_COMM_WORLD);
        MPI_Bcast(&stensor_z,1,MPI_DOUBLE,rank_owner,MPI_COMM_WORLD);

        MPI_Bcast(&entered_loop,1,MPI_INT,rank_owner,MPI_COMM_WORLD);

        MPI_Bcast(&naccepted_moves,1,MPI_INT,rank_owner,MPI_COMM_WORLD);

        if(entered_loop==1 & rank!=rank_owner){ // add stensors back into the system so that 
            stensor.x=stensor_x;
            stensor.y=stensor_y;
            stensor.z=stensor_z;




        }


        // Update maximum diplacement to target a 50% acceptance rate        
        fraction_accepted = float(naccepted_moves)/float(i+1);
        max_displacement = update_max_displacement(fraction_accepted, my_system.boxdim.x, max_displacement);




        if (rank==rank_owner){
        // print statistics if ineeded - don't forget to conver the stress tensor to pressure 
        // Compute instantaneous properties
        
        pressure = numden*temp + 1.0/3.0/pow(my_system.boxdim.x,3.0)*(stensor.x + stensor.y + stensor.z);
        Cv       = 0;

        if (i >= nequil) // Compute values for running averages, only using the equilibrated portion of the trajectory
        {
            stat_avgE   += energy;
            stat_avgEsq += energy*energy;        
            stat_avgP   += pressure/my_LJ.epsilon*pow(my_LJ.sigma,3.0); // Convert to reduced units! ********
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
            cout << " E*/N:  " << setw(10) << left << fixed << setprecision(5) << energy/natoms/my_LJ.epsilon;
            cout << " P*:     " << setw(10) << left << fixed << setprecision(5) << (pressure + temp * redden / pow(my_LJ.sigma,3.0))/my_LJ.epsilon*pow(my_LJ.sigma,3.0)/*write this - this is the reduced pressure*/;
            cout << " P*cold: " << setw(10) << left << fixed << setprecision(5) << (pressure)/my_LJ.epsilon*pow(my_LJ.sigma,3.0) /*write this - this is the reduced Virial component of pressure*/;
            // cout << " Mu*_xs: " << setw(10) << left << fixed << setprecision(5) << (-temp) * log(widom_avg / widom_trials)/*write this - this is excess chemical potential*/; 
            cout << " Cv*/N_xs:  " << setw(15) << left << fixed << setprecision(5) << Cv/my_system.natoms*kB/*write this - this is excess heat capacity per atom*/;
            cout << " E(kJ/mol): " << setw(10) << left << fixed << setprecision(3) << energy * 0.008314; // KJ/mol per K 
            cout << " P(bar):    " << setw(10) << left << fixed << setprecision(3) << pressure * 0.008314 * 10.0e30 * 1000/(6.02*10.0e23)*1.0e-5; // KJ/mol/A^3 to bar
            cout << endl;
        }
        
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
















