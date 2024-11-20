/*
This code runs Monte Carlo simulations for Lennard Jones particles of a single type, in the canonical (NVT) ensemble.

The simulation cell is assumed to be cubic.

R. K. Lindsey (2023)
*/


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "MersenneTwiser.h"
#include "Toolbox.hpp"
#include "system_coordintaes.hpp"
#include "LJ_model.hpp"
#include "inputs.hpp"

using namespace std;

void write_initial(int natoms, string atmtyp, double molmass, double density, double sigma, double epsilon, double rcut, double temp, double nsteps, int iofrq, double redden);

double calculate_initial_energy(system_coordinates &system, LJ_model &LJ, double &stensor_x, double &stensor_y, double &stensor_z);
void run_simulation(system_coordinates &system, LJ_model &LJ, MTRand &mtrand, int nsteps, int nequil, int iofrq, double redden, double temp, double &energy, double &stensor_x, double &stensor_y, double &stensor_z, double &widom_avg, double &widom_trials, int &nrunningav_moves, double natoms, double LJ_epsilon);
void print_statistics(int i, int naccepted_moves, double fraction_accepted, double max_displacement, double energy, double natoms, double LJ_epsilon, double pressure, double temp, double redden, double LJ_sigma, double widom_avg, double widom_trials, double Cv, system_coordinates &system);


int main(int argc, char* argv[]) {

    seed    = stoi(argv[1]);// Seed for random number generator - read from commandline 
    redden  = stod(argv[2]);// Reduced density - read from commandline
    
    // Static variables
    natoms  = 500;          // Number of atoms in the system
    atmtyp  = "Ar";         // Atom type (atomic symbol)
    molmass = 39.948;       // Mass of an particle of atmtyp (Atomic mass in the case of an argon atom)
    density = redden / pow(A2m,3.0) * pow(cm2m,3.0) / nav * molmass / pow(sigma,3.0);     // System density; Units: g/cm^3

    MTRand mtrand(seed); // Seed the (psuedo) random number generator

    system_coordinates system(natoms, atmtyp, density, molmass);
    system.generate_coords();

    write_initial(natoms, atmtyp, molmass, density, sigma, epsilon, rcut, temp, nsteps, iofrq, redden);

    LJ_model LJ(sigma, epsilon, rcut);

    double stensor_x = 0, stensor_y = 0, stensor_z = 0;
    double energy = calculate_initial_energy(system, LJ, stensor_x, stensor_y, stensor_z);

    double widom_avg = 0, widom_trials = 0;
    int nrunningav_moves = 0;

    run_simulation(system, LJ, mtrand, nsteps, nequil, iofrq, redden, temp, energy, stensor_x, stensor_y, stensor_z, widom_avg, widom_trials, nrunningav_moves, natoms, LJ.epsilon);

    return 0;
}

void write_initial(int natoms, string atmtyp, double molmass, double density, double sigma, double epsilon, double rcut, double temp, double nsteps, int iofrq, double redden) {
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
    cout << "# reduc. density (rho*): " << redden     << endl;
    cout << "# reduc. temp (T*):      " << temp / epsilon << endl;
}

double calculate_initial_energy(system_coordinates &system, LJ_model &LJ, double &stensor_x, double &stensor_y, double &stensor_z) {
    double energy = 0;
    xyz rij_vec;
    for (int i=0; i<system.natoms; i++) {
        for (int j=i+1; j<system.natoms; j++) {
            double rij = get_dist(system.coords[i], system.coords[j], system.boxdim, rij_vec);
            if (rij < rcut) {
                energy += LJ.get_eij(rij);
                xyz force;
                LJ.get_fij(rij, rij_vec, force);
                stensor_x += force.x * rij_vec.x;
                stensor_y += force.y * rij_vec.y;
                stensor_z += force.z * rij_vec.z;
            }
        }
    }
    return energy;
}

void run_simulation(system_coordinates &system, LJ_model &LJ, MTRand &mtrand, int nsteps, int nequil, int iofrq, double redden, double temp, double &energy, double &stensor_x, double &stensor_y, double &stensor_z, double &widom_avg, double &widom_trials, int &nrunningav_moves, double natoms, double LJ_epsilon) {
    double max_displacement = 0.5 * system.boxdim.x;
    int naccepted_moves = 0;
    double fraction_accepted;

    for (int i=0; i<nsteps; i++) {
        int selected_atom = int(mtrand()*system.natoms);
        double eold_selected, enew_selected, delta_energy;
        xyz sold_selected, snew_selected, trial_displacement, trial_position;

        LJ.get_single_particle_contributions(system.coords, selected_atom, system.coords[selected_atom], system.boxdim, eold_selected, sold_selected);
        trial_displacement.x = (mtrand()*2.0-1.0)*max_displacement;
        trial_displacement.y = (mtrand()*2.0-1.0)*max_displacement;
        trial_displacement.z = (mtrand()*2.0-1.0)*max_displacement;
        trial_position.x = trial_displacement.x + system.coords[selected_atom].x;
        trial_position.y = trial_displacement.y + system.coords[selected_atom].y;
        trial_position.z = trial_displacement.z + system.coords[selected_atom].z;
        trial_position.x -= floor(trial_position.x/system.boxdim.x)*system.boxdim.x;
        trial_position.y -= floor(trial_position.y/system.boxdim.y)*system.boxdim.y;
        trial_position.z -= floor(trial_position.z/system.boxdim.z)*system.boxdim.z;
        LJ.get_single_particle_contributions(system.coords, selected_atom, trial_position, system.boxdim, enew_selected, snew_selected);

        if (i >= nequil) {
            double ewidom;
            xyz swidom, widom_position;
            widom_position.x = mtrand()*system.boxdim.x;
            widom_position.y = mtrand()*system.boxdim.y;
            widom_position.z = mtrand()*system.boxdim.z;
            LJ.get_single_particle_contributions(system.coords, -1, widom_position, system.boxdim, ewidom, swidom);
            widom_avg   += exp((-1.0 / temp) * ewidom);
            widom_trials++;
        } else {
            widom_avg = 1;
            widom_trials = 1;
        }
        
        delta_energy = enew_selected - eold_selected;
        if (mtrand() < exp(-(1.0/temp)*delta_energy)) {
            system.coords[selected_atom] = trial_position;
            energy += delta_energy;
            stensor_x = stensor_x + snew_selected.x - sold_selected.x;
            stensor_y = stensor_y + snew_selected.y - sold_selected.y;
            stensor_z = stensor_z + snew_selected.z - sold_selected.z;
            naccepted_moves++;
        }

        fraction_accepted = float(naccepted_moves)/float(i+1);
        max_displacement = update_max_displacement(fraction_accepted, system.boxdim.x, max_displacement);

        double numden = redden / pow(A2m,3.0) * pow(cm2m,3.0) / nav * molmass / pow(sigma,3.0);
        double pressure = numden*temp + 1.0/3.0/pow(system.boxdim.x,3.0)*(stensor.x + stensor.y + stensor.z);


        double stat_avgE   = 0;
        double stat_avgEsq = 0;        
        double stat_avgP   = 0;
        double Cv = 0;
        if (i >= nequil) {
            stat_avgE   += energy;
            stat_avgEsq += energy * energy;        
            stat_avgP   += pressure / LJ.epsilon * pow(LJ.sigma, 3.0);
            nrunningav_moves++;
            double avgE   = stat_avgE / float(nrunningav_moves);
            double avgEsq = stat_avgEsq / float(nrunningav_moves);
            Cv = (avgEsq - pow(avgE,2)) / (pow(temp,2)*kB);
        }

        if ((i+1) % iofrq == 0) {
            print_statistics(i, naccepted_moves, fraction_accepted, max_displacement, energy, natoms, LJ.epsilon, pressure, temp, redden, LJ.sigma, widom_avg, widom_trials, Cv, system);
            system.write_frame(i);
        }
    }
}

void print_statistics(int i, int naccepted_moves, double fraction_accepted, double max_displacement, double energy, double natoms, double LJ_epsilon, double pressure, double temp, double redden, double LJ_sigma, double widom_avg, double widom_trials, double Cv, system_coordinates &system) {
    
    cout << "Debug: " << endl;
    cout << LJ_sigma << endl;
    cout << pressure << endl;
    cout << temp << endl; 
    cout << redden << endl; 

    cout << "Step:  " << setw(10) << left << i;
    cout << " NAcc:  " << setw(10) << left << setprecision(3) <<  naccepted_moves;
    cout << " fAcc:  " << setw(10) << left << fixed << setprecision(3) << fraction_accepted;
    cout << " Maxd:  " << setw(10) << left << fixed << setprecision(5) << max_displacement;
    cout << " E*/N:  " << setw(10) << left << fixed << setprecision(5) << energy / natoms / LJ_epsilon;
    cout << " P*:     " << setw(10) << left << fixed << setprecision(5) << (pressure + temp * redden / pow(LJ_sigma, 3.0)) / LJ_epsilon * pow(LJ_sigma, 3.0);
    cout << " P*cold: " << setw(10) << left << fixed << setprecision(5) << pressure / LJ_epsilon * pow(LJ_sigma, 3.0);
    cout << " Mu*_xs: " << setw(10) << left << fixed << setprecision(5) << (-temp) * log(widom_avg / widom_trials);
    cout << " Cv*/N_xs:  " << setw(15) << left << fixed << setprecision(5) << Cv / system.natoms * kB;
    cout << " E(kJ/mol): " << setw(10) << left << fixed << setprecision(3) << energy * 0.008314;
    cout << " P(bar):    " << setw(10) << left << fixed << setprecision(3) << pressure * 0.008314 * 10.0e30 * 1000 / (6.02 * 10.0e23) * 1.0e-5;
    cout << endl;
}
