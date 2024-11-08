/* 
This code computes the radial pair distribution function for a simulation trajectory.

The simulation cell is assumed to be cubic and assumes a .lammpstrj trajectory format with the following fields:
id, type, element, xu, yu, zu

Assumes the configuration has only a single atom type

R. K. Lindsey (2023)
*/


#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<cmath>


using namespace std;

static const double pi      = 3.14159265359;

// Helper data structures and functions

struct xyz
{
    double x;
    double y;
    double z;
};

bool get_next_line(ifstream & instream, string & line)
// Read a line and return it by reference. Returns true if successful, false if error encountered (should detect eof), with error checking.
{
    getline(instream, line);
    if (! instream.good())
        return false;
    else
        return true;
}

string get_next_line(ifstream & instream)
// Read a line and return it by reference. Returns true if successful, false if error encountered (should detect eof), with error checking.
{
    string line;
    
    getline(instream, line) ;
    if (! instream.good())
    {
        cout << "Encountered unexpected error while reading file." << endl;
        exit(1);
        return line;
    }
    else
        return line;
}

int tokenize(string line, vector<string>& tokens)
// Break a line up into tokens based on space separators. Returns the number of tokens parsed.
{
    string buf;
  
    stringstream stream_parser;

    int pos = line.find('\n');   // Strip off terminal new line.
    if ( pos != string::npos ) 
        line.erase(pos, 1);

    stream_parser.str(line);
	 
    tokens.clear();

    while ( stream_parser >> buf ) 
        tokens.push_back(buf);

    return(tokens.size() );
}

string get_token(string line, int field)
// Returns a specific token from a line
{
    vector<string> tokens;
    tokenize(line, tokens);
    
    return tokens[field];
}


// Class for reading and processing system coordinates

class configuration
{
    public:
        
         // General configuration coordinate definitions
        
        int         natoms;     // Number of atoms in the configuration
        xyz         boxdims;    // Box dimensions
        double      density;    // In atoms/Ang^3
        vector<xyz> coords;     // Coordinates for all atoms in our configuration
        
        // Variables used for file i/o
    
        ifstream    trajstream; // Trajectory file - reads the LAMMPS lammpstrj format

        bool    read_frame();      
        double  get_dist(int i, int j);  
        
        // Constructor and deconstructor
        
        configuration(string traj_file);
        ~configuration();
        
};

configuration::configuration(string traj_file)
{
    cout << "Will read from trajectory file: " << traj_file << endl;
    trajstream.open(traj_file);
}

configuration::~configuration()
{
    trajstream.close();
}

bool configuration::read_frame()
{    
    string  line;
    static bool called_before = false;
    
    // Read/ignore unused headers, ensure we haven't reached the end of the trajectory
    
    if (!get_next_line(trajstream, line)) // ITEM: TIMESTEP 
    {
        cout << "End of trajectory file detected." << endl;
            return false;
    }
    get_next_line(trajstream);  // (timestep)
    get_next_line(trajstream);  // ITEM: NUMBER OF ATOMS
    
    // Read the number of atoms
    
    natoms = stoi(get_next_line(trajstream));

    // Read the box dimensions; assume they start from 0,0,0
    
    get_next_line(trajstream);  // ITEM: BOX BOUNDS pp pp pp
    
    boxdims.x = stod(get_token(get_next_line(trajstream), 1));
    boxdims.y = stod(get_token(get_next_line(trajstream), 1));
    boxdims.z = stod(get_token(get_next_line(trajstream), 1));
    
    if(!called_before)
    {
        cout << "Frame has " << natoms << " atoms " << endl;
        cout << "Frame has box dimensions: " << boxdims.x << " " << boxdims.y << " " << boxdims.z << endl;
    }
    
    // Calculate density
    
    density = natoms/boxdims.x/boxdims.y/boxdims.z;
    
    // Read the configuration coordinates
    
    get_next_line(trajstream);  // ITEM: ATOMS id type element xu yu zu
    
    coords.clear();
    xyz     coordinate;
    
    for (int i=0; i<natoms; i++)
    {
        line = get_next_line(trajstream);
        
        //cout << "   ...read line: " << line << endl;
        
        coordinate.x = stod(get_token(line,3));
        coordinate.y = stod(get_token(line,4));
        coordinate.z = stod(get_token(line,5));
        
        coords.push_back(coordinate);
        
        //cout << coords[i].x << " " << coords[i].y << " " << coords[i].z << endl;

    }
    //cout << "=============" << coords.size << endl;
    
    if (!called_before)
        called_before = true;
    
    return true;
}

double configuration::get_dist(int i, int j)
{
    // Compute and return the minimum image convention distance between a pair of particles i and j
    xyz diff;

    // Compute the difference in coordinates between atom i and j
    diff.x = coords[i].x - coords[j].x;
    diff.y = coords[i].y - coords[j].y;
    diff.z = coords[i].z - coords[j].z;

    // Apply minimum image convention to handle periodic boundary conditions
    diff.x -= boxdims.x * round(diff.x / boxdims.x);
    diff.y -= boxdims.y * round(diff.y / boxdims.y);
    diff.z -= boxdims.z * round(diff.z / boxdims.z);

    // Return Euclidean distance
    return sqrt(diff.x * diff.x + diff.y * diff.y + diff.z * diff.z);
}

int get_rdf_bin(int nbins, double binw, double dist)
{
    // Determine the shell bin index the current distance falls in.
    int bin = static_cast<int>(dist / binw);  // Bin index based on distance
    if (bin >= 0 && bin < nbins)
        return bin;
    else
        return -1;
}

int main(int argc, char* argv[])
{
    // Set up the data structure for reading/operating on the trajectory file
    
    string trajfile = argv[1]; //"/Users/becky/Library/CloudStorage/OneDrive-Umich/Virtual_Macbook/Documents/Teaching/CHE_496-696-Simulation/Codes/Monte_Carlo-Inline/indep-1.0.4.MC_traj.lammpstrj";
    configuration frame(trajfile);
    
    
    // Set up the RDF calculation variables
    
    double  binw = stod(argv[2]);   //0.2; // Binwidth for RDF computation - determined based on boxlength
    int     nbins;                  // Number of bins (shells) for RDF compuation
    
    vector<int> num_int;    // Stores histogram of atoms per shell - later used for RDF and number integral calculation
    
    // Read the trajectory file and accumulate RDF data, frame-by-frame
    
    int nframes = 1;
    
    while(frame.read_frame()) // Read and populate the number integral
    {
    
        if (nframes == 1) // Then set up RDF
        {
            nbins = static_cast<int>(frame.boxdims.x / (2.0 * binw));/* complete this */; // Set the number of bins based on the system box length and the requested bin width. 
            num_int.resize(nbins,0);
            
            cout << "Setting nbins, binw to: " << nbins << " " << binw << endl; 
        }
        
        nframes++;
        
        // Populate the single frame
        
        for (int i=0; i<frame.natoms; i++)
        {
            for (int j=i+1; j<frame.natoms; j++)
            {
                double dist = frame.get_dist(i,j);
                int bin = get_rdf_bin(nbins, binw, dist);
                
                // Determine which shell the current atom pair distance is in
                
                if (bin != -1/*complete this criteria for whether the current distance should be added to the num_int histogram*/)
                    num_int[bin] += 2; // Add two to account for both i and j
            }
        }
    }
    
    // Print out the RDF & number integral
    
    
    cout << "# dist, nint g(r)" << endl;
    
    double sanity = 0;
    
    ofstream rdfstream;
    rdfstream.open("rdf.dat");
    
    rdfstream << "# Distance (Ang)      Number Integral (atoms)     RDF (unitless)" << endl;
    
     double avg_nint    = 0;
     double ideal_nint  = 0;
    
     for (int i=0; i<nbins-1; i++)
     {
         double radius              = (i + 0.5) * binw/* complete this */;  
         double shell_vol           = (4.0 / 3.0) * pi * (pow((radius + binw), 3) - pow(radius, 3))/* complete this */;
         double shell_density       = num_int[i] / (frame.natoms) / shell_vol / nframes/* complete this */;
         double shell_density_ideal = frame.natoms / (frame.boxdims.x * frame.boxdims.y * frame.boxdims.z)/* complete this */;
         double avg_rdf             = shell_density / shell_density_ideal/* complete this */;

         avg_nint           += shell_density;

         rdfstream << radius << " " << avg_nint << " " << avg_rdf  << endl;
         
     }
    
    rdfstream.close();
    
    
    
    
    

}


    
    