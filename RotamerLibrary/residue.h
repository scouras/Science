struct minors {
  int     symmetry;         // Degree of symmetry
  int     start;            // Starting angle
  int     stop;             // Ending angle
  int     number;           // 1-Index of bin, taking into account symmetry
  String*  conformation;     // Name of bin
  bool    crosses_0;        // Whether the bin crosses 0
  bool    crosses_180;      // Whether the bin crosses 180
  int     min_start;        // Min start of symmetric bins (mod 360)
}

struct dihedral {
  int       index;          // Index in residue
  int       count;          // Number of bins (reducing by symmetry)
  char[4]*  atoms;          // Four atom names defining DH
  int*      starts;         // Starting angle of each minor bin
  int       symmetry;       // Symmetry for this angle
  minors*   minors;         // Definition of each minor bin
  int       multiplier;     // Multiplier uses when calculating major rotamer index
}

struct residue {
  String  name;               // Name of the amino acid
  String  ab;                 // Amino Acid one letter abbreviation
  int     count;              // Count of this residue
  int     precount;           // Initial count of residues
  int     total_rotamers;     // Count of major rotamers

  String* angle_names;        // Names of Dihedral Angles
  String* all;                // Names of Dihedral Angles inhereted from Utility.pm
  String* mc;                 // ^-- for Main Chain
  String* sc;                 // ^-- for Side Chain
  String* h;                  // ^-- for Hydrogen DH's
  String* o;                  // ^-- for Other DH's
  
  String* mc_atoms;            // Main Chain atoms from Utility.pm
  String* sc_atoms;            // Side Chain atoms from Utility.pm

  hash{String=>int} angle_indecies;       // Index of angles in angle list
  int   omega_index;                      // Index of the Phi angle in angle list
  int   phi_index;                        // Index of the Phi angle in angle list
  int   psi_index;                        // Index of the Psi angle in angle list
  int** rotamer_indecies;                 // List of indecies for rotamer angles
  int** quick_numbers;                    // Lookup for minor number
  int** quick_starts;                     // Lookup for minor start angles
  int*  quick_mults;                      // Lookup for multipliers for major rotamer 


  hash{String=>int}       sasa_count;     // Count of residues in each sasa bin

  hash{String=>dihedral}  dihedrals;      // Dihedral name to dihedral struct

  int**                   major_to_minor; // Lookup by major list of minors
  hash{String=>int}       minor_to_major; // { '1_1' => 1 }
  hash{String=>[String]}  conformations;  // { chi1=>[g+,t,g-] }
  hash{String=>int}       sym;            // { angle_name => symmetry } if > 1 }
  
  bool    initialized_rotamers;   // Rotamers are done being read and initialized
  bool    finalized_rotamers;     // Rotamers are done with optimizations


};
