#include "band_structure_calculations/band_structure_calculations_BHM.h"

#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>

#define PI 3.141592654


int main() {

   // Rb-87 (Immanuel Bloch)
   int    label    = 1;        	// a label for reference
   int    L        = 100;       // L*L*L lattice
   double mass     = 86.99;  // in u
   double a_s      = 98.49;     	// in bohr radius
   double V0_x     = 11.75;      	// in Er_x
   double V0_y     = 11.75;      	// ...
   double V0_z     = 11.75;      	// ...
   double lambda_x = 843.;     // in nm
   double lambda_y = 765.;     // ...
   double lambda_z = 765.;     // ...
   double waist_x  = 10000050.;     	// in um
   double waist_y  = 10000050.;     	// in um
   double waist_z  = 10000050.;     	// in um
   double VT_x     = 75.;			  // in Hz
   double VT_y     = 75.;  			// in Hz
   double VT_z     = 75.;  			// in Hz

/*
   // Cs-133 (Chin Cheng)
   int    label    = 2;        // a label for reference
   int    L        = 51;       // L*L*L lattice
   double mass     = 133.;      // in u
   double a_s      = 310.;     // in bohr radius
   double V0_x     = 10.;      // in Er_x
   double V0_y     = 10.;      // ...
   double V0_z     = 324.;      // ...
   double lambda_x = 1064.;     // in nm
   double lambda_y = 1064.;     // ...
   double lambda_z = 8000.;     // ...
   double waist_x  = 150.;     // in um
   double waist_y  = 150.;     // in um
   double waist_z  = 150.;     // in um
   double VT_x     = 0.01;     // in nK
   double VT_y     = 0.01;     // in nK
   double VT_z     = 0.01;     // in nK
*/

	 typedef alps::applications::band_structure_calculations_BHM<int, double> BH_PARM;

   BH_PARM BH_parm(label, L, V0_x, V0_y, V0_z, a_s, mass, lambda_x, lambda_y, lambda_z, waist_x, waist_y, waist_z, VT_x, VT_y, VT_z);

   return 0;
}



