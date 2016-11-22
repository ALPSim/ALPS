#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/models/measurement_parser.hpp"

using namespace std;

int main()
{

DmrgParameters parms;

parms.set("MEASURE_LOCAL[ XHOYX]" , "exchange" );
parms.set("MEASURE_AVERAGE [dfhdgh]" , "fermion_hop" );
parms.set(" MEASURE_LOCAL_AT [ local at ] " , "Sz:Sz|(3,1),(1,2)" );

std::vector<MeasurementSpecs> vec = parse_measurements(parms);

for( std::vector<MeasurementSpecs>::iterator m = vec.begin() ; m != vec.end() ; m++ )
{ 
  cout << "Measurement type is: " << m->meas_class << endl;
  cout << "Measurement name is: " << m->name << endl;
  cout << "Measurement arguments are: " << m->args << endl;
}

return 0;
}
