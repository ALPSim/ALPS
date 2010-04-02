/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Si_mulations
*
* ALPS Applications
*
* Copyright (C) 2010-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/* $Id: evaluate_observable_observable_xi_coupling_timeseries_2d_hdf5_measurements.cpp 3520 2010-03-19 22:00:00Z tamama $ */


#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <utility>
#include <valarray>
#include <vector>
#include <unistd.h>
#include <alps/numeric/vector_functions.hpp>
#include <alps/numeric/vector_valarray_conversion.hpp>
#include <alps/hdf5.hpp>
#include <alps/alea.h>


template <class T1,class T2>
void generate_timeseries_coupling_obs(std::vector<T1>& o,std::vector<T1>& O, std::vector<T1>& oO, T2 const xi, T2 const L, T2 const L_sq)
{
  O.clear();
  O.resize(L_sq,0);

  for (T2 p=0, j=0; j < L; ++j) {
    for (T2 i=0; i < L; ++i, ++p) {

      for (T2 delta_j=-xi; delta_j <= xi; ++delta_j) {
        for (T2 delta_i=-xi; delta_i <= xi; ++delta_i) {
          if (std::sqrt(delta_i*delta_i + delta_j*delta_j) <= xi) {
            T2 cur_i = i + delta_i;
            T2 cur_j = j + delta_j;
            T2 cur_p = cur_j*L + cur_i;

            if ((cur_i >= 0) && (cur_i < L) && (cur_j >=0) && (cur_j < L))   {  O[p] += o[cur_p];  } 
          }
        }
      }

    }
  }

  using namespace alps::numeric;
  oO = o * O;

/*
  // for debugging purpose, not needed now since totally debugged...
  using std::operator<<;  
  std::cout << "\n\nChecking the validity of this function...\n\n";
  for (T2 i=0; i < L; ++i) {
    T2 p = (L/2)*L + i;
    std::cout << i << "\t" << o[p] << "\t" << o[p-1] << "\t" << o[p+1] << "\t" << o[p-L] << "\t" << o[p+L] << "\t" << O[p] << "\t" << oO[p] << "\n";
  }
  std::cin.get();
*/
}



int main(int argc, char** argv)
{
  std::ifstream inFile;     std::string fileIN;
  std::ofstream outFile;    std::string fileOUT1, fileOUT2;

  std::string thermal_str;   int thermal  = 0;
  std::string skip_str;      int skip     = 1;
  std::string description_str;
  std::string xi_str;        int32_t xi= 0;


  int optchar;

  while ((optchar = getopt (argc, argv, "i:a:b:t:s:d:x:")) != -1)
  {
    switch (optchar)
    {
      case 'i':
        fileIN  = (std::string) strdup (optarg);
        break;

      case 'a':
        fileOUT1 = (std::string) strdup (optarg);
        break;
      case 'b':
        fileOUT2 = (std::string) strdup (optarg);
        break;


      case 't':
        thermal_str = (std::string) strdup (optarg);
        break;
      case 's':
        skip_str = (std::string) strdup (optarg);
        break;

      case 'd':
        description_str = (std::string) strdup (optarg);
        break;

      case 'x':
        xi_str = (std::string) strdup (optarg);
        break;
    }
  }

  if (!thermal_str.empty())
  {
    std::istringstream ss(thermal_str);
    ss >> thermal;
  }

  if (!skip_str.empty())   
  {
    std::istringstream ss(skip_str);
    ss >> skip;
  }

  if (!xi_str.empty())
  {
    std::istringstream ss(xi_str);
    ss >> xi;
  }
  

  uint32_t dummy;  

  std::string cur_filename = fileIN;
  cur_filename += ".h5";

  bool is_hdf5_file_found = false;

  inFile.open(cur_filename.c_str(),std::ios::in);
  if (inFile.good())  {  is_hdf5_file_found = true;  }
  inFile.close();

  uint32_t sweeps =0; 

  fileOUT1 += ".h5";
  alps::hdf5::oarchive oa1(fileOUT1.c_str());

  fileOUT2 += ".h5";
  alps::hdf5::oarchive oa2(fileOUT2.c_str());


  if (is_hdf5_file_found)
  {
    alps::hdf5::iarchive ia(cur_filename.c_str());

    ia >> alps::make_pvp("No_of_datasets",sweeps);

    for (uint32_t counter=1; counter <= sweeps; ++counter)
    {
      std::ostringstream ss2_out;
      ss2_out << counter;
      std::string cur_description_str = description_str + ss2_out.str();

      std::vector<uint32_t> raw_data_elem_o_vec;
      if (ia.is_data(cur_description_str))
      {
        ia >> alps::make_pvp(cur_description_str,raw_data_elem_o_vec);

        int32_t L_sq = raw_data_elem_o_vec.size();
        int32_t L    = std::sqrt(L_sq);

        std::vector<uint32_t> raw_data_elem_O_vec;
        std::vector<uint32_t> raw_data_elem_oO_vec;

        generate_timeseries_coupling_obs<uint32_t,int32_t>(raw_data_elem_o_vec,raw_data_elem_O_vec,raw_data_elem_oO_vec,xi,L,L_sq);


        std::valarray<uint32_t> raw_data_elem_O  = alps::numeric::vector2valarray<uint32_t>(raw_data_elem_O_vec);
        std::valarray<uint32_t> raw_data_elem_oO = alps::numeric::vector2valarray<uint32_t>(raw_data_elem_oO_vec);
        if ((counter > thermal) && ((counter % skip) == 0))
        {
          oa1 << alps::make_pvp(cur_description_str, &raw_data_elem_O[0], L_sq);
          oa1 << alps::make_pvp("No_of_datasets",counter);

          oa2 << alps::make_pvp(cur_description_str, &raw_data_elem_oO[0], L_sq);
          oa2 << alps::make_pvp("No_of_datasets",counter);

          //std::cout << "Dataset (Sweep " << counter << " ) read from hdf5 file... ;  Nsites = " << raw_data_elem_o_vec.size() << std::endl;
        }

      }
    }
  }

  return 0;
}
