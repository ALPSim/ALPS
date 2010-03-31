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

/* $Id: from_3d_to_2d_hdf5_measurements.cpp 3520 2010-03-19 22:00:00Z tamama $ */


#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <utility>
#include <valarray>
#include <vector>
#include <unistd.h>
#include <alps/numeric/special_functions.hpp>
#include <alps/numeric/vector_valarray_conversion.hpp>
#include <alps/hdf5.hpp>
#include <alps/alea.h>


int main(int argc, char** argv)
{
  std::ifstream inFile;     std::string fileIN;
  std::ofstream outFile;    std::string fileOUT;

  std::string len_str;       int length;
  std::string thermal_str;   int thermal  = 0;
  std::string skip_str;      int skip     = 1;
  std::string description_str;
 

  int optchar;

  while ((optchar = getopt (argc, argv, "i:o:t:s:d:")) != -1)
  {
    switch (optchar)
    {
      case 'i':
        fileIN  = (std::string) strdup (optarg);
        break;

      case 'o':
        fileOUT  = (std::string) strdup (optarg);
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


  std::vector<std::valarray<uint8_t> > raw_data;

  uint8_t dummy;  

  std::string cur_filename = fileIN;
  cur_filename += ".h5";

  bool is_hdf5_file_found = false;

  inFile.open(cur_filename.c_str(),std::ios::in);
  if (inFile.good())  {  is_hdf5_file_found = true;  }
  inFile.close();

  uint32_t sweeps =0;

  fileOUT += ".h5";
  alps::hdf5::oarchive oa(fileOUT.c_str());

 
  if (is_hdf5_file_found)
  {
    alps::hdf5::iarchive ia(cur_filename.c_str());

    ia >> alps::make_pvp("No_of_datasets",sweeps);

    for (uint32_t counter=1; counter <= sweeps; ++counter)
    {
      std::ostringstream ss2_out;
      ss2_out << counter;
      std::string cur_description_str = description_str + ss2_out.str();

      std::vector<uint8_t> raw_data_elem_vec;
      if (ia.is_data(cur_description_str))
      {
        ia >> alps::make_pvp(cur_description_str,raw_data_elem_vec);

        uint32_t L    = alps::numeric::cbrt(raw_data_elem_vec.size());
        uint32_t L_sq = alps::numeric::sq(L);
     
        std::valarray<uint8_t> raw_data_elem = alps::numeric::vector2valarray<uint8_t>(raw_data_elem_vec);
        if ((counter > thermal) && ((counter % skip) == 0))  
        {  
          raw_data.push_back(raw_data_elem); 
          oa << alps::make_pvp(cur_description_str, &raw_data_elem[0], raw_data_elem_vec.size());
          oa << alps::make_pvp("No_of_datasets",counter);
          std::cout << "Dataset (Sweep " << counter << " ) read from hdf5 file... ;  Nsites = " << raw_data_elem_vec.size() << " ; Nsites (projected) = " << L_sq << std::endl;
        }
      }
    }
  }


  return 0;
}


