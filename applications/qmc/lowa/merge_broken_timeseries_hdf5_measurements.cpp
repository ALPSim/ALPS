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

/* $Id: evaluate_statistics_of_timeseries_hdf5_density_measurements.cpp 3520 2010-03-19 22:00:00Z tamama $ */


#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <utility>
#include <valarray>
#include <vector>
#include <unistd.h>
#include <alps/numeric/vector_valarray_conversion.hpp>
#include <alps/hdf5.hpp>
#include <alps/alea.h>


int main(int argc, char** argv)
{
  std::ifstream inFile;     std::string fileIN1, fileIN2;
  std::ofstream outFile;    std::string fileOUT;

  std::string len_str;       int length;
  std::string thermal_str;   int thermal  = 0;
  std::string skip_str;      int skip     = 1;
  std::string description_str;
 

  int optchar;

  while ((optchar = getopt (argc, argv, "a:b:o:t:s:d:")) != -1)
  {
    switch (optchar)
    {
      case 'a':
        fileIN1  = (std::string) strdup (optarg);
        break;
      case 'b':
        fileIN2  = (std::string) strdup (optarg);
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

  std::string cur_filename1 = fileIN1;
  cur_filename1 += ".h5";

  std::string cur_filename2 = fileIN2;
  cur_filename2 += ".h5";

  bool is_hdf5_file1_found = false;
  bool is_hdf5_file2_found = false;

  inFile.open(cur_filename1.c_str(),std::ios::in);
  if (inFile.good())  {  is_hdf5_file1_found = true;  }
  inFile.close();

  inFile.open(cur_filename2.c_str(),std::ios::in);
  if (inFile.good())  {  is_hdf5_file2_found = true;  }
  inFile.close();

  uint32_t sweeps1 =0;
  uint32_t sweeps2 =0;

  fileOUT += ".h5";
  alps::hdf5::oarchive oa(fileOUT.c_str());

 
  if (is_hdf5_file1_found)
  {
    alps::hdf5::iarchive ia1(cur_filename1.c_str());

    ia1 >> alps::make_pvp("No_of_datasets",sweeps1);

    for (uint32_t counter=1; counter <= sweeps1; ++counter)
    {
      std::ostringstream ss2_out;
      ss2_out << counter;
      std::string cur_description_str = description_str + ss2_out.str();

      std::vector<uint8_t> raw_data_elem_vec;
      if (ia1.is_data(cur_description_str))
      {
        ia1 >> alps::make_pvp(cur_description_str,raw_data_elem_vec);
        std::valarray<uint8_t> raw_data_elem = alps::numeric::vector2valarray<uint8_t>(raw_data_elem_vec);
        if ((counter > thermal) && ((counter % skip) == 0))  
        {  
          raw_data.push_back(raw_data_elem); 
          oa << alps::make_pvp(cur_description_str, &raw_data_elem[0], raw_data_elem_vec.size());
          oa << alps::make_pvp("No_of_datasets",counter);
          std::cout << "Dataset (Sweep " << counter << " ) read from hdf5 file... ;  Nsites = " << raw_data_elem_vec.size() << std::endl;
        }
      }
    }
  }

  if (is_hdf5_file2_found)
  {
    alps::hdf5::iarchive ia2(cur_filename2.c_str());

    ia2 >> alps::make_pvp("No_of_datasets",sweeps2);

    for (uint32_t counter=sweeps1+1; counter <= sweeps2; ++counter)
    {   
      std::ostringstream ss2_out;
      ss2_out << counter;
      std::string cur_description_str = description_str + ss2_out.str();

      std::vector<uint8_t> raw_data_elem_vec;
      if (ia2.is_data(cur_description_str))
      {   
        ia2 >> alps::make_pvp(cur_description_str,raw_data_elem_vec);
        std::valarray<uint8_t> raw_data_elem = alps::numeric::vector2valarray<uint8_t>(raw_data_elem_vec);        if ((counter > thermal) && ((counter % skip) == 0))  
        {   
          raw_data.push_back(raw_data_elem); 
          oa << alps::make_pvp(cur_description_str, &raw_data_elem[0], raw_data_elem_vec.size());
          oa << alps::make_pvp("No_of_datasets",counter);
          std::cout << "Dataset (Sweep " << counter << " ) read from hdf5 file... ;  Nsites = " << raw_data_elem_vec.size() << std::endl;
        }   
      }   
    }   
  }


  return 0;
}


