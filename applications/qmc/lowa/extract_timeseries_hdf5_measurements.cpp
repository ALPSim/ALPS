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
  std::ifstream inFile;     std::string fileIN;
  std::ofstream outFile;    std::string fileOUT;

  std::string dim_str;       int dim;
  std::string len_str;       int length;
  std::string label_lb_str;  int label_lb;
  std::string label_ub_str;  int label_ub;
  std::string thermal_str;   int thermal  = 0;
  std::string skip_str;      int skip     = 1;
  std::string description_str;
 

  int optchar;

  while ((optchar = getopt (argc, argv, "i:o:l:u:t:s:d:")) != -1)
  {
    switch (optchar)
    {
      case 'i':
        fileIN  = (std::string) strdup (optarg);
        break;
      case 'o':
        fileOUT  = (std::string) strdup (optarg);
        break;

      case 'l':
        label_lb_str = (std::string) strdup (optarg);
        break;
      case 'u':
        label_ub_str = (std::string) strdup (optarg);
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

  if (!label_lb_str.empty())
  {
    std::istringstream ss(label_lb_str);
    ss >> label_lb;
  }

  if (!label_ub_str.empty())
  {
    std::istringstream ss(label_ub_str);
    ss >> label_ub;
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


  bool is_data_inputed = false;


  std::vector<std::valarray<double> > raw_data;

  double dummy;  

  for (int label=label_lb; label <= label_ub; ++label)
  {
    std::ostringstream ss_out;
    ss_out << label;
    std::string cur_filename = fileIN;
    if (label != 0)  {  cur_filename += ss_out.str();  }
    cur_filename += ".h5";

    bool is_hdf5_file_found = false;
    inFile.open(cur_filename.c_str(),std::ios::in);
    if (inFile.good())  {  is_hdf5_file_found = true;  }
    inFile.close();
 
    if (is_hdf5_file_found)
    {
      alps::hdf5::iarchive ia(cur_filename.c_str());
  
      uint32_t sweeps = 0;
      ia >> alps::make_pvp("No_of_datasets",sweeps);
  
      for (uint32_t counter=1; counter <= sweeps; ++counter)
      {
        std::ostringstream ss2_out;
        ss2_out << counter;
        std::string cur_description_str = description_str + ss2_out.str();
  
        std::vector<double> raw_data_elem_vec;
        if (ia.is_data(cur_description_str))
        {
          ia >> alps::make_pvp(cur_description_str,raw_data_elem_vec);
          std::valarray<double> raw_data_elem = alps::numeric::vector2valarray<double>(raw_data_elem_vec);
          if ((counter > thermal) && ((counter % skip) == 0))  
          {  
            raw_data.push_back(raw_data_elem); 
            is_data_inputed = true; 
            std::cout << "Dataset (Label " << label << " , Sweep " << counter << " ) read from hdf5 file... ;  Nsites = " << raw_data_elem_vec.size() << " , total no of particles = " << std::accumulate(raw_data_elem_vec.begin(),raw_data_elem_vec.end(),0.) << std::endl;
          }
        }
      }
    }
  }


  if (is_data_inputed)
  {
    uint32_t Nsites = raw_data[0].size();
    uint32_t count  = raw_data.size();
    alps::RealObservable* data = new alps::RealObservable [Nsites];
    for (uint32_t index=0; index < Nsites; ++index)
    {
      for (uint32_t count_index=0; count_index < count; ++count_index)
      {
        data[index] << raw_data[count_index][index];
      }
    }
    outFile.open(fileOUT.c_str(),std::ios::out);
    for (std::size_t index=0; index < Nsites; ++index)
    {
      outFile << index << "\t" << data[index].mean() << "\t" << data[index].error() << "\n";
    }
    outFile.close();


    delete [] data;
  }

  return 0;
}


