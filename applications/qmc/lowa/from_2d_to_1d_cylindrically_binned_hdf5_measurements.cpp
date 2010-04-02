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
#include <algorithm>
#include <numeric>
#include <unistd.h>
#include <alps/numeric/special_functions.hpp>
#include <alps/numeric/vector_valarray_conversion.hpp>
#include <alps/hdf5.hpp>
#include <alps/alea.h>



template <class T1,class T2>
void from_2dvec_to_1dvec(std::vector<T1>& from, std::vector<T2>& to, std::vector<T1>& to_index, std::vector<T1>& freq, T1 const from_size, T1 const to_size)
{
  to.clear();
  to.resize(to_size,0.);     

  for (T1 p=0; p < from_size; ++p)  {  to[to_index[p]] += from[p];  }
  for (T1 r=0; r < to_size; ++r)    {  to[r] /= freq[r];  }
}



int main(int argc, char** argv)
{
  std::ifstream inFile;     std::string fileIN;
  std::ofstream outFile;    std::string fileOUT;

  std::string thermal_str;   int thermal  = 0;
  std::string skip_str;      int skip     = 1;
  std::string description_str;
  std::string binsize_str;   double   binsize = 1.;
  std::string len_str;       uint32_t L;
  bool is_print = false;
 

  int optchar;

  while ((optchar = getopt (argc, argv, "i:o:t:s:d:b:L:p")) != -1)
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

      case 'b':
        binsize_str = (std::string) strdup (optarg);
        break;
      case 'L':
        len_str = (std::string) strdup (optarg);
        break;
 
      case 'p':
        is_print = true;
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

  if (!len_str.empty())
  {
    std::istringstream ss(len_str);
    ss >> L;
  }

  if (!binsize_str.empty())
  {
    std::istringstream ss(binsize_str);
    ss >> binsize;
  }


  uint32_t L_sq  = L*L;  
  double   xmid  = (static_cast<double>(L-1))/2;
  double   ymid  = (static_cast<double>(L-1))/2;
  uint32_t nrbin = std::ceil(std::sqrt((L-1-xmid)*(L-1-xmid) + (L-1-ymid)*(L-1-ymid))/binsize);

  std::vector<uint32_t> proj_binnr;
  std::vector<uint32_t> proj_binfreq(nrbin,0);

  for (uint32_t j=0; j < L; ++j) {
    for (uint32_t i=0; i < L; ++i) {
      uint32_t cur_binnr = std::floor(std::sqrt((i-xmid)*(i-xmid) + (j-ymid)*(j-ymid))/binsize);
      proj_binnr.push_back(cur_binnr);
      proj_binfreq[cur_binnr] += 1;
    }
  }

  if (is_print) {
    std::cout << "\n\n";
    std::cout << "r" << "\t" << "freq" << "\n";
    for (uint32_t r=0; r < nrbin; ++r) {
      std::cout << ((static_cast<double>(r)+0.5)*binsize) << "\t" << proj_binfreq[r] << "\n";
    }
    std::cout << "\n";
    std::cout << "Cummulative freq = " << std::accumulate(proj_binfreq.begin(),proj_binfreq.end(),0) << "\n";
    std::cout << "\n";
  }


  uint32_t dummy;  

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

      std::vector<uint32_t> raw_data_elem_vec;
      if (ia.is_data(cur_description_str))
      {
        ia >> alps::make_pvp(cur_description_str,raw_data_elem_vec);

        std::vector<double>   raw_data_elem_1d_vec;

        from_2dvec_to_1dvec<uint32_t,double>(raw_data_elem_vec,raw_data_elem_1d_vec,proj_binnr,proj_binfreq,L_sq,nrbin);
 
        std::valarray<double> raw_data_elem_1d = alps::numeric::vector2valarray<double>(raw_data_elem_1d_vec);
        if ((counter > thermal) && ((counter % skip) == 0))  
        {  
          oa << alps::make_pvp(cur_description_str, &raw_data_elem_1d[0], nrbin);
          oa << alps::make_pvp("No_of_datasets",counter);

          //std::cout << "Dataset (Sweep " << counter << " ) read from hdf5 file... ;  Nsites (projected) = " << raw_data_elem_vec.size() << " ; L^2 = " << L_sq << std::endl;
        }
      }
    }
  }

  return 0;
}


