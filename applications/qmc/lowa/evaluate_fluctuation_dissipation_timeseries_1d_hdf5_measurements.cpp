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


#define VC 0.04272                // put here because it is


int main(int argc, char** argv)
{
  std::ifstream inFile;     std::string fileIN1, fileIN2, fileIN3;
  std::ofstream outFile;    std::string fileOUT1, fileOUT2, fileOUT3, fileOUT4;

  std::string dim_str;       int dim;
  std::string len_str;       int length;
  std::string label_lb_str;  int label_lb;
  std::string label_ub_str;  int label_ub;
  std::string thermal_str;   int thermal  = 0;
  std::string skip_str;      int skip     = 1;
  std::string description_str;
  std::string binsize_str;   double   binsize = 1.;
  std::string temp_str;      double temp = 0.;
  bool is_print = false;
 

  int optchar;

  while ((optchar = getopt (argc, argv, "i:j:k:e:f:g:h:l:u:t:s:d:b:T:p")) != -1)
  {
    switch (optchar)
    {
      case 'i':
        fileIN1  = (std::string) strdup (optarg);
        break;
      case 'j':
        fileIN2  = (std::string) strdup (optarg);
        break;
      case 'k':
        fileIN3  = (std::string) strdup (optarg);
        break;

      case 'e':
        fileOUT1  = (std::string) strdup (optarg);
        break;
      case 'f':
        fileOUT2  = (std::string) strdup (optarg);
        break;
      case 'g':
        fileOUT3  = (std::string) strdup (optarg);
        break;
      case 'h':
        fileOUT4  = (std::string) strdup (optarg);
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

      case 'b':
        binsize_str = (std::string) strdup (optarg);
        break;

      case 'T':
        temp_str = (std::string) strdup (optarg);
        break;

      case 'p':
        is_print = true;
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

  if (!binsize_str.empty())
  {
    std::istringstream ss(binsize_str);
    ss >> binsize;
  }

  if (!temp_str.empty())
  {
    std::istringstream ss(temp_str);
    ss >> temp;
  }


  bool is_data1_inputed = false;
  bool is_data2_inputed = false;
  bool is_data3_inputed = false;

  std::vector<std::valarray<double> > raw_data_dns;
  std::vector<std::valarray<double> > raw_data_N;
  std::vector<std::valarray<double> > raw_data_dnsN;

  double dummy;  

  for (int label=label_lb; label <= label_ub; ++label)
  {
    std::ostringstream ss_out;
    ss_out << label;

    std::string cur_filename1 = fileIN1;
    if (label != 0)  {  cur_filename1 += ss_out.str();  }
    cur_filename1 += ".h5";

    bool is_hdf5_file1_found = false;
    inFile.open(cur_filename1.c_str(),std::ios::in);
    if (inFile.good())  {  is_hdf5_file1_found = true;  }
    inFile.close();
 
    if (is_hdf5_file1_found)
    {
      alps::hdf5::iarchive ia1(cur_filename1.c_str());
  
      uint32_t sweeps = 0;
      ia1 >> alps::make_pvp("No_of_datasets",sweeps);
  
      for (uint32_t counter=1; counter <= sweeps; ++counter)
      {
        std::ostringstream ss2_out;
        ss2_out << counter;
        std::string cur_description_str = description_str + ss2_out.str();
  
        std::vector<double> raw_data_dns_elem_vec;
        if (ia1.is_data(cur_description_str))
        {
          ia1 >> alps::make_pvp(cur_description_str,raw_data_dns_elem_vec);
          std::valarray<double> raw_data_dns_elem = alps::numeric::vector2valarray<double>(raw_data_dns_elem_vec);
          if ((counter > thermal) && ((counter % skip) == 0))  
          {  
            raw_data_dns.push_back(raw_data_dns_elem); 
            is_data1_inputed = true; 
            //std::cout << "Dataset (Label " << label << " , Sweep " << counter << " ) read from hdf5 file... ;  Nsites = " << raw_data_dns_elem_vec.size() << " , total no of particles = " << std::accumulate(raw_data_dns_elem_vec.begin(),raw_data_dns_elem_vec.end(),0.) << std::endl;
          }
        }
      }
    }
  }

  for (int label=label_lb; label <= label_ub; ++label)
  {
    std::ostringstream ss_out;
    ss_out << label;

    std::string cur_filename2 = fileIN2;
    if (label != 0)  {  cur_filename2 += ss_out.str();  }
    cur_filename2 += ".h5";

    bool is_hdf5_file2_found = false;    inFile.open(cur_filename2.c_str(),std::ios::in);
    if (inFile.good())  {  is_hdf5_file2_found = true;  }
    inFile.close();
 
    if (is_hdf5_file2_found)
    {
      alps::hdf5::iarchive ia2(cur_filename2.c_str());
  
      uint32_t sweeps = 0;
      ia2 >> alps::make_pvp("No_of_datasets",sweeps);
  
      for (uint32_t counter=1; counter <= sweeps; ++counter)
      {
        std::ostringstream ss2_out;
        ss2_out << counter;
        std::string cur_description_str = description_str + ss2_out.str();

        std::vector<double> raw_data_N_elem_vec;
        if (ia2.is_data(cur_description_str))
        {
          ia2 >> alps::make_pvp(cur_description_str,raw_data_N_elem_vec);
          std::valarray<double> raw_data_N_elem = alps::numeric::vector2valarray<double>(raw_data_N_elem_vec);
          if ((counter > thermal) && ((counter % skip) == 0))
          {
            raw_data_N.push_back(raw_data_N_elem);
            is_data2_inputed = true;
            //std::cout << "Dataset (Label " << label << " , Sweep " << counter << " ) read from hdf5 file... ;  Nsites = " << raw_data_N_elem_vec.size() << " , total no of particles = " << std::accumulate(raw_data_N_elem_vec.begin(),raw_data_N_elem_vec.end(),0.) << std::endl;
          }
        }
      }
    }
  }

  for (int label=label_lb; label <= label_ub; ++label)
  {
    std::ostringstream ss_out;
    ss_out << label;

    std::string cur_filename3 = fileIN3;
    if (label != 0)  {  cur_filename3 += ss_out.str();  }
    cur_filename3 += ".h5";

    bool is_hdf5_file3_found = false;    inFile.open(cur_filename3.c_str(),std::ios::in);
    if (inFile.good())  {  is_hdf5_file3_found = true;  }
    inFile.close();
 
    if (is_hdf5_file3_found)
    {
      alps::hdf5::iarchive ia3(cur_filename3.c_str());
  
      uint32_t sweeps = 0;
      ia3 >> alps::make_pvp("No_of_datasets",sweeps);
  
      for (uint32_t counter=1; counter <= sweeps; ++counter)
      {
        std::ostringstream ss2_out;
        ss2_out << counter;
        std::string cur_description_str = description_str + ss2_out.str();

        std::vector<double> raw_data_dnsN_elem_vec;
        if (ia3.is_data(cur_description_str))
        {
          ia3 >> alps::make_pvp(cur_description_str,raw_data_dnsN_elem_vec);
          std::valarray<double> raw_data_dnsN_elem = alps::numeric::vector2valarray<double>(raw_data_dnsN_elem_vec);
          if ((counter > thermal) && ((counter % skip) == 0))
          {
            raw_data_dnsN.push_back(raw_data_dnsN_elem);
            is_data3_inputed = true;
            //std::cout << "Dataset (Label " << label << " , Sweep " << counter << " ) read from hdf5 file... ;  Nsites = " << raw_data_dnsN_elem_vec.size() << " , total no of particles = " << std::accumulate(raw_data_dnsN_elem_vec.begin(),raw_data_dnsN_elem_vec.end(),0.) << std::endl;
          }
        }
      }
    }
  }

 
  if (!is_data1_inputed || !is_data2_inputed || !is_data3_inputed)   {  return (-1);  }

  uint32_t Nsites = raw_data_dns[0].size();
  uint32_t count  = raw_data_dns.size();

  double* coord = new double [Nsites];

  alps::RealObservable* dns  = new alps::RealObservable [Nsites];
  alps::RealObservable* N    = new alps::RealObservable [Nsites];
  alps::RealObservable* dnsN = new alps::RealObservable [Nsites];

  for (uint32_t index=0; index < Nsites; ++index)
  {
    coord[index] = (index+0.5)*binsize;
  }

  for (uint32_t index=0; index < Nsites; ++index)
  {
    for (uint32_t count_index=0; count_index < count; ++count_index)
    {
      dns[index]  << raw_data_dns[count_index][index];
      N[index]    << raw_data_N[count_index][index];
      dnsN[index] << raw_data_dnsN[count_index][index];
    }
  }
  for (std::size_t index=0; index < Nsites; ++index)  {  if (is_print)  {  std::cout << index << " " << dns[index];  }  }
  for (std::size_t index=0; index < Nsites; ++index)  {  if (is_print)  {  std::cout << index << " " << N[index];    }  }
  for (std::size_t index=0; index < Nsites; ++index)  {  if (is_print)  {  std::cout << index << " " << dnsN[index]; }  }

  alps::RealObsevaluator* cur_dns  = new alps::RealObsevaluator [Nsites];
  alps::RealObsevaluator* cur_N    = new alps::RealObsevaluator [Nsites];
  alps::RealObsevaluator* cur_dnsN = new alps::RealObsevaluator [Nsites];

  for (int index=0; index < Nsites; ++index)
  {
    cur_dns[index]  = dns[index];
    cur_N[index]    = N[index];
    cur_dnsN[index] = dnsN[index];
  }

  alps::RealObsevaluator* dn_dr  = new alps::RealObsevaluator [Nsites];
  alps::RealObsevaluator* diss   = new alps::RealObsevaluator [Nsites];
  
  dn_dr[0]  = ((cur_dns[1] - cur_dns[0]) / (2*(coord[1]-coord[0])));   // dns[-1] = dns[0]   (by symmetry)
  for (int index=1; index < Nsites-1; ++index)
  {
    dn_dr[index]  = ((cur_dns[index+1] - cur_dns[index-1]) / (coord[index+1]-coord[index-1]));
  }
  dn_dr[Nsites-1]  = ((cur_dns[Nsites-1] - cur_dns[Nsites-2]) / (2*(coord[Nsites-1]-coord[Nsites-2])));


  alps::RealObsevaluator* fluct = new alps::RealObsevaluator [Nsites];

  for (int index=0; index < Nsites; ++index)  
  {   
    diss[index].rename("");
    diss[index]  = -((dn_dr[index]) / (2 * VC * coord[index]));

    fluct[index].rename("");
    fluct[index]  = cur_dnsN[index] - cur_N[index] * cur_dns[index];  
  }

  std::string cur_filename4 = fileOUT1;
  std::string cur_filename5 = fileOUT2;
  std::string cur_filename6 = fileOUT3;

  cur_filename4 += ".dat";
  cur_filename5 += ".dat";
  cur_filename6 += ".dat";

  outFile.open(cur_filename4.c_str(),std::ios::out);
  for (int index=0; index < Nsites; ++index)
  {
    outFile << coord[index] << "\t" << diss[index].mean() << "\t" << diss[index].error() << "\n";
    std::cout << coord[index] << "\t" << diss[index];
  }
  outFile.close();

  outFile.open(cur_filename5.c_str(),std::ios::out);
  for (int index=0; index < Nsites; ++index)  
  {   
    outFile << coord[index] << "\t" << fluct[index].mean() << "\t" << fluct[index].error() << "\n";
    std::cout << coord[index] << "\t" << fluct[index];  
  }
  outFile.close();

  outFile.open(cur_filename6.c_str(),std::ios::out);
  for (int index=0; index < Nsites; ++index)
  {
    outFile << diss[index].mean() << "\t" << fluct[index].mean() << "\t" << diss[index].error() << "\t" << fluct[index].error() << "\n";
  }
  outFile.close();

  outFile.open(fileOUT4.c_str(),std::ios::out);
  for (int index=0; index < Nsites; ++index)  {  outFile << diss[index].mean() << "\t" << temp * diss[index].mean() << "\n";  }
  outFile.close();


  delete [] coord;

  delete [] dns;
  delete [] N;
  delete [] dnsN;
  delete [] cur_dns;
  delete [] cur_N;
  delete [] cur_dnsN;


  delete [] dn_dr;
  delete [] diss;
  delete [] fluct;

  return 0;
}


