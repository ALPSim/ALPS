#include "types.h"
#include "green_function.h"
#include <fstream>
#include <sstream>
std::ostream &operator<<(std::ostream &os, const multiple_vector_type &v){
	for(uint i=0;i<v.first.size();++i)
		os<<i<<"\t"<<v.first[i]<<"\t"<<v.second[i]<<std::endl;
	return os;
}
std::ostream &operator<<(std::ostream &os, const multiple_complex_vector_type &v){
	for(uint i=0;i<v.first.size();++i)
		//os<<i<<"\t"<<v.first[i].real()<<"\t"<<v.first[i].imag()<<"\t"<<v.second[i].real()<<"\t"<<v.second[i].imag()<<std::endl;
		os<<i<<"\t"<<v.first[i].imag()<<"\t"<<v.second[i].imag()<<std::endl;
	return os;
}
std::ostream &operator<<(std::ostream &os, const double_vector &v){
	for(uint i=0;i<v.size();++i)
		//os<<i<<"\t"<<v.first[i].real()<<"\t"<<v.first[i].imag()<<"\t"<<v.second[i].real()<<"\t"<<v.second[i].imag()<<std::endl;
		os<<i<<"\t"<<v[i]<<std::endl;
	return os;
}
void print_imag_green_matsubara(std::ostream &os, const multiple_complex_vector_type &v, const double beta){
	for(uint i=0;i<v.first.size();++i)
                os<<(2*i+1)*M_PI/beta<<"\t"<<v.first[i].imag()<<"\t"<<v.second[i].imag()<<std::endl;
}
void print_real_green_matsubara(std::ostream &os, const multiple_complex_vector_type &v, const double beta){
	for(uint i=0;i<v.first.size();++i)
                os<<(2*i+1)*M_PI/beta<<"\t"<<v.first[i].real()<<"\t"<<v.second[i].real()<<std::endl;
}
void print_green_itime(std::ostream &os, const multiple_vector_type &v, const double beta){
        for(uint i=0;i<v.first.size();++i)
                os<<i*beta/(v.first.size()-1)<<"\t"<<v.first[i]<<"\t"<<v.second[i]<<std::endl;

}
void print_all_green_functions(const int iteration_ctr, const multiple_complex_vector_type &G0_omega, const multiple_complex_vector_type &G_omega, const multiple_vector_type &G0_tau, const multiple_vector_type &G_tau, const double beta){
	std::ostringstream G0omega_name, G0tau_name, Gomega_name, Gomegareal_name, Gtau_name, Gtau_name_2, Gomega_name_2;
	G0omega_name<<"G0_omega_"<<iteration_ctr;//<<"_"<<process_id;
	G0tau_name<<"G0_tau_"<<iteration_ctr;//<<"_"<<process_id;
	Gomega_name<<"G_omega_"<<iteration_ctr;//<<"_"<<process_id;
	Gomegareal_name<<"G_omegareal_"<<iteration_ctr;//<<"_"<<process_id;
	Gtau_name<<"G_tau_"<<iteration_ctr;//<<"_"<<process_id;
	Gomega_name_2<<"G_omega";//<<"_"<<process_id;
	Gtau_name_2<<"G_tau";//<<"_"<<process_id;
	std::ofstream G0_omega_file(G0omega_name.str().c_str());
	std::ofstream G0_tau_file(G0tau_name.str().c_str());
	std::ofstream G_omega_file(Gomega_name.str().c_str());
	std::ofstream G_omegareal_file(Gomegareal_name.str().c_str());
	std::ofstream G_tau_file(Gtau_name.str().c_str());
	std::ofstream G_omega_file_2(Gomega_name_2.str().c_str());
	std::ofstream G_tau_file_2(Gtau_name_2.str().c_str());
	assert(G0_omega_file.is_open() && G0_tau_file.is_open() && G_omega_file.is_open() && G_tau_file.is_open());
	print_imag_green_matsubara(G0_omega_file, G0_omega,beta);
	print_green_itime(G0_tau_file,G0_tau,beta);
	print_imag_green_matsubara(G_omega_file,G_omega,beta);
	print_real_green_matsubara(G_omegareal_file,G_omega,beta);
	print_green_itime(G_tau_file,G_tau,beta);
	print_imag_green_matsubara(G_omega_file_2,G_omega,beta);
	print_green_itime(G_tau_file_2,G_tau,beta);
}
void print_tau_green_functions(const int iteration_ctr, const multiple_vector_type &G0_tau, const multiple_vector_type &G_tau, const double beta){
	std::ostringstream G0tau_name, Gtau_name, Gtau_name_2;
	G0tau_name<<"G0_tau_"<<iteration_ctr;//<<"_"<<process_id;
	Gtau_name<<"G_tau_"<<iteration_ctr;//<<"_"<<process_id;
	Gtau_name_2<<"G_tau";//<<"_"<<process_id;
	std::ofstream G0_tau_file(G0tau_name.str().c_str());
	std::ofstream G_tau_file(Gtau_name.str().c_str());
	std::ofstream G_tau_file_2(Gtau_name_2.str().c_str());
	print_green_itime(G0_tau_file,G0_tau,beta);
	print_green_itime(G_tau_file,G_tau,beta);
	print_green_itime(G_tau_file_2,G_tau,beta);
}
std::ostream &operator<<(std::ostream &os, const green_function<double> &v){
        for(uint o=0;o<v.nfreq();++o){
         os<<o<<"\t";
         for(uint i=0;i<v.nsite();++i){
          for(uint j=0;j<v.nsite();++j){
           for(uint z=0;z<v.nflavor();++z){
                os<<v(o,i,j,z)<<"\t";
           }
          }
         }
         std::cout<<std::endl;
        }
        return os;
}
std::ostream &operator<<(std::ostream &os, const green_function<std::complex<double> > &v){
        for(uint o=0;o<v.nfreq();++o){
         os<<o<<"\t";
         for(uint i=0;i<v.nsite();++i){
          for(uint j=0;j<v.nsite();++j){
           for(uint z=0;z<v.nflavor();++z){
                os<<v(o,i,j,z).imag()<<"\t";
           }
          }
         }
         std::cout<<std::endl;
        }
        return os;
}

