#include "types.h"
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
void print_imag_green_matsubara(std::ostream &os, const multiple_complex_vector_type &v, const double beta){
	for(uint i=0;i<v.first.size();++i)
                os<<(2*i+1)*M_PI/beta<<"\t"<<v.first[i].imag()<<"\t"<<v.second[i].imag()<<std::endl;
}
void print_green_itime(std::ostream &os, const multiple_vector_type &v, const double beta){
        for(uint i=0;i<v.first.size();++i)
                os<<i*beta/v.first.size()<<"\t"<<v.first[i]<<"\t"<<v.second[i]<<std::endl;

}
