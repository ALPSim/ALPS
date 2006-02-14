#include<iostream>
#include<fstream>
#include<cmath>
#include<cassert>

int main(int argc, char** argv){
	
	assert(argc==4); //call: ./difftest file1 file2 maxdiff

	double maxdiff=atof(argv[3]);
	std::ifstream file_1(argv[1]);
	std::ifstream file_2(argv[2]);

	do{
		double d1;
		double d2;
		file_1>>d1;
		file_2>>d2;
		if(fabs(d1-d2) > maxdiff){	//test failed
			std::cerr<<"test failed - found difference of "<<maxdiff<<std::endl;
			return(1);	//exit with error
		}
	}while(!file_1.eof());
	std::cout<<"test succeeded"<<std::endl;
	return 0;
}
