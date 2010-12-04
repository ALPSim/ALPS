#ifndef LP_SSE
#define LP_SSE

#include <iostream>
#include <vector>

//#define DEBUGMODE

using namespace std;

class LP_Solve_SSE {
public:
  void init(const vector<double>& );
  void calc_solution(vector<double>& );
	double calc_loc_opt_bounce(); // NOT NEEDED FOR Linear Programming !!!!
private:
  void do_pivot(const int, const int);
  bool check_optimal();
	void print_tableau() const;
	void copy_solution(vector<double>& );
	int find_pivot_column() const;
	int find_pivot_row(const int) const;
  int Nvar;
	int Nslack;
	vector<int> mVar;
	vector<double> mWeight;
  vector<double> C; // optimization expression
  vector<vector<double> > A; // Matrix of variables
  vector<double> B; // right hand side of inequalies;
  vector<double> mSolution;  
  double mBounce;
};

#endif
