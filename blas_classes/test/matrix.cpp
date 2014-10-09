#include "../matrix.h"


int main(int argc, char **argv){
  std::cout<<"matrix test: "<<std::endl;

  blas::matrix M(2);
  M(0,0)=1;
  M(0,1)=2;
  M(1,0)=3;
  M(1,1)=4;
  std::cout<<"starting with matrix: "<<M<<std::endl;
  
  blas::vector row(2); row(0)=10; row(1)=11;
  blas::vector col(2); col(0)=20; col(1)=21;
  M.insert_row_column_last(row, col, 30);
  std::cout<<"adding row: "<<row<<" and column "<<col<<" "<<std::endl<<M<<std::endl;

  M.swap_row_column(0,2);
  std::cout<<"swapping first and third column and row: "<<M<<std::endl;

  M.remove_row_column_last();
  std::cout<<"removing last row/col "<<M<<std::endl;

  M.getrow(0,&(row(0)));
  M.getcol(1,&(col(0)));
  std::cout<<"row 0 is: "<<row<<" col 1 is: "<<col<<" of matrix "<<M<<std::endl;
}
  
