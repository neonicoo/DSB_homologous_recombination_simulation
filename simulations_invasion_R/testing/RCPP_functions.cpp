//RCPP_functions.cpp
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector find_occupancies_remastered(const IntegerVector occupiedRAD51, const IntegerVector additional_removals, 
                                          int lower_window = 1,  int upper_window = 0, int n = 2069){
  IntegerVector indices ;
  IntegerVector remove (n);
  int x1;
  int x2;
  int i;
  int j;
  
  for (i = 0; i < occupiedRAD51.size(); i++){
    x1 = occupiedRAD51[i]-7;
    x2 = occupiedRAD51[i]+7;
    if(x1 < 1){x1 = 1;}
    for(j = x1; j<= x2; j++){
      remove[j] = j;
    }
  }
  
  if(additional_removals[0]!=0){
    for (i = 0; i < additional_removals.size(); i++){
      x1 = additional_removals[i]-7;
      x2 = additional_removals[i]+7;
      if(x1 < 1){x1 = 1;}
      for(j = x1; j<= x2; j++){
        remove[j] = j;
      }
    }
  }
  
  if(upper_window == 0 || upper_window > n-7){upper_window = n-7;}
  for(int w = lower_window; w <= upper_window; w++){
      if(remove[w] == 0){
        indices.push_back(w);
    }
  }
  return(indices);
}