#ifndef ALIAS_H_
#define ALIAS_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

typedef struct cell {
    int alias;
    int index;
    double cutoff;
} Cell;

class Alias {
private:
  vector<Cell> cells;
  void computeTwoPointDistributions();
public:
  Alias() {}
  Alias(vector<double>&, vector<int>&);
  virtual ~Alias() {}
  void initiate(vector<double>&, vector<int>&);
  int pick();
  void printTwoPointDistributions();
};

#endif /* ALIAS_H_ */
