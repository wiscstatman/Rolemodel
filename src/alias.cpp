#include <iostream>
#include<iomanip>

#include <R.h>
#include <R_ext/Print.h>
#include "alias.h"

using namespace std;

//// ps - non-zero probabilities for topologies
//// is - corresponding indices
Alias::Alias(vector<double>& ps, vector<int>& is) {
  initiate(ps,is);
}

void Alias::initiate(vector<double>& ps, vector<int>& is) {
  int index = 0;
  int nCells = ps.size();
  cells.resize(nCells);
  for(vector<double>::iterator pit = ps.begin(); pit != ps.end(); ++pit, ++index) {
    cells[index].cutoff = *pit * nCells;
    cells[index].alias = is[index];
    cells[index].index = is[index];
  }
  computeTwoPointDistributions();
}

int Alias::pick() {
  int k = (int)(unif_rand()*(double)cells.size());
  return ( unif_rand() < cells[k].cutoff ? cells[k].index : cells[k].alias );
}

void Alias::computeTwoPointDistributions() {
    vector<Cell>::iterator litr = cells.end() - 1, ritr = cells.begin();
   
    // sort the vector of multiplied probabilities such that all cells
    // having probability less than 1 are at the end, cells with probabilities
    // greater than 1 are at the beginning
    // e.g., if initial cells are 0.8, 0.9, 1.1, 1.2
    // after this while loop, cells become 1.1, 1.2, 0.8, 0.9
    while (true) {
        while (ritr != cells.end() && ritr->cutoff >= 1.00000) {
            ritr++;
        }
        if (ritr == cells.end()) {
            break;
        }

        while (litr != cells.begin() && litr->cutoff <= 1.00000) {
            litr--;
        }

        if (litr == cells.begin() && litr->cutoff <= 1.00000) {
            break;
        }

        if (litr < ritr)
            break;

        Cell temp = *ritr;
        *ritr = *litr;
        *litr = temp;
        litr--; ritr++;
    }

    vector<Cell>::iterator li = cells.end() - 1, ri = cells.end() - 1;
    while (true) {
        //find a cell with probability initial cutoff greater than 1.0
        while (ri != cells.begin() && ri->cutoff <= 1.00000) {
            ri--;
        }
        if (ri == cells.begin() && ri->cutoff <= 1.00000) {
            break;
        }

        //find a cell with probability initial cutoff less than 1.0
        while (li != cells.begin() && li->cutoff >= 1.00000) {
            li--;
        }

        if (li == cells.begin() && li->cutoff >= 1.00000) {
            break;
        }

        li->alias = ri->index;
        ri->cutoff -= (1 - li->cutoff);
        li--;
    }
}

void Alias::printTwoPointDistributions() {
  REprintf("Two point distribution\n");
  REprintf("Cutoff:\t");
  for (int i = 0; i < cells.size(); i++) {
    Rprintf("%d\t", cells[i].cutoff);
  }
  REprintf("\n");
  REprintf("Index:\t");
  for (int i = 0; i < cells.size(); i++) {
    REprintf("%d\t", cells[i].index);
  }
  REprintf("\n");
  REprintf("Alias:\t");
  for (int i = 0; i < cells.size(); i++) {
    REprintf("%d\t", cells[i].alias);
  }
  REprintf("\n");
}


