#include "constraint.h"
#include <fstream>
#include <iostream>

namespace std;
/* note: only up to 3 dimensions are allowed */

constraint::constraint(int num_dim){
     construct_numVarMatrix_num_dim(num_dim);
}

void constraint::construct_numVarMatrix_num_dim(int num_dim){
    if (num_dim > 0 && num_dim <= 3){

    }

    else{
        cerr<<"invalid constraint dimensions given"<< endl;
        return;
    }

}