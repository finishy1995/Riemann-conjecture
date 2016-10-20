/**
 *  main.cpp
 *  the main procedure to use c++ get the data and use matlab to plot the image
 *
 *  Created by David Wang on 16/10/17.
 *
 */

#include <iostream>
#include "Roe.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    Roe test;
    cout<<test.output()[1][101]<<endl;
    
    return 0;
}
