/**
 *  main.cpp
 *  the main procedure to use c++ get the data and use matlab to plot the image
 *
 *  Created by David Wang on 16/10/17.
 *
 */

#include <iostream>
#include "Roe.hpp"

const string FILEPATH="/Users/cncuser/Desktop/c_code/CFD2/CFD2/";

int main(int argc, const char * argv[]) {
    Roe roe;
    roe.setFilePath(FILEPATH);
    roe.solve();
    roe.output("Roe.txt");
    //long i;
    //for(i=0;i<201;i++) cout<<result[0][i]<<" ";
    //cout<<endl;
    //cout<<result[0][37]<<endl;
    
    return 0;
}
