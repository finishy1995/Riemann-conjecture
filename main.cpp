/**
 *  main.cpp
 *  the main procedure to use c++ get the data and use matlab to plot the image
 *
 *  Created by David Wang on 16/10/17.
 *
 */

#include <iostream>
#include "Roe.hpp"
#include "StegerWarming.hpp"

const string FILEPATH="/Users/cncuser/Desktop/c_code/CFD2/CFD2/";

int main(int argc, const char * argv[]) {
    Roe roe;
    roe.setFilePath(FILEPATH);
    roe.solve();        //normal tol=0.000001
    roe.output("Roe1.txt");
    roe.setTol(0.0001);
    roe.solve();
    roe.output("Roe2.txt");
    roe.setTol(0.00000001);
    roe.solve();
    roe.output("Roe3.txt");
    roe.setTol(1);
    roe.solve();
    roe.output("Roe4.txt");
    roe.setTol(0);
    roe.solve();
    roe.output("Roe5.txt");
    
    StegerWarming sw;
    sw.setFilePath(FILEPATH);
    sw.solve();
    sw.output("StegerWarming.txt");
    
    return 0;
}
