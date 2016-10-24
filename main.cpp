/**
 *  main.cpp
 *  the main procedure to use c++ get the data and use matlab to plot the image
 *  all the number format using three-order Runge-Kutta
 *
 *  Created by David Wang on 16/10/17.
 *
 */

#include "Roe.hpp"
#include "StegerWarming.hpp"
#include "AUSM.hpp"

// The output full file path. Please change to your code path before running.
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
    AUSM ausm;
    ausm.setFilePath(FILEPATH);
    ausm.solve();
    ausm.output("AUSM.txt");
    
    return 0;
}
