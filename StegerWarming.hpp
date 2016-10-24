/**
 *  StegerWarming.hpp
 *  StegerWarming class to solve Riemann Problem by using StegerWarming Format
 *
 *  Created by David Wang on 16/10/22.
 *
 */

#ifndef StegerWarming_hpp
#define StegerWarming_hpp

#include "NumFormat.hpp"

class StegerWarming : public NumFormat
{
public:
    
    StegerWarming(){setM();}
    
    ~StegerWarming(){}
    
    virtual void output(string) const;
    
    virtual void solve();
    
    /**
     *  name    setDeltaT   set the step length of t axis
     *  return  void
     *  param   double  the value
     *
     **/
    void setDeltaT(double);
    
private:
    
    /**
     *  name    setM   set the step number of t axis
     *  return  void
     *
     **/
    void setM();
    
    /**
     *  name    RiemannSolve   use StegerWarming to solve Riemann Problem
     *  return  void
     *  param   long   the index of array
     *  param   double**   the array(u & u1)
     *  param   double**   the target array(fP)
     *  param   double**   the target array(fN)
     *
     **/
    void RiemannSolve(long, double**&, double**&, double**&);
    
    /**
     *  step length of t axis
     */
    double deltaT = 0.0001;
    
    long m;
};


void StegerWarming::output(string filename) const
{
    long i,j;
    ofstream out(filePath+filename);
    if (out.is_open())
    {
        out<<x1<<" "<<deltaX<<" "<<x2<<" "<<setprecision(12)<<fixed<<"\n";
        for(i=0;i<3;i++)
        {
            for(j=0;j<n;j++)
                out<<p[i][j]<<" ";
            out<<"\n";
        }
        out.close();
    }
    cout<<"StegerWarming successful"<<endl;
}

void StegerWarming::setDeltaT(double value)
{
    if (value<=0) deltaT = 0.0001;
    else if (value>t) deltaT = t;
    else deltaT = value;
    setM();
}

void StegerWarming::setM()
{
    m = ceil(t/deltaT);
}

void StegerWarming::solve()
{
    double** u1;            // Runge-Kutta uTemp
    newArray(u1);
    double** fP;            // Positive convective numberical flux
    newArray(fP);
    double** fN;            // Negative convective numberical flux
    newArray(fN);
    long i,j,k;
    
    for(i=0;i<3;i++)
        for(j=0;j<n;j++)
            u1[i][j] = u[i][j];
    
    for(k=0;k<m;k++)
    {
        // Use there-order Runge-Kutta to solve the Riemann Problem
        for(i=0;i<n;i++) RiemannSolve(i, u, fP, fN);
        for(i=1;i<(n-1);i++)
            for(j=0;j<3;j++)
                u1[j][i] = u[j][i]+deltaT/deltaX*(fP[j][i-1]-fP[j][i]+fN[j][i]-fN[j][i+1]);
        for(i=0;i<n;i++) RiemannSolve(i, u1, fP, fN);
        for(i=1;i<(n-1);i++)
            for(j=0;j<3;j++)
                u1[j][i] = 0.75*u[j][i]+0.25*(u1[j][i]+deltaT/deltaX*(fP[j][i-1]-fP[j][i]+fN[j][i]-fN[j][i+1]));
        for(i=0;i<n;i++) RiemannSolve(i, u1, fP, fN);
        for(i=1;i<(n-1);i++)
            for(j=0;j<3;j++)
                u[j][i] = 1.0/3.0*u[j][i]+2.0/3.0*(u1[j][i]+deltaT/deltaX*(fP[j][i-1]-fP[j][i]+fN[j][i]-fN[j][i+1]));
    }
    
    // Get the value of final Velocity(p[0]), Density(p[1]), Pressure(p[2])
    // and Internal Energy(e)
    for(i=0;i<n;i++)
    {
        p[1][i] = u[0][i];
        p[0][i] = u[1][i]/u[0][i];
        e[i] = u[2][i];
        p[2][i] = (e[i]-1/2*p[1][i]*p[0][i]*p[0][i])*(gama-1);
    }
    
    deleteArray(u1);
    deleteArray(fP);
    deleteArray(fN);
}

void StegerWarming::RiemannSolve(long index, double**& uTemp, double**& fP, double**& fN)
{
    // Set the minimum value to avoid lamda(+,-) = 0
    double ee = 0.00000001;
    
    double lamda[3],lamdaP[3],lamdaN[3],uTran[3];
    long i;
    double c;

    // Translation uTemp to Velocity, Density and Internal Energy
    uTran[0] = uTemp[1][index]/uTemp[0][index];
    uTran[1] = uTemp[0][index];
    uTran[2] = uTemp[2][index];
    
    // Computational feature splitting
    c = sqrt(gama*((gama-1.0)*(uTran[2]-0.5*uTran[0]*uTran[0]*uTran[1]))/uTran[1]);
    lamda[0] = uTran[0];
    lamda[1] = uTran[0]-c;
    lamda[2] = uTran[0]+c;
    for(i=0;i<3;i++)
    {
        lamdaP[i] = 0.5*(lamda[i]+sqrt(lamda[i]*lamda[i]+ee*ee));
        lamdaN[i] = 0.5*(lamda[i]-sqrt(lamda[i]*lamda[i]+ee*ee));
    }
    
    // Calculated positive and negative flux
    fP[0][index] = 0.5/gama*uTran[1]*(2.0*(gama-1)*lamdaP[0]+lamdaP[1]+lamdaP[2]);
    fP[1][index] = 0.5/gama*uTran[1]*(2.0*(gama-1)*lamdaP[0]*uTran[0]+lamdaP[1]*(uTran[0]-c)+lamdaP[2]*(uTran[0]+c));
    fP[2][index] = 0.5/gama*uTran[1]*((gama-1)*lamdaP[0]*uTran[0]*uTran[0]+0.5*lamdaP[1]*(uTran[0]-c)*(uTran[0]-c)+0.5*lamdaP[2]*(uTran[0]+c)*(uTran[0]+c)+(0.5*(3-gama)/(gama-1)*(lamdaP[1]+lamdaP[2])*c*c));
    
    // Get convective numberical flux
    fN[0][index] = 0.5/gama*uTran[1]*(2.0*(gama-1)*lamdaN[0]+lamdaN[1]+lamdaN[2]);
    fN[1][index] = 0.5/gama*uTran[1]*(2.0*(gama-1)*lamdaN[0]*uTran[0]+lamdaN[1]*(uTran[0]-c)+lamdaN[2]*(uTran[0]+c));
    fN[2][index] = 0.5/gama*uTran[1]*((gama-1)*lamdaN[0]*uTran[0]*uTran[0]+0.5*lamdaN[1]*(uTran[0]-c)*(uTran[0]-c)+0.5*lamdaN[2]*(uTran[0]+c)*(uTran[0]+c)+(0.5*(3-gama)/(gama-1)*(lamdaN[1]+lamdaN[2])*c*c));
}

#endif /* StegerWarming_hpp */
