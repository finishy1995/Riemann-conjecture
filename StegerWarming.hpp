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
     *  name    RiemannSolve   use AUSM to solve Riemann Problem
     *  return  void
     *  param   long   the index of array
     *  param   double**   the array(u & u1)
     *  param   double**   the target array(f)
     *
     **/
    void RiemannSolve(long, double**&, double**&);
    
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
    double ee = 0.00000001; // Set the minimum value to avoid lamda(+,-) = 0
    long i,j,k;
    double temp1,temp2;
    double epso[n],c[n];
    double lamda[3][n],lamdaP[3][n],lamdaN[3][n],fP[3][n],fN[3][n];
    double rv[n];
    for(i=0;i<n;i++) rv[n] = 0;
    
    for(i=0;i<m;i++)
    {
        //Computational feature splitting
        for(j=0;j<n;j++) epso[j] = ee;
        for(j=0;j<n;j++) c[j] = sqrt(gama*p[2][j]/p[1][j]);
        for(j=0;j<n;j++)
        {
            lamda[0][j] = p[0][j];
            lamda[1][j] = p[0][j]-c[j];
            lamda[2][j] = p[0][j]+c[j];
        }
        for(j=0;j<3;j++)
            for(k=0;k<n;k++)
            {
                lamdaP[j][k] = 0.5*(lamda[j][k]+sqrt(lamda[j][k]*lamda[j][k]+epso[k]*epso[k]));
                lamdaN[j][k] = 0.5*(lamda[j][k]-sqrt(lamda[j][k]*lamda[j][k]+epso[k]*epso[k]));
            }
        
        //Calculated positive and negative flux
        for(j=0;j<n;j++)
        {
            fP[0][j] = 0.5/gama*p[1][j]*(2.0*(gama-1)*lamdaP[0][j]+lamdaP[1][j]+lamdaP[2][j]);
            fP[1][j] = 0.5/gama*p[1][j]*(2.0*(gama-1)*lamdaP[0][j]*p[0][j]+lamdaP[1][j]*(p[0][j]-c[j])+lamdaP[2][j]*(p[0][j]+c[j]));
            fP[2][j] = 0.5/gama*p[1][j]*((gama-1)*lamdaP[0][j]*p[0][j]*p[0][j]+0.5*lamdaP[1][j]*(p[0][j]-c[j])*(p[0][j]-c[j])+0.5*lamdaP[2][j]*(p[0][j]+c[j])*(p[0][j]+c[j])+(0.5*(3-gama)/(gama-1)*(lamdaP[1][j]+lamdaP[2][j])*c[j]*c[j]));
            
            fN[0][j] = 0.5/gama*p[1][j]*(2.0*(gama-1)*lamdaN[0][j]+lamdaN[1][j]+lamdaN[2][j]);
            fN[1][j] = 0.5/gama*p[1][j]*(2.0*(gama-1)*lamdaN[0][j]*p[0][j]+lamdaN[1][j]*(p[0][j]-c[j])+lamdaN[2][j]*(p[0][j]+c[j]));
            fN[2][j] = 0.5/gama*p[1][j]*((gama-1)*lamdaN[0][j]*p[0][j]*p[0][j]+0.5*lamdaN[1][j]*(p[0][j]-c[j])*(p[0][j]-c[j])+0.5*lamdaN[2][j]*(p[0][j]+c[j])*(p[0][j]+c[j])+(0.5*(3-gama)/(gama-1)*(lamdaN[1][j]+lamdaN[2][j])*c[j]*c[j]));
        }
        
        //Calculation parameters
        for(j=1;j<(n-1);j++)
        {
            temp1 = (fP[0][j]-fP[0][j-1])/deltaX;
            temp2 = (fN[0][j+1]-fN[0][j])/deltaX;
            p[1][j] = p[1][j]-deltaT*(temp1+temp2);
            temp1 = (fP[1][j]-fP[1][j-1])/deltaX;
            temp2 = (fN[1][j+1]-fN[1][j])/deltaX;
            rv[j] = rv[j]-deltaT*(temp1+temp2);
            p[0][j] = rv[j]/p[1][j];
            temp1 = (fP[2][j]-fP[2][j-1])/deltaX;
            temp2 = (fN[2][j+1]-fN[2][j])/deltaX;
            e[j] = e[j]-deltaT*(temp1+temp2);
            p[2][j] = (gama-1)*(e[j]-0.5*p[1][j]*p[0][j]*p[0][j]);
        }
    }
}

#endif /* StegerWarming_hpp */
