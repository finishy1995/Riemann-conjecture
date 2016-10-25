/**
 *  AUSM.hpp
 *  AUSM class to solve Riemann Problem by using AUSM Format
 *
 *  Created by David Wang on 16/10/23.
 *
 */

#ifndef AUSM_hpp
#define AUSM_hpp

#include "NumFormat.hpp"

class AUSM : public NumFormat
{
public:
    
    AUSM(){}
    
    ~AUSM(){}
    
    virtual void output(string) const;
    
    virtual void solve();
    
private:
    
    /**
     *  name    RiemannSolve   use AUSM to solve Riemann Problem
     *  return  void
     *  param   long   the index of array
     *  param   double**   the array(u & u1)
     *  param   double**   the target array(f)
     *
     **/
    void RiemannSolve(long, double**&, double**&);
};


void AUSM::output(string filename) const
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
    cout<<"AUSM successful"<<endl;
}

void AUSM::solve()
{
    long count = 1;     // Main repeat times
    double tNow = 0;    // Time now (in main repeat)
    double CFL = 0.4;   // Normal CFL
    
    double** u1;        // Runge-Kutta uTemp
    newArray(u1);
    double** f;         // Convective numberical flux
    newArray(f);
    double** nx;        // Viscous terms
    newArray(nx);
    long i,j;
    // CFL_new actual CFL; max max number to get the delta time
    double temp,deltaT,CFL_new,max;
    
    reset();
    setE();
    setU();
    for(i=0;i<3;i++)
        for(j=0;j<n;j++)
            u1[i][j] = u[i][j];
    
    while (tNow<t)
    {
        // Get the value of delta t
        max = 1e-99;
        for(i=0;i<n;i++)
        {
            temp = sqrt(gama*((gama-1.0)*(u[2][i]-0.5*u[1][i]*u[1][i]/u[0][i]))/u[0][i])+fabs(u[1][i]/u[0][i]);
            if (temp>max) max = temp;
        }
        if (count<=10) CFL_new = 0.5*CFL;
        else CFL_new = CFL;
        deltaT = CFL_new*deltaX/max;
        if (tNow+deltaT>t) deltaT=t-tNow;
        
        // Use there-order Runge-Kutta to solve the Riemann Problem
        for(i=0;i<(n-1);i++)
        {
            RiemannSolve(i, u, f);
            for(j=0;j<3;j++)
            {
                if (i!=0) nx[j][i] = 0.5*0.2*(u[j][i+1]-2.0*u[j][i]+u[j][i-1]);
                if (i!=0) u1[j][i] = u[j][i]+deltaT/deltaX*(f[j][i-1]-f[j][i])+nx[j][i];
            }
        }
        for(i=0;i<(n-1);i++)
        {
            RiemannSolve(i, u1, f);
            for(j=0;j<3;j++)
                if (i!=0) u1[j][i] = 0.75*u[j][i]+0.25*(u1[j][i]+deltaT/deltaX*(f[j][i-1]-f[j][i])+nx[j][i]);
        }
        for(i=0;i<(n-1);i++)
        {
            RiemannSolve(i, u1, f);
            for(j=0;j<3;j++)
                if (i!=0) u[j][i] = 1.0/3.0*u[j][i]+2.0/3.0*(u1[j][i]+deltaT/deltaX*(f[j][i-1]-f[j][i])+nx[j][i]);
        }
        
        count++;
        tNow = tNow+deltaT;
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
    deleteArray(nx);
    deleteArray(f);
}

void AUSM::RiemannSolve(long index, double**& uTemp, double**& f)
{
    double l[4],r[4];
    double mPos,pPos,mNeg,pNeg,mMid,pMid;
    
    // Left and Right x position status
    // 0 Velocity; 1 Pressure; 2 Intermediate Parameter; 3 Mach Number
    l[0] = u[1][index]/u[0][index];
    l[1] = (gama-1)*(u[2][index]-0.5*u[1][index]*u[1][index]/u[0][index]);
    l[2] = sqrt(gama*l[1]/u[0][index]);
    l[3] = l[0]/l[2];
    r[0] = u[1][index+1]/u[0][index+1];
    r[1] = (gama-1)*(u[2][index+1]-0.5*u[1][index+1]*u[1][index+1]/u[0][index+1]);
    r[2] = sqrt(gama*r[1]/u[0][index+1]);
    r[3] = r[0]/r[2];
    
    // Get actual mach number and pressure
    if (l[3]>=1)
    {
        mPos = l[3];
        pPos = l[1];
    } else if (fabs(l[3])<1) {
        mPos = 0.25*(l[3]+1.0)*(l[3]+1.0);
        pPos = 0.25*l[1]*(l[3]+1.0)*(l[3]+1.0)*(2.0-l[3]);
    } else {
        mPos = 0;
        pPos = 0;
    }
    if (r[3]>=1)
    {
        mNeg = 0;
        pNeg = 0;
    } else if (fabs(r[3])<1) {
        mNeg = -0.25*(r[3]-1.0)*(r[3]-1.0);
        pNeg = 0.25*r[1]*(r[3]-1.0)*(r[3]-1.0)*(2.0+r[3]);
    } else {
        mNeg = r[3];
        pNeg = r[1];
    }
    mMid = mPos+mNeg;
    pMid = pPos+pNeg;
    
    // Get convective numberical flux
    f[0][index] = 0.5*mMid*(u[0][index]*l[2]+u[0][index+1]*r[2])-0.5*fabs(mMid)*(u[0][index+1]*r[2]-u[0][index]*l[2]);
    f[1][index] = 0.5*mMid*(u[1][index]*l[2]+u[1][index+1]*r[2])-0.5*fabs(mMid)*(u[1][index+1]*r[2]-u[1][index]*l[2])+pMid;
    f[2][index] = 0.5*mMid*((u[2][index]+l[1])*l[2]+(u[2][index+1]+r[1])*r[2])-0.5*fabs(mMid)*((u[2][index+1]+r[1])*r[2]-(u[2][index]+l[1])*l[2]);
}

#endif /* AUSM_hpp */
