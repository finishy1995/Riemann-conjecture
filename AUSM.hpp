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
    
    AUSM();
    
    ~AUSM();
    
    virtual void output(string) const;
    
    virtual void setN(long);
    
    /**
     *  name    solve   solve the Riemann by StegerWarming
     *  return  void
     *
     **/
    virtual void solve();
    
private:
    
    /**
     *  name    setU   set the time temporary value
     *  return  void
     *
     **/
    void setU();
    
    /**
     *  time temporary value
     */
    double** u = NULL;
};

AUSM::AUSM()
{
    construction();
    newArray(u);
}

AUSM::~AUSM()
{
    deleteArray(u);
}

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

void AUSM::setN(long value)
{
    if (n<2) n = 2;
    else n = value;
    setDeltaX();
    reset();
    setE();
}

void AUSM::solve()
{
    long i,j;
    double max = 1e-100;
    double temp,deltaT,tNow;
    double** f;
    double** nx;
    newArray(f);
    newArray(nx);
    
    for(i=0;i<n;i++)
    {
        temp = sqrt(gama*((gama-1)*(u[2][i]-0.5*u[1][i]*u[1][i]/u[0][i]))/u[0][i])+fabs(u[1][i]/u[0][i]);
        if (temp>max) max = temp;
    }
    deltaT = 0.4*deltaX/max;
    tNow = 0;
    
    while (tNow<t)
    {
        tNow = tNow+deltaT;
        for(i=0;i<(n-1);i++)
        {
            RiemannSolve(i, f);
            for(j=0;j<3;j++)
                if (i!=0) nx[j][i] = 0.5*0.2*(u[j][i+1]-2.0*u[j][i]+u[j][i-1]);
        }
        for(i=1;i<n;i++)
        {
            for(j=0;j<3;j++)
                u[j][i]=u[j][i]-deltaT/deltaX*(f[j][i]-f[j][i-1])+nx[j][i];
        }
    }
    for(i=0;i<n;i++)
    {
        p[1][i] = u[0][i];
        p[0][i] = u[1][i]/u[0][i];
        e[i] = u[2][i];
        p[2][i] = (e[i]-1/2*p[1][i]*p[0][i]*p[0][i])*(gama-1);
    }
    
    deleteArray(nx);
    deleteArray(f);
    
    /*//RiemannSolve//call AUSM(Q(i,1), Q(i,2), Q(i,3), Q(i+1,1), Q(i+1,2), Q(i+1,3), F(i,1), F(i,2), F(i,3))
     */
}

void AUSM::setU()
{
    long i;
    deleteArray(u);
    newArray(u);
    for(i=0;i<n;i++)
    {
        u[0][i] = p[1][i];
        u[1][i] = p[0][i]*p[1][i];
        u[2][i] = e[i];
    }
}

#endif /* AUSM_hpp */
