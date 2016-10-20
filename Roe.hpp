/**
 *  Roe.hpp
 *  Roe class to solve Riemann Problem by using Roe Format
 *
 *  Created by David Wang on 16/10/17.
 *
 */

#ifndef Roe_hpp
#define Roe_hpp

#include "NumFormat.hpp"
#include <math.h>

class Roe : public NumFormat
{
public:
    
    Roe();
    
    ~Roe();
    
    virtual double** output() const;
    
    virtual void setN(long);
    
    /**
     *  name    setGama   set the gama
     *  return  void
     *  param   double   the value
     *
     **/
    void setGama(double);
    
    /**
     *  name    setCFL   set the CFL
     *  return  void
     *  param   double   the value
     *
     **/
    void setCFL(double);
    
private:
    
    /**
     *  name    setE   set the internal energy
     *  return  void
     *
     **/
    void setE();
    
    /**
     *  name    setU   set the time temporary value
     *  return  void
     *
     **/
    void setU();
    
    /**
     *  name    solve   solve the Riemann
     *  return  void
     *
     **/
    void solve();
    
    /**
     *  gama
     */
    double gama = 1.4;
    
    /**
     *  CFL
     */
    double CFL = 0.9;
    
    /**
     *  internal energy
     */
    double* e = NULL;
    
    /**
     *  time temporary value
     */
    double** u = NULL;
};

Roe::Roe()
{
    construction();
    setE();
    setU();
}

Roe::~Roe()
{
    deleteArray(u);
    delete [] e;
    e = NULL;
}

double** Roe::output() const
{
    return p;
}

void Roe::setN(long value)
{
    if (n<2) n = 2;
    else n = value;
    setDeltaX();
    reset();
    setE();
    setU();
}

void Roe::setGama(double value)
{
    if (value<0) gama = 0;
    else gama = value;
}

void Roe::setCFL(double value)
{
    if (value<0) CFL = 0;
    else CFL = value;
}

void Roe::setE()
{
    long i;
    if (e!=NULL)
    {
        delete [] e;
    }
    e = new double[n];
    for(i=0;i<n;i++) e[i] = p[2][i]/(gama-1)+0.5*p[1][i]*p[0][i]*p[0][i];
}

void Roe::setU()
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

void Roe::solve()
{
    long i,j;
    double** u1;
    double* a = new double[n];
    double** f;
    newArray(f);
    long count = 1;
    double tNow = 0;
    double max = 0;
    double CFL_new,deltaT;
    
    for(i=0;i<3;i++)
        for(j=0;j<n;j++)
            u1[i][j] = u[i][j];
    while (tNow<t) {
        for(i=0;i<n;i++) a[i]=sqrt(gama*p[2][i]/p[1][i]);
        for(i=0;i<n;i++)
            if ((fabs(p[0][i])+a[i])>max) max = (fabs(p[0][i])+a[i])>max;
        if (count<=10) CFL_new = 0.2*CFL;
        else CFL_new = CFL;
        deltaT = CFL_new*deltaX/max;
        if (tNow+deltaT>t) deltaT=t-tNow;
    }
    /*
     %Riemann Reoblem:compute flux at 1+1/2
     F(:,1)=Riemann_Roe(U(:,1),U(:,2));
     
     %compute U in next step
     %3阶龙阁库塔方法
     for i=2:N-1
     %Riemann Reoblem:compute flux at i+1/2
     F(:,i)=Riemann_Roe(U(:,i),U(:,i+1));
     U1(:,i)=U(:,i)+d_t/d_x*(F(:,i-1)-F(:,i));
     end
     F(:,1)=Riemann_Roe(U1(:,1),U1(:,2));
     for i=2:N-1
     %Riemann Reoblem:compute flux at i+1/2
     F(:,i)=Riemann_Roe(U1(:,i),U1(:,i+1));
     U1(:,i)=U(:,i)*3/4+1/4*(U1(:,i)+d_t/d_x*(F(:,i-1)-F(:,i)));
     end
     F(:,1)=Riemann_Roe(U1(:,1),U1(:,2));
     for i=2:N-1
     %Riemann Reoblem:compute flux at i+1/2
     F(:,i)=Riemann_Roe(U1(:,i),U1(:,i+1));
     U(:,i)=U(:,i)*1/3+2/3*(U1(:,i)+d_t/d_x*(F(:,i-1)-F(:,i)));
     r(i)=U(1,i);
     v(i)=U(2,i)/U(1,i);
     E(i)=U(3,i);
     p(i)=(E(i)-1/2*r(i)*v(i)^2)*(gama-1);
     end
     Time=Time+d_t
     icount=icount+1;
     end  */
    
    deleteArray(f);
    delete [] a;
    deleteArray(u1);
}

#endif /* Roe_hpp */
