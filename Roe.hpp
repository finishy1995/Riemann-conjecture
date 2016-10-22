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

class Roe : public NumFormat
{
public:
    
    Roe();
    
    ~Roe();
    
    virtual void output(string) const;
    
    virtual void setN(long);
    
    /**
     *  name    setCFL   set the CFL
     *  return  void
     *  param   double   the value
     *
     **/
    void setCFL(double);
    
    /**
     *  name    setTol   set the TOL
     *  return  void
     *  param   double   the value
     *
     **/
    void setTol(double);
    
    /**
     *  name    solve   solve the Roe
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
     *  name    RiemannSolve   solve the Riemann matrix
     *  return  void
     *  param   long   the index numbers in the array
     *  param   double**   the array(u || u1)
     *  param   double**   the target array(f)
     *
     **/
    void RiemannSolve(long, double**&, double**&);
    
    /**
     *  CFL
     */
    double CFL = 0.9;
    
    /**
     *  time temporary value
     */
    double** u = NULL;
    
    /**
     *  Contrast entropy correction 对比熵修正
     */
    double tol = 0.000001;
};

Roe::Roe()
{
    construction();
    setU();
}

Roe::~Roe()
{
    deleteArray(u);
}

void Roe::output(string filename) const
{
    long i,j;
    ofstream out(filePath+filename);
    if (out.is_open())
    {
        out<<x1<<" "<<deltaX<<" "<<x2<<" "<<setprecision(12)<<fixed<<tol<<"\n";
        for(i=0;i<3;i++)
        {
            for(j=0;j<n;j++)
                out<<p[i][j]<<" ";
            out<<"\n";
        }
        out.close();
    }
    cout<<"Roe successful  tol="<<tol<<endl;
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

void Roe::setCFL(double value)
{
    if (value<0) CFL = 0;
    else CFL = value;
}

void Roe::setTol(double value)
{
    if (value<0) tol = 0;
    else tol = value;
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

void Roe::RiemannSolve(long index, double**& uTemp, double**& f)
{
    double lamda[3],l[6],r[6],s[3],tol,sm[3],sFinal[3];
    
    l[0] = uTemp[0][index];
    l[1] = uTemp[1][index]/uTemp[0][index];
    l[2] = uTemp[2][index];
    l[3] = (l[2]-0.5*l[0]*l[1]*l[1])*(gama-1);
    r[0] = uTemp[0][index+1];
    r[1] = uTemp[1][index+1]/uTemp[0][index+1];
    r[2] = uTemp[2][index+1];
    r[3] = (r[2]-0.5*r[0]*r[1]*r[1])*(gama-1);
    s[0] = (sqrt(l[0])*l[1]+sqrt(r[0])*r[1])/(sqrt(l[0])+sqrt(r[0]));;
    l[4] = l[3]/(gama-1)+0.5*l[0]*l[1]*l[1];
    r[4] = r[3]/(gama-1)+0.5*r[0]*r[1]*r[1];
    l[5] = (l[4]+l[3])/l[0];
    r[5] = (r[4]+r[3])/r[0];
    s[2] = (sqrt(l[0])*l[5]+sqrt(r[0])*r[5])/(sqrt(l[0])+sqrt(r[0]));
    s[1] = sqrt((gama-1)*(s[2]-0.5*s[1]*s[1]));
    
    //Get the eigenvalues
    if (fabs(s[0])>=tol) lamda[0] = fabs(s[0]);
    else lamda[0]=(s[0]*s[0]+tol*tol)/2.0/tol;
    if (fabs(s[0]+s[1])>=tol) lamda[1] = fabs(s[0]+s[1]);
    else lamda[1]=((s[0]+s[1])*(s[0]+s[1])+tol*tol)/2.0/tol;
    if (fabs(s[0]-s[1])>=tol) lamda[2] = fabs(s[0]-s[1]);
    else lamda[2]=((s[0]-s[1])*(s[0]-s[1])+tol*tol)/2.0/tol;
    
    double S[3][3] = {1,1,1,s[0],s[0]+s[1],s[0]-s[1],0.5*s[0]*s[0],s[2]+s[0]*s[1],s[2]-s[0]*s[1]};
    sm[0] = (gama-1)/(s[1]*s[1])*((r[0]-l[0])*(s[2]-s[0]*s[0])+s[0]*(r[1]*r[0]-l[1]*l[0])-(r[2]-l[2]));
    sm[1]=0.5/s[1]*((r[0]-l[0])*(-s[0]+s[1])+(r[1]*r[0]-l[1]*l[0])-s[1]*sm[0]);
    sm[2]=(r[0]-l[0])-sm[0]-sm[1];
    lamda[0] = lamda[0]*sm[0];
    lamda[1] = lamda[1]*sm[1];
    lamda[2] = lamda[2]*sm[2];
    sFinal[0] = S[0][0]*lamda[0]+S[0][1]*lamda[1]+S[0][2]*lamda[2];
    sFinal[1] = S[1][0]*lamda[0]+S[1][1]*lamda[1]+S[1][2]*lamda[2];
    sFinal[2] = S[2][0]*lamda[0]+S[2][1]*lamda[1]+S[2][2]*lamda[2];
    f[0][index] = 0.5*(l[0]*l[1]+r[0]*r[1]-sFinal[0]);
    f[1][index] = 0.5*(l[0]*l[1]*l[1]+l[3]+r[0]*r[1]*r[1]+r[3]-sFinal[1]);
    f[2][index] = 0.5*((l[2]+l[3])*l[1]+(r[2]+r[3])*r[1]-sFinal[2]);
}

void Roe::solve()
{
    long i,j;
    double** u1;
    newArray(u1);
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
            if ((fabs(p[0][i])+a[i])>max) max = fabs(p[0][i])+a[i];
        if (count<=10) CFL_new = 0.2*CFL;
        else CFL_new = CFL;
        deltaT = CFL_new*deltaX/max;
        if (tNow+deltaT>t) deltaT=t-tNow;
    
        // there-order Runge-Kutta
        RiemannSolve(0, u, f);
        for(i=1;i<(n-1);i++)
        {
            RiemannSolve(i, u, f);
            for(j=0;j<3;j++) u1[j][i] = u[j][i]+deltaT/deltaX*(f[j][i-1]-f[j][i]);
        }
        RiemannSolve(0, u1, f);
        for(i=1;i<(n-1);i++)
        {
            RiemannSolve(i, u1, f);
            for(j=0;j<3;j++) u1[j][i] = 0.75*u[j][i]+0.25*(u1[j][i]+deltaT/deltaX*(f[j][i-1]-f[j][i]));
        }
        RiemannSolve(0, u1, f);
        for(i=1;i<(n-1);i++)
        {
            RiemannSolve(i, u1, f);
            for(j=0;j<3;j++) u[j][i] = 1.0/3.0*u[j][i]+2.0/3.0*(u1[j][i]+deltaT/deltaX*(f[j][i-1]-f[j][i]));
            p[0][i] = u[1][i]/u[0][i];
            p[1][i] = u[0][i];
            e[i] = u[2][i];
            p[2][i] = (e[i]-1/2*p[1][i]*p[0][i]*p[0][i])*(gama-1);
        }
        tNow = tNow+deltaT;
        count++;
    }
    
    deleteArray(f);
    delete [] a;
    deleteArray(u1);
}

#endif /* Roe_hpp */
