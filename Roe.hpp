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
    
    Roe(){}
    
    ~Roe(){}
    
    virtual void output(string) const;
    
    virtual void solve();
    
    /**
     *  name    setTol   set the contrast entropy correction 对比熵修正
     *  return  void
     *  param   double   the value
     *
     **/
    void setTol(double);
    
private:
    
    /**
     *  name    RiemannSolve   solve the Riemann matrix
     *  return  void
     *  param   long   the index numbers in the array
     *  param   double**   the array(u & u1)
     *  param   double**   the target array(f)
     *
     **/
    void RiemannSolve(long, double**&, double**&);
    
    /**
     *  Contrast entropy correction 对比熵修正
     */
    double tol = 0.000001;
};

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

void Roe::setTol(double value)
{
    if (value<0) tol = 0;
    else tol = value;
}

void Roe::solve()
{
    long count = 1;     // Main repeat times
    double tNow = 0;    // Time now (in main repeat)
    double CFL = 0.9;   // Normal CFL
    
    double** u1;        // Runge-Kutta uTemp
    newArray(u1);
    double** f;         // Convective numberical flux
    newArray(f);
    long i,j;
    // CFL_new actual CFL; max max number to get the delta time
    double CFL_new,deltaT,temp,max;
    
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
        if (count<=10) CFL_new = 0.2*CFL;
        else CFL_new = CFL;
        deltaT = CFL_new*deltaX/max;
        if (tNow+deltaT>t) deltaT=t-tNow;
        
        // Use there-order Runge-Kutta to solve the Riemann Problem
        for(i=0;i<(n-1);i++)
        {
            RiemannSolve(i, u, f);
            for(j=0;j<3;j++)
                if (i!=0) u1[j][i] = u[j][i]+deltaT/deltaX*(f[j][i-1]-f[j][i]);
        }
        for(i=0;i<(n-1);i++)
        {
            RiemannSolve(i, u1, f);
            for(j=0;j<3;j++)
                if (i!=0) u1[j][i] = 0.75*u[j][i]+0.25*(u1[j][i]+deltaT/deltaX*(f[j][i-1]-f[j][i]));
        }
        for(i=0;i<(n-1);i++)
        {
            RiemannSolve(i, u1, f);
            if (i!=0) for(j=0;j<3;j++) u[j][i] = 1.0/3.0*u[j][i]+2.0/3.0*(u1[j][i]+deltaT/deltaX*(f[j][i-1]-f[j][i]));
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
    
    deleteArray(f);
    deleteArray(u1);
}

void Roe::RiemannSolve(long index, double**& uTemp, double**& f)
{
    double l[4],r[4],s[3],lamda[3],sm[3],sFinal[3];
    double tempL,tempR;
    
    // Left and Right x position status
    // 0 Density; 1 Velocity; 2 Internal Energy; 3 Pressure
    l[0] = uTemp[0][index];
    l[1] = uTemp[1][index]/uTemp[0][index];
    l[2] = uTemp[2][index];
    l[3] = (l[2]-0.5*l[0]*l[1]*l[1])*(gama-1);
    r[0] = uTemp[0][index+1];
    r[1] = uTemp[1][index+1]/uTemp[0][index+1];
    r[2] = uTemp[2][index+1];
    r[3] = (r[2]-0.5*r[0]*r[1]*r[1])*(gama-1);
    
    // s[] use to get lamda matrix and get a parameter(h) to calculate f
    // lamda matrix:
    //      |  u   0   0  |     lamda[0] = u
    //      |  0  u+a  0  |     lamda[1] = u+a
    //      |  0   0  u-a |     lamda[2] = u-a
    // s[0] u; s[1] a; s[2] H
    s[0] = (sqrt(l[0])*l[1]+sqrt(r[0])*r[1])/(sqrt(l[0])+sqrt(r[0]));;
    tempL = (l[2]+l[3])/l[0];
    tempR = (r[2]+r[3])/r[0];
    s[2] = (sqrt(l[0])*tempL+sqrt(r[0])*tempR)/(sqrt(l[0])+sqrt(r[0]));
    s[1] = sqrt((gama-1)*(s[2]-0.5*s[0]*s[0]));
    
    // Get the lamda by contrast entropy correction
    if (fabs(s[0])>=tol) lamda[0] = fabs(s[0]);
    else lamda[0]=(s[0]*s[0]+tol*tol)/2.0/tol;
    if (fabs(s[0]+s[1])>=tol) lamda[1] = fabs(s[0]+s[1]);
    else lamda[1]=((s[0]+s[1])*(s[0]+s[1])+tol*tol)/2.0/tol;
    if (fabs(s[0]-s[1])>=tol) lamda[2] = fabs(s[0]-s[1]);
    else lamda[2]=((s[0]-s[1])*(s[0]-s[1])+tol*tol)/2.0/tol;
    
    
    /*d_u=UR-UL;
     Sm1=(gama-1)/a_^2*(d_u(1)*(H_-u_^2)+u_*d_u(2)-d_u(3));
     Sm2=1/2/a_*(d_u(1)*(-u_+a_)+d_u(2)-a_*Sm1);
     Sm3=d_u(1)-Sm1-Sm2;
     Sm=[Sm1;Sm2;Sm3];
     %solve flux at i+1/2
     FL=[rL*vL;rL*vL^2+pL;(EL+pL)*vL];
     FR=[rR*vR;rR*vR^2+pR;(ER+pR)*vR];
     out_flux=1/2*(FL+FR)-1/2*S*(abs_lamda.*Sm);*/
    // Matrix operation (expand and simplify)
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
    
    // Get convective numberical flux
    f[0][index] = 0.5*(l[0]*l[1]+r[0]*r[1]-sFinal[0]);
    f[1][index] = 0.5*(l[0]*l[1]*l[1]+l[3]+r[0]*r[1]*r[1]+r[3]-sFinal[1]);
    f[2][index] = 0.5*((l[2]+l[3])*l[1]+(r[2]+r[3])*r[1]-sFinal[2]);
}

#endif /* Roe_hpp */
