/**
 *  NumFormat.hpp
 *  abstract class NumFormat to solve Riemann Problem by num Format
 *
 *  Created by David Wang on 16/10/17.
 *
 */

#ifndef NumFormat_hpp
#define NumFormat_hpp

#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>

using namespace std;

class NumFormat
{
public:
    
    ~NumFormat();
    
    /**
     *  name    output   output results by .txt file
     *  param   string   the file name
     *
     **/
    virtual void output(string) const = 0;
    
    /**
     *  name    solve   solve the Riemann Problem
     *  return  void
     *
     **/
    virtual void solve() = 0;
    
    /**
     *  name    setT   set the end time
     *  return  void
     *  param   double   the value
     *
     **/
    void setT(double);
    
    /**
     *  name    setX1   set the left point of x axis
     *  return  void
     *  param   double   the value
     *
     **/
    void setX1(double);
    
    /**
     *  name    setX2   set the right point of y axis
     *  return  void
     *  param   double   the value
     *
     **/
    void setX2(double);
    
    /**
     *  name    setN   set the step number of x axis
     *  return  void
     *  param   long   the value
     *
     **/
    void setN(long);
    
    /**
     *  name    setMiddleX   set the cut-point in x axis
     *  return  void
     *  param   double   the value
     *
     **/
    void setMiddleX(double);
    
    /**
     *  name    setGama   set the gama γ
     *  return  void
     *  param   double   the value
     *
     **/
    void setGama(double);
    
    /**
     *  name    setFilePath
     *  return  void
     *  param   string   the file path
     *
     **/
    void setFilePath(string);
    
protected:
    
    NumFormat();
    
    /**
     *  name    newArray   new the target array double**
     *  return  void
     *  param   double**   target array
     *
     **/
    void newArray(double**&);
    
    /**
     *  name    deleteArray   delete the target array double**
     *  return  void
     *  param   double**   target array
     *
     **/
    void deleteArray(double**&);
    
    /**
     *  name    setDeltaX   set the step length of x axis
     *  return  void
     *
     **/
    void setDeltaX();
    
    /**
     *  name    reset   reset the answer array p
     *  return  void
     *
     **/
    void reset();
    
    /**
     *  name    setE   set the internal energy
     *  return  void
     *
     **/
    void setE();
    
    /**
     *  name    setU   set the transition matrix
     *  return  void
     *
     **/
    void setU();
    
    /**
     *  end time
     */
    double t = 0.2;
    
    /**
     *  left point of x axis
     */
    double x1 = 0;
    
    /**
     *  right point of x axis
     */
    double x2 = 1;
    
    /**
     *  step number of x axis
     */
    long n = 201;
    
    /**
     *  cut-point in x axis
     */
    double middleX = 0.5;
    
    /**
     *  step length of x axis
     */
    double deltaX;
    
    /**
     *  answer array
     */
    double** p = NULL;
    
    /**
     *  internal energy
     */
    double* e = NULL;
    
    /**
     *  transition matrix
     */
    double** u = NULL;
    
    /**
     *  gama γ
     */
    double gama = 1.4;
    
    /**
     *  file path
     */
    string filePath = "";
};


NumFormat::~NumFormat()
{
    deleteArray(p);
    delete [] e;
    e = NULL;
    deleteArray(u);
}

void NumFormat::setT(double value)
{
    if (t<0) t = 0;
    else t = value;
}

void NumFormat::setX1(double value)
{
    if (value>x2) x1 = x2;
    else x1 = value;
    setDeltaX();
}

void NumFormat::setX2(double value)
{
    if (value<x1) x2 = x1;
    else x2 = value;
    setDeltaX();
}

void NumFormat::setN(long value)
{
    if (n<2) n = 2;
    else n = value;
    setDeltaX();
    reset();
    setE();
    setU();
}

void NumFormat::setMiddleX(double value)
{
    if (value>x2) middleX = x2;
    else if (value<x1) middleX = x1;
    else middleX = value;
}


void NumFormat::setGama(double value)
{
    if (value<=0) gama = 1.4;
    else gama = value;
}

void NumFormat::setFilePath(string str)
{
    if (str != "") filePath = str;
}

NumFormat::NumFormat()
{
    setDeltaX();
    reset();
    setE();
    setU();
}

void NumFormat::newArray(double**& arr)
{
    int i;
    arr = new double*[3];
    for(i=0;i<3;i++) arr[i] = new double[n];
}

void NumFormat::deleteArray(double**& arr)
{
    if (arr!=NULL)
    {
        int i;
        for(i=0;i<3;i++) delete[] arr[i];
        delete[] arr;
        arr = NULL;
    }
}

void NumFormat::setDeltaX()
{
    deltaX = (x2-x1)/(n-1);
}

void NumFormat::reset()
{
    long i;
    double xpos;
    deleteArray(p);
    newArray(p);
    for(i=0;i<n;i++)
    {
        xpos = x1+(i-1)*deltaX;
        if (xpos<middleX)
        {
            p[0][i] = 0;
            p[1][i] = 1;
            p[2][i] = 1;
        } else {
            p[0][i] = 0;
            p[1][i] = 0.125;
            p[2][i] = 0.1;
        }
    }
}

void NumFormat::setE()
{
    long i;
    if (e!=NULL)
    {
        delete [] e;
    }
    e = new double[n];
    for(i=0;i<n;i++) e[i] = p[2][i]/(gama-1)+0.5*p[1][i]*p[0][i]*p[0][i];
}

void NumFormat::setU()
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

#endif /* NumFormat_hpp */
