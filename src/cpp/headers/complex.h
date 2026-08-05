#ifndef COMPLEX_H
#define COMPLEX_H

class complex
{
public:

   double real;
   double imag;

   complex();
   complex(const complex& y);
   complex(const double& r, const double& i);
   ~complex();

   complex& operator = (const complex& y);
   complex& operator = (const double& y);
};

// other functions defined as inlines

complex operator - (const complex& x);
complex conj(const complex& x);
complex operator + (const complex& x, const complex& y);
complex operator + (const complex& x, const double& y);
complex operator + (const double& x, const complex& y);
complex operator - (const complex& x, const complex& y);
complex operator - (const complex& x, const double& y);
complex operator - (const double& x, const complex& y);
complex operator * (const complex& x, const complex& y);
complex operator * (const complex& x, const double& y);
complex operator * (const double& x, const complex& y);
complex operator / (const complex& x, const complex& y);
complex operator / (const complex& x, const double& r);
complex operator / (const double& r, const complex& y);
complex cos(const complex& x);
complex sin(const complex& x);
double abs(const complex& x);

inline complex::complex ()
{
   real = 0.;
   imag = 0.;
}

inline complex::complex (const complex& y)
{
    real = y.real;
    imag = y.imag;
}

inline complex::complex(const double& r, const double& i){
    real = r;
    imag = i;
}

complex::~complex()
{
}

inline complex& complex::operator = (const complex& y) 
{ 
    real = y.real; imag = y.imag; return *this; 
} 

inline complex& complex::operator = (const double& y) 
{ 
    real = y; imag = 0; return *this; 
} 

inline complex operator - (const complex& x)
{
    return complex(-x.real, -x.imag);
}

inline complex conj(const complex& x)
{
    return complex(x.real, -x.imag);
}

inline complex operator + (const complex& x, const complex& y)
{
    return complex(x.real + y.real, x.imag + y.imag);
}

inline complex operator + (const complex& x, const double& y)
{
    return complex(x.real + y, x.imag);
}

inline complex operator + (const double& x, const complex& y)
{
    return complex(x + y.real, y.imag);
}

inline complex operator - (const complex& x, const complex& y)
{
    return complex(x.real - y.real, x.imag - y.imag);
}

inline complex operator - (const complex& x, const double& y)
{
    return complex(x.real - y, x.imag);
}

inline complex operator - (const double& x, const complex& y)
{
    return complex(x - y.real, -y.imag);
}

inline complex operator * (const complex& x, const complex& y)
{
    return complex(x.real * y.real - x.imag * y.imag, 
                   x.real * y.imag + x.imag * y.real);
}

inline complex operator * (const complex& x, const double& y)
{
    return complex(x.real * y, x.imag * y);
}

inline complex operator * (const double& x, const complex& y)
{
    return complex(x * y.real, x * y.imag);
}

inline complex operator / (const complex& x, const complex& y)
{
    double d=y.real*y.real+y.imag*y.imag;
    return(complex((x.real*y.real+x.imag*y.imag)/d,
                   (x.imag*y.real-x.real*y.imag)/d));
}

inline complex operator / (const double& r, const complex& y)
{
    double d=y.real*y.real+y.imag*y.imag;
    return(complex((r*y.real)/d,(-r*y.imag)/d));
}

inline complex operator / (const complex& x, const double& r)
{
    return(complex(x.real/r,x.imag/r));
}

inline complex cos(const complex& x)
{
    return(complex(cos(x.real)*cosh(x.imag),-sin(x.real)*sinh(x.imag)));
}

inline complex sin(const complex& x)
{
    return(complex(sin(x.real)*cosh(x.imag),cos(x.real)*sinh(x.imag)));
}

inline double abs(const complex& x)
{
    return sqrt(x.real*x.real+x.imag*x.imag);
}

#endif
