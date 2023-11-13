///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// LTensor                                                                   //
//                                                                           //
// A Tensor Library with full Indicial Notation                              //
//                                                                           //
// Version 0.1                                                               //
// December 1, 2008                                                          //
//                                                                           //
// Copyright (C) 2007-2009-...                                               //
// Alejandro Limache & Sebastian Rojas Fredini                               //
//                                                                           //
// International Center of Computational Methods in Engineering  (CIMEC)     //
// Santa Fe, Argentina                                                       //
// alejandrolimache@gmail.com                                                //
//                                                                           //
// LTensor is freely available through the websites:                         //
// http://www.cimec.org.ar/alimache/                                         //
// http://code.google.com/p/ltensor/                                         //
// It may be copied, modified, and redistributed for non-commercial use as   // 
// long as the original authors and the Library get proper public credit     // 
// for its use.                                                              //
// Please consult the file LICENSE for the detailed copyright notices.       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#ifndef promote_H
#define promote_H


template<class T1, class T2>
class promote {
public:
  typedef T1 V;
};



template <class T>
class retPromote{
    public:
    typedef T V;

    };

template <class T>
class sqrtTrait{
public:
	static T function(T data){
		return sqrt((double)data);
	}
};




template <class T>
class absTrait{
public:
    static  T function(T data){
        return fabs((T)data);
        }

};

#define DECLARE_ABSFUNCTION(A,C) \
template<> class absTrait<A>{public: static A function(A data){return C((A)data);} }


#define DECLARE_SQRTFUNCTION(A,C) \
template<> class sqrtTrait<A>{public: static A function(A data){return C(data);} }

#define DECLARE_RETPROMOTE(A,C) \
template<> class retPromote<A> { public: typedef C V; }

#define DECLARE_PROMOTE(A,B,C) \
template<> class promote<A,B > { public: typedef C V; }

namespace std { template<typename T> class complex; }

DECLARE_PROMOTE(int,double,double);
DECLARE_PROMOTE(double,int,double);
DECLARE_PROMOTE(int,std::complex<double>,std::complex<double>);
DECLARE_PROMOTE(std::complex<double>,int,std::complex<double>);
DECLARE_PROMOTE(double,std::complex<double>,std::complex<double>);
DECLARE_PROMOTE(std::complex<double>,double,std::complex<double>);

DECLARE_RETPROMOTE(int,double);
DECLARE_RETPROMOTE(float,double);
DECLARE_RETPROMOTE(double,double);
DECLARE_RETPROMOTE(short,double);
DECLARE_RETPROMOTE(std::complex<double>,std::complex<double>);

DECLARE_SQRTFUNCTION(double,sqrt);
DECLARE_SQRTFUNCTION(float,sqrtf);

DECLARE_ABSFUNCTION(double,fabs);
DECLARE_ABSFUNCTION(int,abs);

#undef DECLARE_PROMOTE

#endif



