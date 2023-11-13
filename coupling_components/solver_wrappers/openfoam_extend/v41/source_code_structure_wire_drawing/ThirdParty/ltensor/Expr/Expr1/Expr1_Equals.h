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

#ifndef Expr1_Equals_H
#define Expr1_Equals_H

//#define CHECK_Expr1_Definitions
//#define USE_ASSERT_Expr1
#ifdef USE_ASSERT_Expr1
#include <assert.h>
#endif


//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as T(i) = C(i)
//////////////////////////////////////////////////////////////////////


template<class A, class T, char i>
template<class B, class U>
inline Expr1<A,T,i> &
Expr1<A,T,i>::operator=(const Expr1<B,U,i> &rhs)
{
	Li_equals_Ri((*this),rhs);
	return *this;
}




template<class A, class T, char i>
inline Expr1<A,T,i> &
Expr1<A,T,i>::operator=(const Expr1<A,T,i> &rhs)
{
	Li_equals_Ri((*this),rhs);
	return *this;
}




template<class A, class T, char i>
template<class B,class U>
inline Expr1<A,T,i> &
Expr1<A,T,i>::operator+=(const Expr1<B,U,i> &rhs)
{
	Li_plusequals_Ri((*this),rhs);
	return *this;
}

template<class A, class T, char i>
template<class B,class U>
inline Expr1<A,T,i> &
Expr1<A,T,i>::operator-=(const Expr1<B,U,i> &rhs)
{
	Li_minusequals_Ri((*this),rhs);
	return *this;
}

#endif
