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
/* Declares a wrapper class for the unary plus (+) operator. */
#ifndef plus_Expr2_H
#define plus_Expr2_H


//#define CHECK_plus_Expr2

//#define USE_ASSERT_Expr2
#ifdef USE_ASSERT_Expr2
#include <assert.h>
#endif

#define ExprAij   Expr2<A,T,i,j>

template < class A, class T, char i , char j>
inline const Expr2 < const plus_Aij < ExprAij, T>, T, i, j >
operator+ (const ExprAij  &a)
{
	typedef const plus_Aij < ExprAij, T> ExprObj;
	return Expr2 < ExprObj, T, i, j > ( ExprObj(a) );
}

#undef ExprAij

#endif
