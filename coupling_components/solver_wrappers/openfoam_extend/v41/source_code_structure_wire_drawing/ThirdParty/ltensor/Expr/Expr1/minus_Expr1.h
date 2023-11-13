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
/* Declares a wrapper class for the unary minus (-) operator. */
#ifndef minus_Expr1_H
#define minus_Expr1_H


//#define CHECK_minus_Expr1
//#define USE_ASSERT_Expr1
#ifdef USE_ASSERT_Expr1
#include <assert.h>
#endif

#define ExprAi   Expr1<A,T,i>

template < class A, class T, char i>
inline const Expr1 < const minus_Ai < ExprAi, T>, T, i>
operator- (const ExprAi  &a)
{
	typedef const minus_Ai < ExprAi, T> ExprObj;
	return Expr1 < ExprObj, T, i > ( ExprObj(a) );
}

#undef ExprAi


#endif
