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

/* Adds two Tensor1's together, yielding a Tensor1. */
#ifndef Expr1_minus_Expr1_H
#define Expr1_minus_Expr1_H


//#define CHECK_Expr1_minus_Expr1
//#define USE_ASSERT_Expr1
#ifdef USE_ASSERT_Expr1
#include <assert.h>
#endif

#include "./Expr1.h"

#define ExprAi  Expr1<A,T,i>
#define ExprBi  Expr1<B,U,i>
#define PromotedType promote<T,U>::V


/*!\brief Minus operation between two Expr Objects of rank 1
*	
*	Performs the - operation and returns another Expr of rank 1 containing the result
*/

template < class A, class B, class T,class U, char i >
inline const Expr1 <  const Ai_minus_Bi < ExprAi,ExprBi , typename PromotedType >, typename PromotedType, i >
operator-  (const ExprAi &ExprL, const ExprBi &ExprR)
{
	typedef const Ai_minus_Bi < ExprAi , ExprBi , typename PromotedType > Expr_Obj;
	return Expr1 < Expr_Obj, typename PromotedType , i > (Expr_Obj (ExprL, ExprR));
}

#undef ExprAi
#undef ExprBi
#undef PromotedType
#endif
