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
//TODO check if complete regarding permutations


#ifndef Expr2_plus_Expr2_H
#define Expr2_plus_Expr2_H


//#define CHECK_Expr2_plus_Expr2

#ifdef USE_ASSERT_Expr2
#include <assert.h>
#endif

#include "./Expr2.h"

#define ExprAij   Expr2<A,T,i,j>
#define ExprBij   Expr2<B,U,i,j>
#define ExprBji   Expr2<B,U,j,i>
#define pType promote<T,U>::V


/*!\brief Plus operation between two Expr Objects of rank 2
*	
*	Performs the + operation and returns another Expr of rank 2 containing the result
*	A(i,j)+B(i,j)
*/


template < class A, class B, class T, class U, char i, char j >
inline const Expr2 <  const Aij_plus_Bij < ExprAij , ExprBij , typename pType > ,
				typename pType, i, j >
operator+  (const ExprAij &ExprL, const ExprBij &ExprR)
{
	typedef const Aij_plus_Bij <  ExprAij , ExprBij , typename pType > Expr_Obj;
	return Expr2 < Expr_Obj,typename pType , i, j > (Expr_Obj (ExprL, ExprR));
}

///////////////////////////////////////////////
// Define SUM = A(i,j)+B(j,i)
///////////////////////////////////////////////


/*!\brief Plus operation between two Expr Objects of rank 2
*	
*	Performs the + operation and returns another Expr of rank 2 containing the result
*	\remarks
*	NOTE the permutted indexes A(i,j)+B(j,i)
*/

template < class A, class B, class T, class U, char i, char j >
inline const Expr2 < const Aij_plus_Bji < ExprAij , ExprBji , typename pType > ,
				typename pType, i, j >
operator+  (const ExprAij &ExprL, const ExprBji &ExprR)
{
	typedef const Aij_plus_Bji <  ExprAij , ExprBji , typename pType > Expr_Obj;
	return Expr2 < Expr_Obj,typename pType , i, j > (Expr_Obj (ExprL, ExprR));
}

#undef ExprAij
#undef ExprBij
#undef ExprBji
#undef pType
#endif
