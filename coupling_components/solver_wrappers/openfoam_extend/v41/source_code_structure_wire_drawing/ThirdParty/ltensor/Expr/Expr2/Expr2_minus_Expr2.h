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
/* Adds two Tensor2's together, yielding a Tensor2. */

#ifndef Expr2_minus_Expr2_H
#define Expr2_minus_Expr2_H


#ifdef USE_ASSERT_Expr2
#include <assert.h>
#endif

#include "./Expr2.h"

#define ExprAij   Expr2<A,T,i,j>
#define ExprBij   Expr2<B,U,i,j>
#define ExprBji   Expr2<B,U,j,i>
#define promotedType promote<T,U>::V



/*!\brief Minus operation between two Expr Objects of rank 2
*	
*	Performs the - operation and returns another Expr of rank 2 containing the result
*	A(i,j)-B(i,j)
*/

template < class A, class B, class T,class U, char i, char j >
inline const Expr2 <  const Aij_minus_Bij < ExprAij , ExprBij , typename promotedType > , typename promotedType, i, j >
operator-  (const ExprAij &ExprL, const ExprBij &ExprR)
{
	typedef const Aij_minus_Bij < ExprAij , ExprBij , typename promotedType > Expr_Obj;
	return Expr2 < Expr_Obj, typename promotedType , i, j > (Expr_Obj (ExprL, ExprR));
}

///////////////////////////////////////////////
// Define DIFF = A(i,j)-B(j,i)
///////////////////////////////////////////////

/*!\brief Minus operation between two Expr Objects of rank 2
*	
*	Performs the - operation and returns another Expr of rank 2 containing the result
*	
*	\remarks NOTE the permutted indexes A(i,j)-B(j,i)
*/


template < class A, class B, class T,class U, char i, char j >
inline const Expr2 <  const Aij_minus_Bji < ExprAij , ExprBji, typename promotedType >, typename promotedType, i, j >
operator-  (const ExprAij &ExprL, const ExprBji &ExprR)
{
	typedef const Aij_minus_Bji < ExprAij , ExprBji , T > Expr_Obj;
	return Expr2 < Expr_Obj, typename promotedType , i, j > (Expr_Obj (ExprL, ExprR));
}

#undef ExprAij
#undef ExprBij
#undef ExprBji
#undef promotedType
#endif
