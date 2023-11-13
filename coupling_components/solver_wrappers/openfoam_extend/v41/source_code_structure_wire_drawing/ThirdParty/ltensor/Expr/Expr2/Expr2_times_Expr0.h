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
/* Multipliess a Tensor2 by a generic (or vice versa), yielding a Tensor2.
   Usually used for doubles, but could be used for complex, etc.  All
   that it requires is that you can add an element of the Tensor1 to
   it.  */

#ifndef Expr2_times_Expr0_H
#define Expr2_times_Expr0_H

//#define CHECK_Expr2_times_Expr0

/*!\brief Product operation between an Expr of rank 2 and a scalar 
*	
*	Performs the * operation and returns an Expr of rank 2 containing the result
*	A(i,j) * u0 -> Tensor2
*/

#define ExprAij  Expr2<A,T,i,j>
#define promotedType promote<T,U>::V

template < class A, class T, class U, char i, char j >
inline const Expr2 < const Aij_times_u < ExprAij, typename promotedType , U>,typename promotedType, i, j >
operator* (const ExprAij &a, const U & u0)
{
	typedef const Aij_times_u < ExprAij, typename promotedType, U> ExprObj;
	return Expr2 < ExprObj, typename promotedType, i, j > (ExprObj (a, u0));
}




/*!\brief Product operation between an Expr of rank 2 and a scalar 
*	
*	Performs the * operation and returns an Expr of rank 2 containing the result
*	u0 * A(i,j) -> Tensor2
*/

template < class A, class T, class U, char i, char j >
inline const Expr2 < const Aij_times_u < ExprAij, typename promotedType , U>, typename promotedType , i, j >
operator* (const U & u0, const ExprAij &a)
{
	typedef const Aij_times_u < ExprAij,typename promotedType, U> ExprObj;
	return Expr2 < ExprObj,typename promotedType,i,j > (ExprObj (a, u0));
}


#undef ExprAij
#undef promotedType

#endif
