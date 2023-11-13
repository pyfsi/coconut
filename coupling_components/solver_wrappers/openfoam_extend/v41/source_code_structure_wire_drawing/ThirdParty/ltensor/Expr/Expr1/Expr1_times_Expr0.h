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

/* Multipliess a Tensor1 by a generic (or vice versa), yielding a Tensor1.
   Usually used for doubles, but could be used for complex, etc.  All
   that it requires is that you can add an element of the Tensor1 to
   it.  */

#ifndef Expr1_times_Expr0_H
#define Expr1_times_Expr0_H

#define ExprAi  Expr1<A,T,i>
#define promotedType promote<T,U>::V



/*!\brief Product operation between an Expr of rank 1 and a scalar 
*	
*	Performs the * operation and returns an Expr of rank 1 containing the result
*	A(i) * u0 -> Tensor1
*/

template < class A, class T, class U, char i >
inline const Expr1 < const Ai_times_u < ExprAi, typename promotedType, U>, typename promotedType, i >
operator* (const ExprAi &a, const U & u0)
{
	typedef const Ai_times_u < ExprAi, typename promotedType, U> ExprObj;
	return Expr1 < ExprObj,typename promotedType,i > (ExprObj (a, u0));
}

/* u0 * A(i) -> Tensor1 */


/*!\brief Product operation between an Expr of rank 1 and a scalar 
*	
*	Performs the * operation and returns an Expr of rank 1 containing the result
*	A(i) * u0 -> Tensor1
*/
template < class A, class T, class U, char i >
inline const Expr1 < const Ai_times_u < ExprAi, typename promotedType, U>, typename promotedType, i >
operator* ( const U & u0, const ExprAi &a)
{
	typedef const Ai_times_u < ExprAi,typename promotedType, U> ExprObj;
	return Expr1 < ExprObj, typename promotedType,i > (ExprObj (a, u0));
}

#undef ExprAi
#undef promotedType
#endif
