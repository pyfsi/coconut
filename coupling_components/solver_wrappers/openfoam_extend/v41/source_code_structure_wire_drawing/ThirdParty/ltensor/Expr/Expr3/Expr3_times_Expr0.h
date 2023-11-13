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

#ifndef Expr3_times_Expr0_H
#define Expr3_times_Expr0_H

#define ExprAijk   Expr3<A,T,i,j,k>
#define promotedType promote<T,U>::V

// Define:
// A(i,j,k) * u0 -> Tensor3

/*!\brief Product operation between an Expr of rank 3 and a scalar 
*	
*	Performs the * operation and returns an Expr of rank 3 containing the result
*	A(i,j,k) * u0 -> Tensor2
*/

template < class A, class T, class U, char i , char j, char k>
inline const Expr3 < const Aijk_times_u < ExprAijk, typename promotedType, U>,
				typename promotedType, i, j, k >
operator* (const ExprAijk &a, const U & u0)
{
	typedef const Aijk_times_u < ExprAijk, typename promotedType, U> ExprObj;
	return Expr3 < ExprObj, typename promotedType , i, j, k > (ExprObj (a, u0));
}

// Define
// u0 * A(i,j) -> Tensor3


/*!\brief Product operation between an Expr of rank 3 and a scalar 
*	
*	Performs the * operation and returns an Expr of rank 3 containing the result
*	u0 * A(i,j,k)  -> Tensor2
*/

template < class A, class T, class U, char i , char j, char k>
inline const Expr3 < const Aijk_times_u < ExprAijk, typename promotedType, U>,
				typename promotedType, i, j, k >
operator* (const U & u0, const ExprAijk &a)
{
	typedef const Aijk_times_u < ExprAijk, typename promotedType, U> ExprObj;
	return Expr3 < ExprObj, typename promotedType , i, j, k > (ExprObj (a, u0));
}

#undef ExprAijk
#undef pType
#endif
