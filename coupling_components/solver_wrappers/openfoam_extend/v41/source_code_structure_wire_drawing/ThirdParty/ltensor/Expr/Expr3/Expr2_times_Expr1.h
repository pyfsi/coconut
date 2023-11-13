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

#ifndef Expr2_times_Expr1_H
#define Expr2_times_Expr1_H

#define ExprAij	Expr2<A,T,i,j>
#define ExprBk	Expr1<B,U,k>
#define promotedType promote<T,U>::V


////////////////////////////////////////////
// Define Operations 
///////////////////////////////////////////


/*!\brief Product operation between an Expr of rank 2 and an Expr of rank 1
*	
*	Performs the * operation and returns an Expr of rank 3 containing the result
*	A(i,j) * B(k) -> C(i,j,k)
*/

template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr3 < const Aij_times_Bk < ExprAij, typename promotedType, ExprBk>, typename promotedType, i, j, k >
operator* (const ExprAij &a, const ExprBk &b)
{
	typedef const Aij_times_Bk< ExprAij, typename promotedType, ExprBk > ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, k> (ExprObj (a,b));
}


////////////////////////////////////////////
// Define Operations B(k)*A(i,j)  -> C(i,j,k) 
///////////////////////////////////////////


/*!\brief Product operation between an Expr of rank 2 and an Expr of rank 1
*	
*	Performs the * operation and returns an Expr of rank 3 containing the result
*	B(k)*A(i,j) -> C(i,j,k)
*/



template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr3 < const Aij_times_Bk < ExprAij, typename promotedType, ExprBk>, typename promotedType, i, j, k >
operator* (const ExprBk &b, const ExprAij &a)
{
	typedef const Aij_times_Bk< ExprAij, typename promotedType, ExprBk > ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, k> (ExprObj (a,b));
}


#undef ExprAij
#undef ExprBk
#undef promotedType

#endif
