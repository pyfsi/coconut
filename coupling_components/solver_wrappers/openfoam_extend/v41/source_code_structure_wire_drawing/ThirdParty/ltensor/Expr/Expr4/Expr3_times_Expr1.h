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

#ifndef Expr3_times_Expr1_H
#define Expr3_times_Expr1_H

#define ExprAijk Expr3<A,T,i,j,k>
#define ExprBl Expr1<B,U,l>
#define promotedType promote<T,U>::V


////////////////////////////////////////////
/*!brief A(i,j,k) * B(l) -> C(i,j,k,l) */
///////////////////////////////////////////

/*!\brief Product operation between an Expr of rank 3 and an Expr of rank 1
*	
*	Performs the * operation and returns an Expr of rank 4 containing the result
*	A(i,j,k) * B(l) -> C(i,j,k,l)
*/


template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr4 < const Aijk_times_Bl < ExprAijk, typename promotedType, ExprBl>, typename promotedType, i, j, k, l >
operator* (const ExprAijk &a, const ExprBl &b)
{
	typedef const Aijk_times_Bl< ExprAijk, typename promotedType, ExprBl > ExprObj;
	return Expr4 < ExprObj, typename promotedType, i, j, k,l> (ExprObj (a,b));
}


////////////////////////////////////////////
// Define Operations B(l)*A(i,j,k)  -> C(i,j,k,l)
///////////////////////////////////////////

/*!\brief Product operation between an Expr of rank 3 and an Expr of rank 1
*	
*	Performs the * operation and returns an Expr of rank 4 containing the result
*	B(l) * A(i,j,k)  -> C(i,j,k,l)
*/

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr4 < const Aijk_times_Bl < ExprAijk, typename promotedType, ExprBl>, typename promotedType, i, j, k,l >
operator* (const ExprBl &b, const ExprAijk &a)
{
	typedef const Aijk_times_Bl< ExprAijk, typename promotedType, ExprBl > ExprObj;
	return Expr4 < ExprObj, typename promotedType, i, j, k,l> (ExprObj (a,b));
}


#undef ExprAijk
#undef ExprBl
#undef pType

#endif
