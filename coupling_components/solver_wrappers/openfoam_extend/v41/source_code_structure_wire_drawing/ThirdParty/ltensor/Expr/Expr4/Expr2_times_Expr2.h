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

#ifndef Expr2_times_Expr2_H
#define Expr2_times_Expr2_H

#define ExprAij	Expr2<A,T,i,j>
#define ExprBkl	Expr2<B,U,k,l>
#define promotedType promote<T,U>::V




/*!\brief Product operation between two Expr of rank 2
*	
*	Performs the * operation and returns an Expr of rank 4 containing the result
*	A(i,j) * B(k,l) -> C(i,j,k.l)
*/


template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr4 < const Aij_times_Bkl < ExprAij, typename promotedType, ExprBkl>, typename promotedType, i, j, k, l >
operator* (const ExprAij &a, const ExprBkl &b)
{
	typedef const Aij_times_Bkl< ExprAij, typename promotedType, ExprBkl > ExprObj;
	return Expr4 < ExprObj, typename promotedType, i, j, k,l> (ExprObj (a,b));
}



#undef ExprAij
#undef ExprBkl
#undef pType

#endif
