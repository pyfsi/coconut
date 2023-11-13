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

#ifndef Expr4_divided_Expr0_H
#define Expr4_divided_Expr0_H

#define ExprAijkl   Expr4<A,T,i,j,k,l>
#define promotedType promote<T,U>::V

// Define:
// A(i,j,k,l) / u0 -> Tensor4


template < class A, class T, class U, char i , char j, char k, char l>
inline const Expr4 < const Aijkl_divided_u < ExprAijkl, typename promotedType, U>,
				typename promotedType, i, j, k, l>
operator/ (const ExprAijkl &a, const U & u0)
{
	typedef const Aijkl_divided_u < ExprAijkl, typename promotedType, U> ExprObj;
	return Expr4 < ExprObj, typename promotedType, i, j, k, l> (ExprObj (a, u0));
}


#undef ExprAijkl
#undef promotedType
#endif
