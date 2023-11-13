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
#ifndef minus_Expr4_H
#define minus_Expr4_H

#define ExprAijkl   Expr4<A,T,i,j,k,l>

template < class A, class T, char i , char j, char k, char l>
inline const Expr4 < const minus_Aijkl < ExprAijkl, T>, T, i, j, k, l >
operator- (const ExprAijkl &a)
{
	typedef const minus_Aijkl < ExprAijkl, T>  ExprObj;
	return Expr4 < ExprObj, T, i, j, k, l > ( ExprObj(a) );
}

#undef ExprAijkl

#endif
