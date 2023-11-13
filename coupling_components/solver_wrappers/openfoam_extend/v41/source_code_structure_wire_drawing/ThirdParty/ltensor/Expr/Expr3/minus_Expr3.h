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
#ifndef minus_Expr3_H
#define minus_Expr3_H

#define ExprAijk   Expr3<A,T,i,j,k>

template < class A, class T, char i , char j, char k>
inline const Expr3 < const minus_Aijk < ExprAijk, T>, T, i, j, k >
operator- (const ExprAijk &a)
{
	typedef const minus_Aijk < ExprAijk, T>  ExprObj;
	return Expr3 < ExprObj, T, i, j,k > ( ExprObj(a) );
}

#undef ExprAijk

#endif
