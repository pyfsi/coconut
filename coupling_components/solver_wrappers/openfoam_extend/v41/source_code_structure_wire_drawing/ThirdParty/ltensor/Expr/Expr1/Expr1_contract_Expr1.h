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
/* Adds two Tensor1's together, yielding a Tensor1. */
#ifndef Expr1_contract_Expr1_H
#define Expr1_contract_Expr1_H


#define ExprAi   Expr1<A,T,i>
#define ExprBi   Expr1<B,U,i>
#define promotedType promote<T,U>::V
/*! \brief The operator that defines a contraction
*
*	This operator performs the contraction along i index and returns the value
*	A(i)*B(i)
*/
template < class A, class B, class T, class U, char i >
inline const typename promotedType
operator*  (const ExprAi &ExprL, const ExprBi &ExprR)
{     
	return Ai_contracts_Bi<ExprAi,ExprBi,typename promotedType >(ExprL,ExprR);
}



#undef ExprAi
#undef ExprBi
#undef promotedType


#endif
