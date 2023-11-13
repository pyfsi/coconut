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
/* Contracts a Tensor2 by a Tensor 1, yielding a Tensor1.
 */

#ifndef Expr2_contract_Expr1_H
#define Expr2_contract_Expr1_H

//#define USE_ASSERT_Expr2
#ifdef USE_ASSERT_Expr2
#include <assert.h>
#endif


#define ExprAij   Expr2<A,T,i,j>
#define ExprAji   Expr2<A,T,j,i>
#define ExprBj    Expr1<B,U,j>
#define promotedType promote<T,U>::V

/*! \brief The operator that defines a contraction
*
*	This operator performs the contraction along i index and returns the value
*	A(i,j)*B(j) ==> C(i)
*/

template < class A, class B, class T, class U, char i , char j>
inline const Expr1 < const Aij_contracts_Bj < ExprAij , ExprBj , typename promotedType>,
				typename promotedType, i >
operator* (const ExprAij &a, const ExprBj &b)
{
	typedef const Aij_contracts_Bj < ExprAij , ExprBj , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i> (ExprObj (a,b));
}

/*! \brief The operator that defines a contraction
*
*	This operator performs the contraction along i index and returns the value
*	B(j)*A(i,j) ==> C(i)
*/


template < class A, class B, class T, class U, char i , char j>
inline const Expr1 < const Aij_contracts_Bj < ExprAij , ExprBj , typename promotedType>,
				typename promotedType, i >
operator* (const ExprBj &b , const ExprAij &a)
{
	typedef const Aij_contracts_Bj < ExprAij , ExprBj , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i> (ExprObj (a,b));
}

/*! \brief The operator that defines a contraction
*
*	This operator performs the contraction along j index and returns the value
*	A(j,i)*B(j)==> C(i)
*/


template < class A, class B, class T, class U, char i , char j>
inline const Expr1 < const Aji_contracts_Bj < ExprAji , ExprBj , typename promotedType>,
				typename promotedType, i >
operator* (const ExprAji &a, const ExprBj &b)
{
	typedef const Aji_contracts_Bj < ExprAji , ExprBj , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i> (ExprObj (a,b));
}
/*! \brief The operator that defines a contraction
*
*	This operator performs the contraction along j index and returns the value
*	B(j)*A(j,i)==> C(i)
*/

template < class A, class B, class T, class U, char i , char j>
inline const Expr1 < const Aji_contracts_Bj < ExprAji , ExprBj , typename promotedType>,
				typename promotedType, i >
operator* (const ExprBj &b, const ExprAji &a)
{
	typedef const Aji_contracts_Bj < ExprAji , ExprBj , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i> (ExprObj (a,b));
}


#undef ExprAij
#undef ExprAji
#undef ExprBj
#undef promotedType

#endif
