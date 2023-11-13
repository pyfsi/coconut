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

#ifndef Expr1_times_Expr1_H
#define Expr1_times_Expr1_H

#define ExprAi	Expr1<A,T,i>
#define ExprBj	Expr1<B,U,j>
#define promotedType promote<T,U>::V


////////////////////////////////////////////
/*!brief A(i) * B(j) -> C(i,j) */ 
///////////////////////////////////////////

/*\brief Product operation between two Expr of rank 1 with different indexes
*	
*	Performs the * operation and returns an Expr of Rank 2 with the result
*	A(i) * B(j) -> C(i,j) 
*/

template < class A, class B, class T, class U, char i , char j>
inline const Expr2 < const Ai_times_Bj < ExprAi, typename promotedType, ExprBj>, typename promotedType, i, j >
operator* (const ExprAi &a, const ExprBj &b)
{
	typedef const Ai_times_Bj< ExprAi, typename promotedType, ExprBj > ExprObj;
	return Expr2 < ExprObj, typename promotedType, i, j> (ExprObj (a,b));
}


#undef ExprAi
#undef ExprBj
#undef pType

#endif
