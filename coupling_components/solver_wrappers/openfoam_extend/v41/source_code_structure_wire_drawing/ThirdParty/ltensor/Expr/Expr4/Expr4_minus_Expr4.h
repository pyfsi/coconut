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
#ifndef Expr4_minus_Expr4_H
#define Expr4_minus_Expr4_H


#ifdef USE_ASSERT_Expr4
#include <assert.h>
#endif


#define ExprAijkl   Expr4<A,T,i,j,k,l>

#define ExprBijkl   Expr4<B,U,i,j,k,l>
#define ExprBjikl   Expr4<B,U,j,i,k,l>
#define ExprBijlk   Expr4<B,U,i,j,l,k>
#define ExprBjilk   Expr4<B,U,j,i,l,k>
#define ExprBklij   Expr4<B,U,k,l,i,j>
#define ExprBlkij   Expr4<B,U,l,k,i,j>
#define ExprBlkji   Expr4<B,U,l,k,j,i>
#define ExprBklji   Expr4<B,U,k,l,j,i>
#define ExprBljki   Expr4<B,U,l,j,k,i>

#define promotedType promote<T,U>::V

///////////////////////////////////////////////
// Define MINUS = A(i,j,k,l)+B(i,j,k,l)
///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k , char l>
inline const Expr4< const Aijkl_minus_Bijkl < ExprAijkl, ExprBijkl , typename promotedType>,
				typename promotedType, i, j, k, l >

operator-  (const ExprAijkl &ExprL, const ExprBijkl &ExprR)
{
	typedef const  Aijkl_minus_Bijkl <  ExprAijkl,  ExprBijkl , typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, ExprR));
}

//~ ///////////////////////////////////////////////
//~ // Define MINUS = A(i,j,k,l)+B(j,i,k,l)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_minus_Bijkl < ExprAijkl,
				Ajikl_to_Aijkl < ExprBjikl ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >
operator-  (const ExprAijkl &ExprL, const ExprBjikl &ExprR)
{
	typedef Ajikl_to_Aijkl < ExprBjikl ,U> Permuted_Obj;
	typedef const  Aijkl_minus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

//~ ///////////////////////////////////////////////
//~ // Define MINUS = A(i,j,k,l)+B(i,j,l,k)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_minus_Bijkl < ExprAijkl,
				Aijlk_to_Aijkl < ExprBijlk ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >
operator-  (const ExprAijkl &ExprL, const ExprBijlk &ExprR)
{
	typedef Aijlk_to_Aijkl < ExprBijlk ,U> Permuted_Obj;
	typedef const  Aijkl_minus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

//~ ///////////////////////////////////////////////
//~ // Define MINUS = A(i,j,k,l)+B(j,i,l,k)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_minus_Bijkl < ExprAijkl,
				Ajilk_to_Aijkl < ExprBjilk ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >
operator-  (const ExprAijkl &ExprL, const ExprBjilk &ExprR)
{
	typedef Ajilk_to_Aijkl < ExprBjilk ,U> Permuted_Obj;
	typedef const  Aijkl_minus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

//~ ///////////////////////////////////////////////
//~ // Define MINUS = A(i,j,k,l)+B(k,l,i,j)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_minus_Bijkl < ExprAijkl,
				Aklij_to_Aijkl < ExprBklij ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >
operator-  (const ExprAijkl &ExprL, const ExprBklij &ExprR)
{
	typedef Aklij_to_Aijkl < ExprBklij ,U> Permuted_Obj;
	typedef const  Aijkl_minus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

//~ ///////////////////////////////////////////////
//~ // Define MINUS = A(i,j,k,l)+B(l,k,i,j)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_minus_Bijkl < ExprAijkl,
				Alkij_to_Aijkl < ExprBlkij ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >
operator-  (const ExprAijkl &ExprL, const ExprBlkij &ExprR)
{
	typedef Alkij_to_Aijkl < ExprBlkij ,U> Permuted_Obj;
	typedef const  Aijkl_minus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

//~ ///////////////////////////////////////////////
//~ // Define MINUS = A(i,j,k,l)+B(l,k,j,i)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_minus_Bijkl < ExprAijkl,
				Alkji_to_Aijkl < ExprBlkji ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >
operator-  (const ExprAijkl &ExprL, const ExprBlkji &ExprR)
{
	typedef Alkji_to_Aijkl < ExprBlkji ,U> Permuted_Obj;
	typedef const  Aijkl_minus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

//~ ///////////////////////////////////////////////
//~ // Define MINUS = A(i,j,k,l)+B(k,l,j,i)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_minus_Bijkl < ExprAijkl,
				Aklji_to_Aijkl < ExprBklji ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >
operator-  (const ExprAijkl &ExprL, const ExprBklji &ExprR)
{
	typedef Aklji_to_Aijkl < ExprBklji ,U> Permuted_Obj;
	typedef const  Aijkl_minus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

//~ ///////////////////////////////////////////////
//~ // Define MINUS = A(i,j,k,l)+B(l,j,k,i)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_minus_Bijkl < ExprAijkl,
				Aljki_to_Aijkl < ExprBljki ,typename promotedType>, typename promotedType >,
				typename promotedType, i, j, k, l >
operator-  (const ExprAijkl &ExprL, const ExprBljki &ExprR)
{
	typedef Aljki_to_Aijkl < ExprBljki ,U> Permuted_Obj;
	typedef const  Aijkl_minus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

#undef ExprAijkl

#undef ExprBijkl
#undef ExprBjikl
#undef ExprBijlk
#undef ExprBjilk
#undef ExprBklij
#undef ExprBlkij
#undef ExprBlkji
#undef ExprBklji
#undef ExprBljki

#undef pType

#endif
