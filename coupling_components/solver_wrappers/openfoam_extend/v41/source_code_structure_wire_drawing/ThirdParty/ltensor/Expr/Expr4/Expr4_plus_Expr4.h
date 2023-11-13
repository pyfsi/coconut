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
#ifndef Expr4_plus_Expr4_H
#define Expr4_plus_Expr4_H


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
#define ExprBikjl   Expr4<B,U,i,k,j,l>
#define ExprBiklj   Expr4<B,U,i,k,l,j>
#define ExprBiljk   Expr4<B,U,i,l,j,k>
#define ExprBilkj   Expr4<B,U,i,l,k,j>
#define ExprBjkil   Expr4<B,U,j,k,i,l>
#define ExprBjkli   Expr4<B,U,j,k,l,i>
#define ExprBjlik   Expr4<B,U,j,l,i,k>
#define ExprBjlki   Expr4<B,U,j,l,k,i>
#define ExprBkijl   Expr4<B,U,k,i,j,l>
#define ExprBkilj   Expr4<B,U,k,i,l,j>
#define ExprBkjil   Expr4<B,U,k,j,i,l>
#define ExprBkjli   Expr4<B,U,k,j,l,i>
#define ExprBlijk   Expr4<B,U,l,i,j,k>
#define ExprBlikj   Expr4<B,U,l,i,k,j>
#define ExprBljik   Expr4<B,U,l,j,i,k>







#define promotedType promote<T,U>::V

///////////////////////////////////////////////
// Define SUM = A(i,j,k,l)+B(i,j,k,l)
///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k , char l>
inline const Expr4< const Aijkl_plus_Bijkl < ExprAijkl, ExprBijkl , typename promotedType>,
				typename promotedType, i, j, k, l >
operator+  (const ExprAijkl &ExprL, const ExprBijkl &ExprR)
{
	typedef const  Aijkl_plus_Bijkl <  ExprAijkl,  ExprBijkl , typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, ExprR));
}

//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(j,i,k,l)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Ajikl_to_Aijkl < ExprBjikl ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >
operator+  (const ExprAijkl &ExprL, const ExprBjikl &ExprR)
{
	typedef Ajikl_to_Aijkl < ExprBjikl ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(i,j,l,k)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Aijlk_to_Aijkl < ExprBijlk ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >
operator+  (const ExprAijkl &ExprL, const ExprBijlk &ExprR)
{
	typedef Aijlk_to_Aijkl < ExprBijlk ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(j,i,l,k)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Ajilk_to_Aijkl < ExprBjilk ,typename promotedType>, typename promotedType >,
				typename promotedType, i, j, k, l >
operator+  (const ExprAijkl &ExprL, const ExprBjilk &ExprR)
{
	typedef Ajilk_to_Aijkl < ExprBjilk ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(k,l,i,j)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Aklij_to_Aijkl < ExprBklij ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >
operator+  (const ExprAijkl &ExprL, const ExprBklij &ExprR)
{
	typedef Aklij_to_Aijkl < ExprBklij ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(l,k,i,j)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Alkij_to_Aijkl < ExprBlkij ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >
operator+  (const ExprAijkl &ExprL, const ExprBlkij &ExprR)
{
	typedef Alkij_to_Aijkl < ExprBlkij ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(l,k,j,i)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Alkji_to_Aijkl < ExprBlkji ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >
operator+  (const ExprAijkl &ExprL, const ExprBlkji &ExprR)
{
	typedef Alkji_to_Aijkl < ExprBlkji ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(k,l,j,i)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Aklji_to_Aijkl < ExprBklji ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >
operator+  (const ExprAijkl &ExprL, const ExprBklji &ExprR)
{
	typedef Aklji_to_Aijkl < ExprBklji ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(l,j,k,i)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Aljki_to_Aijkl < ExprBljki ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBljki &ExprR)
{
	typedef Aljki_to_Aijkl < ExprBljki ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}



//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(i,k,j,l)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Aikjl_to_Aijkl < ExprBikjl ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBikjl &ExprR)
{
	typedef Aikjl_to_Aijkl < ExprBikjl ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}


//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(i,k,l,j)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Aiklj_to_Aijkl < ExprBiklj ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBiklj &ExprR)
{
	typedef Aiklj_to_Aijkl < ExprBiklj ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}


//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(i,l,j,k)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Ailjk_to_Aijkl < ExprBiljk ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBiljk &ExprR)
{
	typedef Ailjk_to_Aijkl < ExprBiljk ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}



//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(i,l,k,j)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Ailkj_to_Aijkl < ExprBilkj ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBilkj &ExprR)
{
	typedef Ailkj_to_Aijkl < ExprBilkj ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}


//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(j,k,i,l)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Ajkil_to_Aijkl < ExprBjkil ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBjkil &ExprR)
{
	typedef Ajkil_to_Aijkl < ExprBjkil ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}



//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(j,k,l,i)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Ajkli_to_Aijkl < ExprBjkli ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBjkli &ExprR)
{
	typedef Ajkli_to_Aijkl < ExprBjkli ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}



//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(j,l,i,k)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Ajlik_to_Aijkl < ExprBjlik ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBjlik &ExprR)
{
	typedef Ajlik_to_Aijkl < ExprBjlik ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}


//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(j,l,k,i)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Ajlki_to_Aijkl < ExprBjlki ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBjlki &ExprR)
{
	typedef Ajlki_to_Aijkl < ExprBjlki ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}



//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(k,i,j,l)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Akijl_to_Aijkl < ExprBkijl ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBkijl &ExprR)
{
	typedef Akijl_to_Aijkl < ExprBkijl ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}



//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(k,i,l,j)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Akilj_to_Aijkl < ExprBkilj ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBkilj &ExprR)
{
	typedef Akilj_to_Aijkl < ExprBkilj ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}



//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(k,j,i,l)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Akjil_to_Aijkl < ExprBkjil ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBkjil &ExprR)
{
	typedef Akjil_to_Aijkl < ExprBkjil ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}


//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(k,j,l,i)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Akjli_to_Aijkl < ExprBkjli ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBkjli &ExprR)
{
	typedef Akjli_to_Aijkl < ExprBkjli ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}



//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(l,i,j,k)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Alijk_to_Aijkl < ExprBlijk ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBlijk &ExprR)
{
	typedef Alijk_to_Aijkl < ExprBlijk ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}



//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(l,i,k,j)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Alikj_to_Aijkl < ExprBlikj ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBlikj &ExprR)
{
	typedef Alikj_to_Aijkl < ExprBlikj ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
	return Expr4 < Expr_Obj, typename promotedType, i, j, k, l > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}



//~ ///////////////////////////////////////////////
//~ // Define SUM = A(i,j,k,l)+B(l,j,i,k)
//~ ///////////////////////////////////////////////

template < class A, class B, class T, class U, char i, char j, char k, char l >
inline const Expr4 < const Aijkl_plus_Bijkl < ExprAijkl,
				Aljik_to_Aijkl < ExprBljik ,U>, typename promotedType>,
				typename promotedType, i, j, k, l >

operator+  (const ExprAijkl &ExprL, const ExprBljik &ExprR)
{
	typedef Aljik_to_Aijkl < ExprBljik ,U> Permuted_Obj;
	typedef const  Aijkl_plus_Bijkl < ExprAijkl, Permuted_Obj, typename promotedType> Expr_Obj;
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

#undef ExprBikjl
#undef ExprBiklj
#undef ExprBiljk
#undef ExprBilkj
#undef ExprBjkil
#undef ExprBjkli
#undef ExprBjlik
#undef ExprBjlki
#undef ExprBkijl
#undef ExprBkilj
#undef ExprBkjil
#undef ExprBkjli
#undef ExprBlijk
#undef ExprBlikj
#undef ExprBljik


#undef promotedType

#endif
