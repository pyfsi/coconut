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
#ifndef Expr3_plus_Expr3_H
#define Expr3_plus_Expr3_H


#ifdef USE_ASSERT_Expr3
#include <assert.h>
#endif


#define ExprAijk   Expr3<A,T,i,j,k>

#define ExprBijk   Expr3<B,U,i,j,k>
#define ExprBikj   Expr3<B,U,i,k,j>
#define ExprBjik   Expr3<B,U,j,i,k>
#define ExprBkij   Expr3<B,U,k,i,j>
#define ExprBjki   Expr3<B,U,j,k,i>
#define ExprBkji   Expr3<B,U,k,j,i>

#define promotedType promote<T,U>::V



/*!\brief Plus operation between two Expr Objects of rank 3
*	
*	Performs the + operation and returns another Expr of rank 3 containing the result
*	 A(i,j,k)+B(i,j,k)
*/

template < class A, class B, class T, class U, char i, char j, char k >
inline const Expr3< const Aijk_plus_Bijk < ExprAijk, ExprBijk , typename promotedType >,
				typename promotedType, i, j, k >
operator+  (const ExprAijk &ExprL, const ExprBijk &ExprR)
{
	typedef const  Aijk_plus_Bijk <  ExprAijk,  ExprBijk , typename promotedType> Expr_Obj;
	return Expr3 < Expr_Obj, typename promotedType, i, j, k > (Expr_Obj (ExprL, ExprR));
}


/*!\brief Plus operation between two Expr Objects of rank 3
*	
*	Performs the + operation and returns another Expr of rank 3 containing the result
*	 A(i,j,k)+B(i,k,j)
*/

template < class A, class B, class T, class U, char i, char j, char k >
inline const Expr3 < const Aijk_plus_Bijk < ExprAijk,
				Aikj_to_Aijk < ExprBikj ,U>, typename promotedType >,
				typename promotedType, i, j, k >
operator+  (const ExprAijk &ExprL, const ExprBikj &ExprR)
{
	typedef Aikj_to_Aijk < ExprBikj ,U > Permuted_Obj;
	typedef const   Aijk_plus_Bijk < ExprAijk, Permuted_Obj, typename promotedType > Expr_Obj;
	return Expr3 < Expr_Obj, typename promotedType, i, j, k > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}


/*!\brief Plus operation between two Expr Objects of rank 3
*	
*	Performs the + operation and returns another Expr of rank 3 containing the result
*	 A(i,j,k)+B(j,i,k)
*/
template < class A, class B, class T, class U, char i, char j, char k >
inline const Expr3 < const Aijk_plus_Bijk < ExprAijk,
				Ajik_to_Aijk < ExprBjik ,U>, typename promotedType >,
				typename promotedType, i, j, k >
operator+  (const ExprAijk &ExprL, const ExprBjik &ExprR)
{
	typedef Ajik_to_Aijk < ExprBjik ,U> Permuted_Obj;
	typedef const   Aijk_plus_Bijk < ExprAijk, Permuted_Obj, typename promotedType > Expr_Obj;
	return Expr3 < Expr_Obj, typename promotedType, i, j, k > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}

/*!\brief Plus operation between two Expr Objects of rank 3
*	
*	Performs the + operation and returns another Expr of rank 3 containing the result
*	A(i,j,k)+B(k,i,j)
*/
template < class A, class B, class T, class U, char i, char j, char k >
inline const Expr3 < const Aijk_plus_Bijk < ExprAijk,
				Akij_to_Aijk < ExprBkij ,U>, typename promotedType >,
				typename promotedType, i, j, k >
operator+  (const ExprAijk &ExprL, const ExprBkij &ExprR)
{
	typedef Akij_to_Aijk < ExprBkij ,U> Permuted_Obj;
	typedef const   Aijk_plus_Bijk < ExprAijk, Permuted_Obj, typename promotedType > Expr_Obj;
	return Expr3 < Expr_Obj, typename promotedType, i, j, k > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}


/*!\brief Plus operation between two Expr Objects of rank 3
*	
*	Performs the + operation and returns another Expr of rank 3 containing the result
*	A(i,j,k)+B(j,k,i)
*/
template < class A, class B, class T, class U, char i, char j, char k >
inline const Expr3 < const Aijk_plus_Bijk < ExprAijk,
				Ajki_to_Aijk < ExprBjki ,U>, typename promotedType >,
				typename promotedType, i, j, k >
operator+  (const ExprAijk &ExprL, const ExprBjki &ExprR)
{
	typedef Ajki_to_Aijk < ExprBjki ,U> Permuted_Obj;
	typedef const Aijk_plus_Bijk < ExprAijk, Permuted_Obj, typename promotedType > Expr_Obj;
	return Expr3 < Expr_Obj, typename promotedType, i, j, k > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}


/*!\brief Plus operation between two Expr Objects of rank 3
*	
*	Performs the + operation and returns another Expr of rank 3 containing the result
*	A(i,j,k)+B(k,j,i)
*/

template < class A, class B, class T, class U, char i, char j, char k >
inline const Expr3 < const Aijk_plus_Bijk < ExprAijk,
				Akji_to_Aijk < ExprBkji ,U >, typename promotedType >,
				typename promotedType, i, j, k >
operator+  (const ExprAijk &ExprL, const ExprBkji &ExprR)
{
	typedef Akji_to_Aijk < ExprBkji ,U> Permuted_Obj;
	typedef const Aijk_plus_Bijk < ExprAijk, Permuted_Obj, typename promotedType > Expr_Obj;
	return Expr3 < Expr_Obj, typename promotedType, i, j, k > (Expr_Obj (ExprL, Permuted_Obj(ExprR)));
}



#undef ExprAijk
#undef ExprBijk
#undef ExprBikj
#undef ExprBjik
#undef ExprBkij
#undef ExprBjki
#undef ExprBkji
#undef pType

#endif
