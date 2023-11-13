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
#ifndef Expr4_contract_Expr1_H
#define Expr4_contract_Expr1_H


////////////////////////////////////////////
/*!brief A(i,j,k,l) * B(l) -> C(i,j,k) */
///////////////////////////////////////////

#define ExprAijkl	Expr4<A,T,i,j,k,l>
#define ExprAijlk	Expr4<A,T,i,j,l,k>
#define ExprAiljk	Expr4<A,T,i,l,j,k>
#define ExprAlijk	Expr4<A,T,l,i,j,k>

#define ExprBl	Expr1<B,U,l>
#define promotedType	promote<T,U>::V


/*!brief A(i,j,k,l)*B(l) ==> C(i,j,k) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijkl_contracts_Bl < ExprAijkl , ExprBl , typename promotedType>,
				typename promotedType, i , j, k >
operator* (const ExprAijkl &a, const ExprBl &b)
{
	typedef  const Aijkl_contracts_Bl < ExprAijkl , ExprBl , typename promotedType> ExprObj;
	return Expr3 < ExprObj,typename promotedType , i, j, k > (ExprObj (a,b));
}

/*!brief B(l)*A(i,j,k,l) ==> C(i,j,k) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijkl_contracts_Bl < ExprAijkl , ExprBl , typename promotedType>,
				typename promotedType, i , j, k >
operator* (const ExprBl &b, const ExprAijkl &a)
{
	typedef  const Aijkl_contracts_Bl < ExprAijkl , ExprBl , typename promotedType> ExprObj;
	return Expr3 < ExprObj,typename promotedType , i, j, k > (ExprObj (a,b));
}

/*!brief A(i,j,l,k)*B(l) ==> C(i,j,k) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijkl_contracts_Bl < Aijlk_to_Aijkl < ExprAijlk ,T > ,
			 ExprBl , typename promotedType>, typename promotedType, i, j, k >
operator* (const ExprAijlk &a, const ExprBl &b)
{
	typedef Aijlk_to_Aijkl < ExprAijlk ,T > Permuted_Obj;
	typedef const Aijkl_contracts_Bl < Permuted_Obj , ExprBl , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType , i, j, k > (ExprObj (Permuted_Obj(a),b));
}

/*!brief B(l)*A(i,j,l,k) ==> C(i,j,k) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijkl_contracts_Bl < Aijlk_to_Aijkl < ExprAijlk ,T > ,
			 ExprBl , typename promotedType>, typename promotedType, i, j, k >
operator* (const ExprBl &b, const ExprAijlk &a)
{
	typedef Aijlk_to_Aijkl < ExprAijlk ,T > Permuted_Obj;
	typedef const Aijkl_contracts_Bl < Permuted_Obj , ExprBl , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType , i, j, k > (ExprObj (Permuted_Obj(a),b));
}



/*!brief A(i,l,j,k)*B(l) ==> C(i,j,k) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijkl_contracts_Bl < Ailjk_to_Aijkl < ExprAiljk ,T > ,
			 ExprBl , typename promotedType>, typename promotedType, i, j, k >
operator* (const ExprAiljk &a, const ExprBl &b)
{
	typedef Ailjk_to_Aijkl < ExprAiljk ,T > Permuted_Obj;
	typedef const Aijkl_contracts_Bl < Permuted_Obj , ExprBl , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType , i, j, k > (ExprObj (Permuted_Obj(a),b));
}

/*!brief B(l)*A(i,l,j,k) ==> C(i,j,k) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijkl_contracts_Bl < Ailjk_to_Aijkl < ExprAiljk ,T > ,
			 ExprBl , typename promotedType>, typename promotedType, i, j, k >
operator* (const ExprBl &b,const ExprAiljk &a)
{
	typedef Ailjk_to_Aijkl < ExprAiljk ,T > Permuted_Obj;
	typedef const Aijkl_contracts_Bl < Permuted_Obj , ExprBl , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType , i, j, k > (ExprObj (Permuted_Obj(a),b));
}

/*!brief A(l,i,j,k)*B(l) ==> C(i,j,k) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijkl_contracts_Bl < Alijk_to_Aijkl < ExprAlijk ,T > ,
			 ExprBl , typename promotedType>, typename promotedType, i, j, k >
operator* (const ExprAlijk &a, const ExprBl &b)
{
	typedef Alijk_to_Aijkl < ExprAlijk ,T > Permuted_Obj;
	typedef const Aijkl_contracts_Bl < Permuted_Obj , ExprBl , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType , i, j, k > (ExprObj (Permuted_Obj(a),b));
}

/*!brief B(l)*A(l,i,j,k) ==> C(i,j,k) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijkl_contracts_Bl < Alijk_to_Aijkl < ExprAlijk ,T> ,
			 ExprBl , typename promotedType>, typename promotedType, i, j, k >
operator* (const ExprBl &b, const ExprAlijk &a)
{
	typedef Alijk_to_Aijkl < ExprAlijk ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Bl < Permuted_Obj , ExprBl , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType , i, j, k > (ExprObj (Permuted_Obj(a),b));
}

#undef ExprAijkl
#undef ExprAijlk
#undef ExprAiljk
#undef ExprAlijk
#undef ExprBl
#undef pType

#endif
