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
#ifndef Expr3_contract_Expr1_H
#define Expr3_contract_Expr1_H


////////////////////////////////////////////
/*!\ brief A(i,j,k) * B(k) -> C(i,j) */
///////////////////////////////////////////

#define ExprAijk   Expr3<A,T,i,j,k>
#define ExprAikj   Expr3<A,T,i,k,j>
#define ExprAkij   Expr3<A,T,k,i,j>
#define ExprBk     Expr1<B,U,k>
#define promotedType promote<T,U>::V


template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr2 < const Aijk_contracts_Bk < ExprAijk , ExprBk , typename promotedType>,
				typename promotedType, i , j>
operator* (const ExprAijk &a, const ExprBk &b)
{
	typedef  const Aijk_contracts_Bk < ExprAijk , ExprBk , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , i, j> (ExprObj (a,b));
}

/*!\ brief B(k)*A(i,j,k) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr2 < const Aijk_contracts_Bk < ExprAijk , ExprBk , typename promotedType>,
				typename promotedType, i , j>
operator* (const ExprBk &b, const ExprAijk &a)
{
	typedef  const Aijk_contracts_Bk < ExprAijk , ExprBk , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , i, j> (ExprObj (a,b));
}

/*!\ brief A(i,k,j)*B(k) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr2 < const Aijk_contracts_Bk < Aikj_to_Aijk < ExprAikj ,T > ,
			 ExprBk , typename promotedType>,
				typename promotedType, i , j>
operator* (const ExprAikj &a, const ExprBk &b)
{
	typedef Aikj_to_Aijk < ExprAikj ,typename promotedType > Permuted_Obj;
	typedef const Aijk_contracts_Bk < Permuted_Obj , ExprBk , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , i, j> (ExprObj (Permuted_Obj(a),b));
}

/*!\ brief B(k)*A(i,k,j) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr2 < const Aijk_contracts_Bk < Aikj_to_Aijk < ExprAikj ,T > ,
			 ExprBk , typename promotedType>,
				typename promotedType, i , j>
operator* (const ExprBk &b, const ExprAikj &a)
{
	typedef Aikj_to_Aijk < ExprAikj ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bk < Permuted_Obj , ExprBk , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , i, j> (ExprObj (Permuted_Obj(a),b));
}

/*!\ brief A(k,i,j)*B(k) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr2 < const Aijk_contracts_Bk < Akij_to_Aijk < ExprAkij ,T > ,
			 ExprBk , typename promotedType>,
				typename promotedType, i , j>
operator* (const ExprAkij &a, const ExprBk &b)
{
	typedef Akij_to_Aijk < ExprAkij ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bk < Permuted_Obj , ExprBk , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , i, j> (ExprObj (Permuted_Obj(a),b));
}

/*!\ brief B(k)*A(k,i,j) ==> C(i,j) */
//

template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr2 < const Aijk_contracts_Bk < Akij_to_Aijk < ExprAkij ,T > ,
			 ExprBk , typename promotedType>,
				typename promotedType, i , j>
operator* (const ExprBk &b, const ExprAkij &a )
{
	typedef Akij_to_Aijk < ExprAkij ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bk < Permuted_Obj , ExprBk , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , i, j> (ExprObj (Permuted_Obj(a),b));
}

#undef ExprAijk
#undef ExprAikj
#undef ExprAkij
#undef ExprBk
#undef promotedType

#endif
