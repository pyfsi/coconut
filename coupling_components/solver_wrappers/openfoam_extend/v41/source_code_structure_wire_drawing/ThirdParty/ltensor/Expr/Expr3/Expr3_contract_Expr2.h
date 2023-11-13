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
#ifndef Expr3_contract_Expr2_H
#define Expr3_contract_Expr2_H


////////////////////////////////////////////
/*!\ brief A(i,j,k) * B(j,k) -> C(i) */
///////////////////////////////////////////

#define ExprAijk   Expr3<A,T,i,j,k>
#define ExprAikj   Expr3<A,T,i,k,j>
#define ExprAjik   Expr3<A,T,j,i,k>
#define ExprAkij   Expr3<A,T,k,i,j>
#define ExprAjki   Expr3<A,T,j,k,i>
#define ExprAkji   Expr3<A,T,k,j,i>

#define ExprBjk     Expr2<B,U,j,k>
#define ExprBkl     Expr2<B,U,k,l>
#define ExprBlk     Expr2<B,U,l,k>
#define promotedType promote<T,U>::V


//////////////////////////////////////////////////////////////////////////
// Contraction of two indices
//////////////////////////////////////////////////////////////////////////

/*!\ brief A(i,j,k)*B(j,k) ==> C(i) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr1 < const Aijk_contracts_Bjk < ExprAijk , ExprBjk , typename promotedType>,
				typename promotedType, i >
operator* (const ExprAijk &a, const ExprBjk &b)
{
	typedef  const Aijk_contracts_Bjk < ExprAijk , ExprBjk , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i> (ExprObj (a,b));
}

/*!\ brief B(j,k)*A(i,j,k) ==> C(i) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr1 < const Aijk_contracts_Bjk < ExprAijk , ExprBjk , typename promotedType>,
				typename promotedType, i >
operator* (const ExprBjk &b, const ExprAijk &a)
{
	typedef  const Aijk_contracts_Bjk < ExprAijk , ExprBjk , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i> (ExprObj (a,b));
}

/*!\ brief A(i,k,j)*B(j,k) ==> C(i) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr1 < const Aijk_contracts_Bjk < Aikj_to_Aijk < ExprAikj ,T > ,
			 ExprBjk , typename promotedType>, typename promotedType, i>
operator* (const ExprAikj &a, const ExprBjk &b)
{
	typedef Aikj_to_Aijk < ExprAikj ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bjk < Permuted_Obj , ExprBjk , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i > (ExprObj (Permuted_Obj(a),b));
}

/*!\ brief B(j,k)*A(i,k,j) ==> C(i) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr1 < const Aijk_contracts_Bjk < Aikj_to_Aijk < ExprAikj ,T > ,
			 ExprBjk , typename promotedType>, typename promotedType, i>
operator* (const ExprBjk &b, const ExprAikj &a)
{
	typedef Aikj_to_Aijk < ExprAikj ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bjk < Permuted_Obj , ExprBjk , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i > (ExprObj (Permuted_Obj(a),b));
}

/*!\ brief A(j,i,k)*B(j,k) ==> C(i) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr1 < const Aijk_contracts_Bjk < Ajik_to_Aijk < ExprAjik ,T > ,
			 ExprBjk , typename promotedType>, typename promotedType, i>
operator* (const ExprAjik &a, const ExprBjk &b)
{
	typedef Ajik_to_Aijk < ExprAjik ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bjk < Permuted_Obj , ExprBjk , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i > (ExprObj (Permuted_Obj(a),b));
}

/*!\ brief B(j,k)*A(j,i,k) ==> C(i) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr1 < const Aijk_contracts_Bjk < Ajik_to_Aijk < ExprAjik ,T > ,
			 ExprBjk , typename promotedType>, typename promotedType, i>
operator* (const ExprBjk &b, const ExprAjik &a)
{
	typedef Ajik_to_Aijk < ExprAjik ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bjk < Permuted_Obj , ExprBjk , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i > (ExprObj (Permuted_Obj(a),b));
}

/*!\ brief A(k,i,j)*B(j,k) ==> C(i) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr1 < const Aijk_contracts_Bjk < Akij_to_Aijk < ExprAkij ,T > ,
			 ExprBjk , typename promotedType>, typename promotedType, i>
operator* (const ExprAkij &a, const ExprBjk &b)
{
	typedef Akij_to_Aijk < ExprAkij ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bjk < Permuted_Obj , ExprBjk , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i > (ExprObj (Permuted_Obj(a),b));
}

/*!\ brief B(j,k)*A(k,i,j) ==> C(i) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr1 < const Aijk_contracts_Bjk < Akij_to_Aijk < ExprAkij ,T > ,
			 ExprBjk , typename promotedType>, typename promotedType, i>
operator* (const ExprBjk &b , const ExprAkij &a)
{
	typedef Akij_to_Aijk < ExprAkij ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bjk < Permuted_Obj , ExprBjk , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i > (ExprObj (Permuted_Obj(a),b));
}

/*!\ brief A(j,k,i)*B(j,k) ==> C(i) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr1 < const Aijk_contracts_Bjk < Ajki_to_Aijk < ExprAjki ,T > ,
			 ExprBjk , typename promotedType>, typename promotedType, i>
operator* (const ExprAjki &a, const ExprBjk &b)
{
	typedef Ajki_to_Aijk < ExprAjki ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bjk < Permuted_Obj , ExprBjk , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i > (ExprObj (Permuted_Obj(a),b));
}

/*!\ brief B(j,k)*A(j,k,i) ==> C(i) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr1 < const Aijk_contracts_Bjk < Ajki_to_Aijk < ExprAjki ,T > ,
			 ExprBjk , typename promotedType>, typename promotedType, i>
operator* (const ExprBjk &b, const ExprAjki &a )
{
	typedef Ajki_to_Aijk < ExprAjki ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bjk < Permuted_Obj , ExprBjk , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i > (ExprObj (Permuted_Obj(a),b));
}

/*!\ brief A(k,j,i)*B(j,k) ==> C(i) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr1 < const Aijk_contracts_Bjk < Akji_to_Aijk < ExprAkji ,T > ,
			 ExprBjk , typename promotedType>, typename promotedType, i>
operator* (const ExprAkji &a, const ExprBjk &b)
{
	typedef Akji_to_Aijk < ExprAkji ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bjk < Permuted_Obj , ExprBjk , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i > (ExprObj (Permuted_Obj(a),b));
}

/*!\ brief B(j,k)*A(k,j,i) ==> C(i) */
//
template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr1 < const Aijk_contracts_Bjk < Akji_to_Aijk < ExprAkji ,T > ,
			 ExprBjk , typename promotedType>, typename promotedType, i>
operator* (const ExprBjk &b, const ExprAkji &a)
{
	typedef Akji_to_Aijk < ExprAkji ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bjk < Permuted_Obj , ExprBjk , typename promotedType> ExprObj;
	return Expr1 < ExprObj,typename promotedType , i > (ExprObj (Permuted_Obj(a),b));
}

//////////////////////////////////////////////////////////////////////////
// Contraction of two indices
//////////////////////////////////////////////////////////////////////////

/*!\ brief A(i,j,k)*B(k,l) ==> C(i,j,l) */
//


template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijk_contracts_Bkl < ExprAijk , ExprBkl , typename promotedType>,
				typename promotedType, i, j, l >
operator* (const ExprAijk &a, const ExprBkl &b)
{
	typedef  const Aijk_contracts_Bkl < ExprAijk , ExprBkl , typename promotedType> ExprObj;
	return Expr3 < ExprObj,typename promotedType , i, j, l> (ExprObj (a,b));
}

/*!\ brief B(k,l)*A(i,j,k) ==> C(i,j,l) */
//


template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijk_contracts_Bkl < ExprAijk , ExprBkl , typename promotedType>,
				typename promotedType, i, j, l >
operator* (const ExprBkl &b,const ExprAijk &a )
{
	typedef  const Aijk_contracts_Bkl < ExprAijk , ExprBkl , typename promotedType> ExprObj;
	return Expr3 < ExprObj,typename promotedType , i, j, l> (ExprObj (a,b));
}


/*!\ brief A(i,j,k)*B(l,k) ==> C(i,j,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijk_contracts_Bkl < ExprAijk , Aji_to_Aij<ExprBlk,U> , typename promotedType>,
				typename promotedType, i, j, l >
operator* (const ExprAijk &a, const ExprBlk &b)
{
	typedef Aji_to_Aij<ExprBlk,U> Permuted_Obj;
	typedef  const Aijk_contracts_Bkl < ExprAijk , Permuted_Obj , typename promotedType> ExprObj;
	return Expr3 < ExprObj,typename promotedType , i, j, l> (ExprObj (a,Permuted_Obj(b)));
}


/*!\ brief B(l,k)*A(i,j,k) ==> C(i,j,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijk_contracts_Bkl < ExprAijk , Aji_to_Aij<ExprBlk,U> , typename promotedType>,
				typename promotedType, i, j, l >
operator* (const ExprBlk &b, const ExprAijk &a)
{
	typedef Aji_to_Aij<ExprBlk,U> Permuted_Obj;
	typedef  const Aijk_contracts_Bkl < ExprAijk , Permuted_Obj , typename promotedType> ExprObj;
	return Expr3 < ExprObj,typename promotedType , i, j, l> (ExprObj (a,Permuted_Obj(b)));
}

/*!\ brief A(i,k,j)*B(k,l) ==> C(i,j,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijk_contracts_Bkl < Aikj_to_Aijk<ExprAikj,T> , ExprBkl , typename promotedType>,
				typename promotedType, i, j, l >
operator* (const ExprAikj &a, const ExprBkl &b)
{
	typedef Aikj_to_Aijk<ExprAikj,T> Permuted_Obj;
	typedef  const Aijk_contracts_Bkl < Permuted_Obj, ExprBkl , typename promotedType> ExprObj;
	return Expr3 < ExprObj,typename promotedType , i, j, l> (ExprObj (Permuted_Obj(a),b));
}

/*!\ brief B(k,l)*A(i,k,j) ==> C(i,j,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijk_contracts_Bkl < Aikj_to_Aijk<ExprAikj,T> , ExprBkl , typename promotedType>,
				typename promotedType, i, j, l >
operator* (const ExprBkl &b, const ExprAikj &a)
{
	typedef Aikj_to_Aijk<ExprAikj,T> Permuted_Obj;
	typedef  const Aijk_contracts_Bkl < Permuted_Obj, ExprBkl , typename promotedType> ExprObj;
	return Expr3 < ExprObj,typename promotedType , i, j, l> (ExprObj (Permuted_Obj(a),b));
}


/*!\ brief A(i,k,j)*B(l,k) ==> C(i,j,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijk_contracts_Bkl < Aikj_to_Aijk<ExprAikj,T> , Aji_to_Aij<ExprBlk,U> , typename promotedType>,
				typename promotedType, i, j, l >
operator* (const ExprAikj &a, const ExprBlk &b)
{
	typedef Aikj_to_Aijk<ExprAikj,T> Permuted_Obj1;
	typedef Aji_to_Aij<ExprBlk,U>  Permuted_Obj2;
	typedef  const Aijk_contracts_Bkl < Permuted_Obj1, Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj,typename promotedType , i, j, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!\ brief B(l,k)*A(i,k,j) ==> C(i,j,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijk_contracts_Bkl < Aikj_to_Aijk<ExprAikj,T> , Aji_to_Aij<ExprBlk,U> , typename promotedType>,
				typename promotedType, i, j, l >
operator* (const ExprBlk &b, const ExprAikj &a)
{
	typedef Aikj_to_Aijk<ExprAikj,T> Permuted_Obj1;
	typedef Aji_to_Aij<ExprBlk,U>  Permuted_Obj2;
	typedef  const Aijk_contracts_Bkl < Permuted_Obj1, Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj,typename promotedType , i, j, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!\ brief A(k,i,j)*B(k,l) ==> C(i,j,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijk_contracts_Bkl < Akij_to_Aijk<ExprAkij,T> , ExprBkl , typename promotedType> ,
				typename promotedType, i, j, l >
operator* (const ExprAkij &a, const ExprBkl &b)
{
	typedef Akij_to_Aijk<ExprAkij,T> Permuted_Obj1;
	typedef  const Aijk_contracts_Bkl < Permuted_Obj1, ExprBkl, typename promotedType> ExprObj;
	return Expr3 < ExprObj,typename promotedType , i, j, l> (ExprObj (Permuted_Obj1(a),b));
}


/*!\ brief B(k,l)*A(k,i,j)==> C(i,j,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijk_contracts_Bkl < Akij_to_Aijk<ExprAkij,T> , ExprBkl ,typename promotedType> ,
				typename promotedType, i, j, l >
operator* (const ExprBkl &b,const ExprAkij &a)
{
	typedef Akij_to_Aijk<ExprAkij,T> Permuted_Obj1;
	typedef  const Aijk_contracts_Bkl < Permuted_Obj1, ExprBkl, typename promotedType> ExprObj;
	return Expr3 < ExprObj,typename promotedType , i, j, l> (ExprObj (Permuted_Obj1(a),b));
}



/*!\ brief A(k,i,j)*B(l,k) ==> C(i,j,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijk_contracts_Bkl < Akij_to_Aijk<ExprAkij,T> , Aji_to_Aij<ExprBlk,U> , typename promotedType>,
				typename promotedType, i, j, l >
operator* (const ExprAkij &a, const ExprBlk &b)
{
	typedef Akij_to_Aijk<ExprAkij,T> Permuted_Obj1;
	typedef Aji_to_Aij<ExprBlk,U> Permuted_Obj2;
	typedef  const Aijk_contracts_Bkl < Permuted_Obj1, Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj,typename promotedType , i, j, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!\ brief B(l,k)*A(k,i,j)==> C(i,j,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr3 < const Aijk_contracts_Bkl < Akij_to_Aijk<ExprAkij,T> , Aji_to_Aij<ExprBlk,U> , typename promotedType>,
				typename promotedType, i, j, l >
operator* (const ExprBlk &b,const ExprAkij &a)
{
	typedef Akij_to_Aijk<ExprAkij,T> Permuted_Obj1;
	typedef Aji_to_Aij<ExprBlk,U> Permuted_Obj2;
	typedef  const Aijk_contracts_Bkl < Permuted_Obj1, Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj,typename promotedType , i, j, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


#undef ExprAijk
#undef ExprAikj
#undef ExprAjik
#undef ExprAkij
#undef ExprAjki
#undef ExprAkji

#undef ExprBjk
#undef promotedType

#endif
