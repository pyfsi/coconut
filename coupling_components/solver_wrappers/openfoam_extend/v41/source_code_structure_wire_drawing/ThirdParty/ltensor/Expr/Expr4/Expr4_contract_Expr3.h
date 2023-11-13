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
#ifndef Expr4_contract_Expr3_H
#define Expr4_contract_Expr3_H


////////////////////////////////////////////
/*!brief A(i,j,k,l) * B(k,l,m) -> C(i,j,m) */
///////////////////////////////////////////

#define ExprAijkl	Expr4<A,T,i,j,k,l>
#define ExprAijlk	Expr4<A,T,i,j,l,k>
#define ExprAklij	Expr4<A,T,k,l,i,j>
#define ExprAlkij	Expr4<A,T,l,k,i,j>
#define ExprAikjl	Expr4<A,T,i,k,j,l>
#define ExprAjikl	Expr4<A,T,j,i,k,l>
#define ExprAiklj   Expr4<A,T,i,k,l,j>
#define ExprAkijl   Expr4<A,T,k,i,j,l>
#define ExprAkilj   Expr4<A,T,k,i,l,j>
#define ExprAjkil	Expr4<A,T,j,k,i,l>
#define ExprAjkli	Expr4<A,T,j,k,l,i>

#define ExprBklm	Expr3<B,U,k,l,m>
#define ExprBmkl	Expr3<B,U,m,k,l>
#define ExprBmlk	Expr3<B,U,m,l,k>
#define ExprBkml	Expr3<B,U,k,m,l>
#define ExprBlkm	Expr3<B,U,l,k,m>
#define ExprBlmk	Expr3<B,U,l,m,k>
#define ExprBjkl	Expr3<B,U,j,k,l>
#define ExprBjlk	Expr3<B,U,j,l,k>
#define ExprBkjl	Expr3<B,U,k,j,l>
#define ExprBklj	Expr3<B,U,k,l,j>
#define ExprBljk	Expr3<B,U,l,j,k>
#define ExprBlkj	Expr3<B,U,l,k,j>

#define promotedType	promote<T,U>::V

/////////////////////////////////////////////////////////////////////
// Double Contraction with first 2 indices of Third Order Tensor
////////////////////////////////////////////////////////////////////

/*!brief A(i,j,k,l)*B(k,l,m) ==> C(i,j,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm < ExprAijkl , ExprBklm , typename promotedType>,
				typename promotedType, i , j, m>
operator* (const ExprAijkl &a, const ExprBklm &b)
{
	typedef  const Aijkl_contracts_Bklm < ExprAijkl , ExprBklm , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType , i, j, m > (ExprObj (a,b));
}//CHECKED 1

/*!brief B(k,l,m)*A(i,j,k,l) ==> C(i,j,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm < ExprAijkl , ExprBklm , typename promotedType>,
				typename promotedType, i , j, m>
operator* ( const ExprBklm &b, const ExprAijkl &a)
{
	typedef  const Aijkl_contracts_Bklm < ExprAijkl , ExprBklm , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType , i, j, m > (ExprObj (a,b));
}


/*!brief A(i,j,l,k)*B(k,l,m) ==> C(i,j,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm < Aijlk_to_Aijkl < ExprAijlk ,T> ,
			 ExprBklm , typename promotedType>, typename promotedType, i, j, m >
operator* (const ExprAijlk &a, const ExprBklm &b)
{
	typedef Aijlk_to_Aijkl < ExprAijlk ,T > Permuted_Obj;
	typedef const Aijkl_contracts_Bklm < Permuted_Obj , ExprBklm , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType , i, j, m > (ExprObj (Permuted_Obj(a),b));
}//CHECKED 2

/*!brief B(k,l,m)*A(i,j,l,k) ==> C(i,j,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm < Aijlk_to_Aijkl < ExprAijlk ,T> ,
			 ExprBklm , typename promotedType>, typename promotedType, i, j, m >
operator* (const ExprBklm &b, const ExprAijlk &a)
{
	typedef Aijlk_to_Aijkl < ExprAijlk ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Bklm < Permuted_Obj , ExprBklm , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType , i, j, m > (ExprObj (Permuted_Obj(a),b));
}

/*!brief A(k,l,i,j)*B(k,l,m) ==> C(i,j,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm <
				Aklij_to_Aijkl < ExprAklij ,T>,
			 	ExprBklm,
				typename promotedType>,
		 typename promotedType, i, j, m>
operator* (const ExprAklij &a, const ExprBklm &b)
{
	typedef Aklij_to_Aijkl < ExprAklij ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Bklm < Permuted_Obj ,  ExprBklm , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType , i, j, m > (ExprObj (Permuted_Obj(a),b));
}//CHECKED 3

/*!brief B(k,l,m)*A(k,l,i,j) ==> C(i,j,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm <
				Aklij_to_Aijkl < ExprAklij ,T>,
			 	ExprBklm,
				typename promotedType>,
		 typename promotedType, i, j, m>
operator* (const ExprBklm &b, const ExprAklij &a)
{
	typedef Aklij_to_Aijkl < ExprAklij ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Bklm < Permuted_Obj ,  ExprBklm , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType , i, j, m > (ExprObj (Permuted_Obj(a),b));
}

/*!brief A(l,k,i,j)*B(k,l,m) ==> C(i,j,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm <
				Alkij_to_Aijkl < ExprAlkij ,T>,
			 	ExprBklm,
				typename promotedType>,
		 typename promotedType, i, j, m>
operator* (const ExprAlkij &a, const ExprBklm &b)
{
	typedef Alkij_to_Aijkl < ExprAlkij ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Bklm < Permuted_Obj ,  ExprBklm , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj(a),b));
}//CHECKED 4

/*!brief B(k,l,m)*A(l,k,i,j) ==> C(i,j,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm <
				Alkij_to_Aijkl < ExprAlkij ,T>,
			 	ExprBklm,
				typename promotedType>,
		 typename promotedType, i, j, m>
operator* (const ExprBklm &b, const ExprAlkij &a)
{
	typedef Alkij_to_Aijkl < ExprAlkij ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Bklm < Permuted_Obj ,  ExprBklm , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj(a),b));
}

/////////////////////////////////////////////////////////////////////
// Double Contraction with last 2 indices of Third Order Tensor
////////////////////////////////////////////////////////////////////


/*!brief A(i,j,k,l)*B(m,k,l) ==> C(i,j,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm < ExprAijkl,
			 	Akij_to_Aijk < ExprBmkl, U>,
				typename promotedType>, typename promotedType, i, j, m>
operator* (const ExprAijkl &a, const ExprBmkl &b)
{
	typedef Akij_to_Aijk < ExprBmkl ,U> Permuted_Obj1;
	typedef const Aijkl_contracts_Bklm < ExprAijkl , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m> (ExprObj (a,Permuted_Obj1(b)));
}//CHECKED 5

/*!brief B(m,k,l)*A(i,j,k,l) ==> C(i,j,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm < ExprAijkl,
			 	Akij_to_Aijk < ExprBmkl, U>,
				typename promotedType>, typename promotedType, i, j, m>
operator* (const ExprBmkl &b, const ExprAijkl &a)
{
	typedef Akij_to_Aijk < ExprBmkl ,U> Permuted_Obj1;
	typedef const Aijkl_contracts_Bklm < ExprAijkl , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m> (ExprObj (a,Permuted_Obj1(b)));
}

/*!brief A(i,j,k,l)*B(m,l,k) ==> C(i,j,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm < ExprAijkl,
			 	Akji_to_Aijk < ExprBmlk, U>,
				typename promotedType>, typename promotedType, i, j, m>
operator* (const ExprAijkl &a, const ExprBmlk &b)
{
	typedef Akji_to_Aijk < ExprBmlk ,U> Permuted_Obj1;
	typedef const Aijkl_contracts_Bklm < ExprAijkl , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m> (ExprObj (a,Permuted_Obj1(b)));
}//CHECKED 6

/*!brief B(m,l,k)*A(i,j,k,l) ==> C(i,j,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm < ExprAijkl,
			 	Akji_to_Aijk < ExprBmlk, U>,
				typename promotedType>, typename promotedType, i, j, m>

operator* (const ExprBmlk &b, const ExprAijkl &a)
{
	typedef Akji_to_Aijk < ExprBmlk ,U> Permuted_Obj1;
	typedef const Aijkl_contracts_Bklm < ExprAijkl , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m> (ExprObj (a,Permuted_Obj1(b)));
}


/*!brief A(k,l,i,j)*B(m,k,l) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm <
				Aklij_to_Aijkl < ExprAklij ,T>,
			 	Akij_to_Aijk < ExprBmkl, U>,
				typename promotedType>,
		 		typename promotedType, i, j, m>
operator* (const ExprAklij &a, const ExprBmkl &b)
{
	typedef Aklij_to_Aijkl < ExprAklij ,T> Permuted_Obj;
	typedef Akij_to_Aijk < ExprBmkl ,U> Permuted_Obj1;
	typedef const Aijkl_contracts_Bklm < Permuted_Obj , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj(a),Permuted_Obj1(b)));
}//CHECKED 7

/*!brief B(m,k,l)*A(k,l,i,j) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm <
				Aklij_to_Aijkl < ExprAklij ,T>,
			 	Akij_to_Aijk < ExprBmkl, U>,
				typename promotedType>,
		 		typename promotedType, i, j, m>
operator* (const ExprBmkl &b, const ExprAklij &a)
{
	typedef Aklij_to_Aijkl < ExprAklij ,T> Permuted_Obj;
	typedef Akij_to_Aijk < ExprBmkl ,U> Permuted_Obj1;
	typedef const Aijkl_contracts_Bklm < Permuted_Obj , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj(a),Permuted_Obj1(b)));
}

/*!brief A(k,l,i,j)*B(m,l,k) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm <
				Aklij_to_Aijkl < ExprAklij ,T>,
			 	Akji_to_Aijk < ExprBmlk, U>,
				typename promotedType>,
		 		typename promotedType, i, j, m>
operator* (const ExprAklij &a, const ExprBmlk &b)
{
	typedef Aklij_to_Aijkl < ExprAklij ,T> Permuted_Obj;
	typedef Akji_to_Aijk < ExprBmlk ,U > Permuted_Obj1;
	typedef const Aijkl_contracts_Bklm < Permuted_Obj , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj(a),Permuted_Obj1(b)));
}//CHECKED 8

/*!brief B(m,l,k)*A(k,l,i,j) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3 < const Aijkl_contracts_Bklm <
				Aklij_to_Aijkl < ExprAklij ,T>,
			 	Akji_to_Aijk < ExprBmlk, U>,
				typename promotedType>,
		 		typename promotedType, i, j, m>
operator* (const ExprBmlk &b, const ExprAklij &a)
{
	typedef Aklij_to_Aijkl < ExprAklij ,T> Permuted_Obj;
	typedef Akji_to_Aijk < ExprBmlk ,U > Permuted_Obj1;
	typedef const Aijkl_contracts_Bklm < Permuted_Obj , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj(a),Permuted_Obj1(b)));
}


/*!brief A(i,j,k,l)*B(k,m,l) */
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            ExprAijkl,
            Aikj_to_Aijk<ExprBkml,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAijkl &a,const ExprBkml &b){
    typedef Aikj_to_Aijk<ExprBkml,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklm < ExprAijkl , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (a,Permuted_Obj1(b)));
}//CHECKED 9

/*!brief B(k,m,l)*A(i,j,k,l) */
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aikj_to_Aijk<ExprBkml,U>,
            ExprAijkl,
            typename promotedType>,
            typename promotedType, i, j,m >
operator*(const ExprBkml &b , const ExprAijkl &a){
    typedef Aikj_to_Aijk<ExprBkml,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklm < ExprAijkl , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i , j , m > (ExprObj (a,Permuted_Obj1(b)));
}


/*!brief A(i,j,k,l)*B(l,m,k) */
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            ExprAijkl,
            Ajki_to_Aijk<ExprBlmk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAijkl &a,const ExprBlmk &b){
    typedef Ajki_to_Aijk<ExprBlmk,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklm < ExprAijkl , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (a,Permuted_Obj1(b)));
}//CHECKED 10

/*!brief B(l,m,k)*A(i,j,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            ExprAijkl,
            Ajki_to_Aijk<ExprBlmk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBlmk &b,const ExprAijkl &a){
    typedef Ajki_to_Aijk<ExprBlmk,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklm < ExprAijkl , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (a,Permuted_Obj1(b)));
}
/*!brief A(i,k,j,l)*B(k,l,m) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aikjl_to_Aijkl<ExprAikjl,T>,
            ExprBklm,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAikjl &a,const ExprBklm &b){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , ExprBklm , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),b));
}//CHECKED 11

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aikjl_to_Aijkl<ExprAikjl,T>,
            ExprBklm,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBklm &b,const ExprAikjl &a){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , ExprBklm , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),b));
}

/*!brief A(i,k,j,l)*B(k,m,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aikjl_to_Aijkl<ExprAikjl,T>,
            Aikj_to_Aijk<ExprBkml,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAikjl &a,const ExprBkml &b){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Aikj_to_Aijk<ExprBkml,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 12


template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aikjl_to_Aijkl<ExprAikjl,T>,
            Aikj_to_Aijk<ExprBkml,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBkml &b,const ExprAikjl &a){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Aikj_to_Aijk<ExprBkml,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}

/*!brief A(i,k,j,l)*B(l,k,m) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aikjl_to_Aijkl<ExprAikjl,T>,
            Ajik_to_Aijk<ExprBlkm,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAikjl &a,const ExprBlkm &b){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Ajik_to_Aijk<ExprBlkm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 13



template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aikjl_to_Aijkl<ExprAikjl,T>,
            Ajik_to_Aijk<ExprBlkm,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBlkm &b,const ExprAikjl &a){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Ajik_to_Aijk<ExprBlkm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(i,k,j,l)*B(l,m,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aikjl_to_Aijkl<ExprAikjl,T>,
            Ajki_to_Aijk<ExprBlmk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAikjl &a,const ExprBlmk &b){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Ajki_to_Aijk<ExprBlmk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 14

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aikjl_to_Aijkl<ExprAikjl,T>,
            Ajki_to_Aijk<ExprBlmk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBlmk &b,const ExprAikjl &a){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Ajki_to_Aijk<ExprBlmk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(i,k,j,l)*B(m,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aikjl_to_Aijkl<ExprAikjl,T>,
            Akij_to_Aijk<ExprBmkl,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAikjl &a,const ExprBmkl &b){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Akij_to_Aijk<ExprBmkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 15

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aikjl_to_Aijkl<ExprAikjl,T>,
            Akij_to_Aijk<ExprBmkl,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBmkl &b,const ExprAikjl &a){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Akij_to_Aijk<ExprBmkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief A(i,k,j,l)*B(m,l,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aikjl_to_Aijkl<ExprAikjl,T>,
            Akji_to_Aijk<ExprBmlk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAikjl &a,const ExprBmlk &b){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Akji_to_Aijk<ExprBmlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 16

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aikjl_to_Aijkl<ExprAikjl,T>,
            Akji_to_Aijk<ExprBmlk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBmlk &b,const ExprAikjl &a){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Akji_to_Aijk<ExprBmlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief A(i,k,l,j)*B(k,l,m) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aiklj_to_Aijkl<ExprAiklj,T>,
            ExprBklm,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAiklj &a,const ExprBklm &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , ExprBklm , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),b));
}//CHECKED 17

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aiklj_to_Aijkl<ExprAiklj,T>,
            ExprBklm,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBklm &b,const ExprAiklj &a){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , ExprBklm , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),b));
}

/*!brief A(i,k,l,j)*B(k,m,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aiklj_to_Aijkl<ExprAiklj,T>,
            Aikj_to_Aijk<ExprBkml,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAiklj &a,const ExprBkml &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Aikj_to_Aijk<ExprBkml,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 18

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aiklj_to_Aijkl<ExprAiklj,T>,
            Aikj_to_Aijk<ExprBkml,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBkml &b,const ExprAiklj &a){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Aikj_to_Aijk<ExprBkml,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief A(i,k,l,j)*B(l,k,m) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aiklj_to_Aijkl<ExprAiklj,T>,
            Ajik_to_Aijk<ExprBlkm,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAiklj &a,const ExprBlkm &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Ajik_to_Aijk<ExprBlkm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 19

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aiklj_to_Aijkl<ExprAiklj,T>,
            Ajik_to_Aijk<ExprBlkm,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBlkm &b,const ExprAiklj &a){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Ajik_to_Aijk<ExprBlkm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}

/*!brief A(i,k,l,j)*B(l,m,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aiklj_to_Aijkl<ExprAiklj,T>,
            Ajki_to_Aijk<ExprBlmk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAiklj &a,const ExprBlmk &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Ajki_to_Aijk<ExprBlmk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 20

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aiklj_to_Aijkl<ExprAiklj,T>,
            Ajki_to_Aijk<ExprBlmk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBlmk &b,const ExprAiklj &a){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Ajki_to_Aijk<ExprBlmk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief A(i,k,l,j)*B(m,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aiklj_to_Aijkl<ExprAiklj,T>,
            Akij_to_Aijk<ExprBmkl,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAiklj &a,const ExprBmkl &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Akij_to_Aijk<ExprBmkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 21

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aiklj_to_Aijkl<ExprAiklj,T>,
            Akij_to_Aijk<ExprBmkl,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBmkl &b,const ExprAiklj &a ){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Akij_to_Aijk<ExprBmkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief A(i,k,l,j)*B(m,l,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aiklj_to_Aijkl<ExprAiklj,T>,
            Akji_to_Aijk<ExprBmlk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAiklj &a,const ExprBmlk &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Akji_to_Aijk<ExprBmlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 22


template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aiklj_to_Aijkl<ExprAiklj,T>,
            Akji_to_Aijk<ExprBmlk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBmlk &b,const ExprAiklj &a){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Akji_to_Aijk<ExprBmlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}
/*!brief A(k,i,j,l)*B(k,l,m) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akijl_to_Aijkl<ExprAkijl,T>,
            ExprBklm,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAkijl &a,const ExprBklm &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , ExprBklm , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),b));
}//CHECKED 24

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akijl_to_Aijkl<ExprAkijl,T>,
            ExprBklm,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBklm &b,const ExprAkijl &a){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , ExprBklm , typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),b));
}



/*!brief A(k,i,j,l)*B(k,m,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akijl_to_Aijkl<ExprAkijl,T>,
            Aikj_to_Aijk<ExprBkml,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAkijl &a,const ExprBkml &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Aikj_to_Aijk<ExprBkml, U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 25

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akijl_to_Aijkl<ExprAkijl,T>,
            Aikj_to_Aijk<ExprBkml,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBkml &b,const ExprAkijl &a){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Aikj_to_Aijk<ExprBkml, U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief A(k,i,j,l)*B(l,k,m) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akijl_to_Aijkl<ExprAkijl,T>,
            Ajik_to_Aijk<ExprBlkm,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAkijl &a,const ExprBlkm &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Ajik_to_Aijk<ExprBlkm, U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 26


template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akijl_to_Aijkl<ExprAkijl,T>,
            Ajik_to_Aijk<ExprBlkm,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBlkm &b,const ExprAkijl &a){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Ajik_to_Aijk<ExprBlkm, U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief A(k,i,j,l)*B(l,m,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akijl_to_Aijkl<ExprAkijl,T>,
            Ajki_to_Aijk<ExprBlmk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAkijl &a,const ExprBlmk &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Ajki_to_Aijk<ExprBlmk, U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 27

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akijl_to_Aijkl<ExprAkijl,T>,
            Ajki_to_Aijk<ExprBlmk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBlmk &b,const ExprAkijl &a){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Ajki_to_Aijk<ExprBlmk, U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief A(k,i,j,l)*B(m,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akijl_to_Aijkl<ExprAkijl,T>,
            Akij_to_Aijk<ExprBmkl,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAkijl &a,const ExprBmkl &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Akij_to_Aijk<ExprBmkl, U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 28


template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akijl_to_Aijkl<ExprAkijl,T>,
            Akij_to_Aijk<ExprBmkl,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBmkl &b,const ExprAkijl &a){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Akij_to_Aijk<ExprBmkl, U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(k,i,j,l)*B(m,l,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akijl_to_Aijkl<ExprAkijl,T>,
            Akji_to_Aijk<ExprBmlk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAkijl &a,const ExprBmlk &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Akji_to_Aijk<ExprBmlk, U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 29

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akijl_to_Aijkl<ExprAkijl,T>,
            Akji_to_Aijk<ExprBmlk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBmlk &b,const ExprAkijl &a){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Akji_to_Aijk<ExprBmlk, U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief A(k,i,l,j)*B(k,l,m) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akilj_to_Aijkl<ExprAkilj,T>,
            ExprBklm,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAkilj &a,const ExprBklm &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , ExprBklm, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),b));
}//CHECKED 30


template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akilj_to_Aijkl<ExprAkilj,T>,
            ExprBklm,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBklm &b,const ExprAkilj &a){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , ExprBklm, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),b));
}



/*!brief A(k,i,l,j)*B(l,k,m) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akilj_to_Aijkl<ExprAkilj,T>,
            Ajik_to_Aijk<ExprBlkm,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAkilj &a,const ExprBlkm &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Ajik_to_Aijk<ExprBlkm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED



/*!brief B(l,k,m)*A(k,i,l,j) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akilj_to_Aijkl<ExprAkilj,T>,
            Ajik_to_Aijk<ExprBlkm,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBlkm &b,const ExprAkilj &a){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Ajik_to_Aijk<ExprBlkm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED RIGHT?


/*!brief A(k,i,l,j)*B(k,m,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akilj_to_Aijkl<ExprAkilj,T>,
            Aikj_to_Aijk<ExprBkml,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAkilj &a,const ExprBkml &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Aikj_to_Aijk<ExprBkml,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 31

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akilj_to_Aijkl<ExprAkilj,T>,
            Aikj_to_Aijk<ExprBkml,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBkml &b,const ExprAkilj &a){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Aikj_to_Aijk<ExprBkml,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}

/*!brief A(k,i,l,j)*B(l,m,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akilj_to_Aijkl<ExprAkilj,T>,
            Ajki_to_Aijk<ExprBlmk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAkilj &a,const ExprBlmk &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Ajki_to_Aijk<ExprBlmk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 32


template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akilj_to_Aijkl<ExprAkilj,T>,
            Ajki_to_Aijk<ExprBlmk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBlmk &b,const ExprAkilj &a){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Ajki_to_Aijk<ExprBlmk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}




/*!brief A(k,i,l,j)*B(m,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akilj_to_Aijkl<ExprAkilj,T>,
            Akij_to_Aijk<ExprBmkl,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAkilj &a,const ExprBmkl &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Akij_to_Aijk<ExprBmkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 33

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akilj_to_Aijkl<ExprAkilj,T>,
            Akij_to_Aijk<ExprBmkl,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBmkl &b,const ExprAkilj &a){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Akij_to_Aijk<ExprBmkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}

/*!brief A(k,i,l,j)*B(m,l,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akilj_to_Aijkl<ExprAkilj,T>,
            Akji_to_Aijk<ExprBmlk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAkilj &a,const ExprBmlk &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Akji_to_Aijk<ExprBmlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 34


template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Akilj_to_Aijkl<ExprAkilj,T>,
            Akji_to_Aijk<ExprBmlk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBmlk &b,const ExprAkilj &a){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Akji_to_Aijk<ExprBmlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(k,l,i,j)*B(k,m,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aklij_to_Aijkl<ExprAklij,T>,
            Aikj_to_Aijk<ExprBkml,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAklij &a,const ExprBkml &b){
    typedef Aklij_to_Aijkl<ExprAklij,T> Permuted_Obj1;
    typedef Aikj_to_Aijk<ExprBkml,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 35

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aklij_to_Aijkl<ExprAklij,T>,
            Aikj_to_Aijk<ExprBkml,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBkml &b, const ExprAklij &a){
    typedef Aklij_to_Aijkl<ExprAklij,T> Permuted_Obj1;
    typedef Aikj_to_Aijk<ExprBkml,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}





/*!brief A(k,l,i,j)*B(l,m,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aklij_to_Aijkl<ExprAklij,T>,
            Ajki_to_Aijk<ExprBlmk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprAklij &a,const ExprBlmk &b){
    typedef Aklij_to_Aijkl<ExprAklij,T> Permuted_Obj1;
    typedef Ajki_to_Aijk<ExprBlmk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 37

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr3<const Aijkl_contracts_Bklm<
            Aklij_to_Aijkl<ExprAklij,T>,
            Ajki_to_Aijk<ExprBlmk,U>,
            typename promotedType>,
            typename promotedType,i,j,m>
operator*(const ExprBlmk &b,const ExprAklij &a){
    typedef Aklij_to_Aijkl<ExprAklij,T> Permuted_Obj1;
    typedef Ajki_to_Aijk<ExprBlmk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklm < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr3 < ExprObj, typename promotedType, i, j, m > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


////////////////////////////
//Three contracting indexes
//makes sense right?
///////////////////////////


/*!brief A(i,j,k,l)*B(j,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            ExprAijkl,
            ExprBjkl,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAijkl &a,const ExprBjkl &b){


    typedef const Aijkl_contracts_Bjkl < ExprAijkl , ExprBjkl, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (a,b));
}


/*!brief B(j,k,l)*A(i,j,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            ExprAijkl,
            ExprBjkl,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBjkl &b, const ExprAijkl &a){


    typedef const Aijkl_contracts_Bjkl < ExprAijkl , ExprBjkl, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (a,b));
}


/*!brief A(i,j,k,l)*B(j,l,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            ExprAijkl,
            Aikj_to_Aijk<ExprBjlk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAijkl &a,const ExprBjlk &b){

    typedef Aikj_to_Aijk<ExprBjlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjkl < ExprAijkl , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (a,Permuted_Obj2(b)));
}



/*!brief B(j,l,k)*A(i,j,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            ExprAijkl,
            Aikj_to_Aijk<ExprBjlk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBjlk &b,const ExprAijkl &a){

    typedef Aikj_to_Aijk<ExprBjlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjkl < ExprAijkl , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (a,Permuted_Obj2(b)));
}


/*!brief A(i,j,k,l)*B(k,j,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            ExprAijkl,
            Ajik_to_Aijk<ExprBkjl,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAijkl &a,const ExprBkjl &b){

    typedef Ajik_to_Aijk<ExprBkjl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjkl < ExprAijkl , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (a,Permuted_Obj2(b)));
}



/*!brief B(k,j,l)*A(i,j,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            ExprAijkl,
            Ajik_to_Aijk<ExprBkjl,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBkjl &b, const ExprAijkl &a){

    typedef Ajik_to_Aijk<ExprBkjl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjkl < ExprAijkl , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (a,Permuted_Obj2(b)));
}




/*!brief A(i,j,k,l)*B(k,l,j) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            ExprAijkl,
            Ajki_to_Aijk<ExprBklj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAijkl &a,const ExprBklj &b){

    typedef Ajki_to_Aijk<ExprBklj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjkl < ExprAijkl , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (a,Permuted_Obj2(b)));
}


/*!brief B(k,l,j)*A(i,j,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            ExprAijkl,
            Ajki_to_Aijk<ExprBklj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBklj &b, const ExprAijkl &a){

    typedef Ajki_to_Aijk<ExprBklj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjkl < ExprAijkl , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (a,Permuted_Obj2(b)));
}



/*!brief A(i,j,k,l)*B(l,j,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            ExprAijkl,
            Akij_to_Aijk<ExprBljk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAijkl &a,const ExprBljk &b){

    typedef Akij_to_Aijk<ExprBljk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjkl < ExprAijkl , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (a,Permuted_Obj2(b)));
}



/*!brief B(l,j,k)*A(i,j,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            ExprAijkl,
            Akij_to_Aijk<ExprBljk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBljk &b, const ExprAijkl &a){

    typedef Akij_to_Aijk<ExprBljk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjkl < ExprAijkl , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (a,Permuted_Obj2(b)));
}



/*!brief A(i,j,k,l)*B(l,k,j) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            ExprAijkl,
            Akji_to_Aijk<ExprBlkj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAijkl &a,const ExprBlkj &b){

    typedef Akji_to_Aijk<ExprBlkj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjkl < ExprAijkl , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (a,Permuted_Obj2(b)));
}



/*!brief B(l,k,j)* A(i,j,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            ExprAijkl,
            Akji_to_Aijk<ExprBlkj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBlkj &b, const ExprAijkl &a){

    typedef Akji_to_Aijk<ExprBlkj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjkl < ExprAijkl , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (a,Permuted_Obj2(b)));
}


/*!brief A(j,i,k,l)*B(j,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajikl_to_Aijkl < ExprAjikl ,T> ,
            ExprBjkl ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjikl &a,const ExprBjkl &b){

	typedef Ajikl_to_Aijkl < ExprAjikl ,T> Permuted_Obj1;

    typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , ExprBjkl, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),b));
}



/*!brief B(j,k,l) * A(j,i,k,l)*/

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajikl_to_Aijkl < ExprAjikl ,T> ,
            ExprBjkl ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBjkl &b, const ExprAjikl &a){

	typedef Ajikl_to_Aijkl < ExprAjikl ,T> Permuted_Obj1;

    typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , ExprBjkl, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),b));
}




/*!brief A(j,i,k,l)*B(j,l,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajikl_to_Aijkl < ExprAjikl ,T> ,
            Aikj_to_Aijk< ExprBjlk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjikl &a,const ExprBjlk &b){

	typedef Ajikl_to_Aijkl < ExprAjikl ,T> Permuted_Obj1;
	typedef Aikj_to_Aijk< ExprBjlk,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief B(j,l,k) * A(j,i,k,l)*/

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajikl_to_Aijkl < ExprAjikl ,T> ,
            Aikj_to_Aijk< ExprBjlk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBjlk &b, const ExprAjikl &a){

	typedef Ajikl_to_Aijkl < ExprAjikl ,T> Permuted_Obj1;
	typedef Aikj_to_Aijk< ExprBjlk,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(j,i,k,l)*B(k,j,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajikl_to_Aijkl < ExprAjikl ,T> ,
            Ajik_to_Aijk< ExprBkjl,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjikl &a,const ExprBkjl &b){

	typedef Ajikl_to_Aijkl < ExprAjikl ,T> Permuted_Obj1;
	typedef Ajik_to_Aijk< ExprBkjl,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief B(k,j,l)*A(j,i,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajikl_to_Aijkl < ExprAjikl ,T> ,
            Ajik_to_Aijk< ExprBkjl,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBkjl &b, const ExprAjikl &a){

	typedef Ajikl_to_Aijkl < ExprAjikl ,T> Permuted_Obj1;
	typedef Ajik_to_Aijk< ExprBkjl,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(j,i,k,l)*B(k,l,j) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajikl_to_Aijkl < ExprAjikl ,T> ,
            Ajki_to_Aijk< ExprBklj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjikl &a,const ExprBklj &b){

	typedef Ajikl_to_Aijkl < ExprAjikl ,T> Permuted_Obj1;
	typedef Ajki_to_Aijk< ExprBklj,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief B(k,l,j)*A(j,i,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajikl_to_Aijkl < ExprAjikl ,T> ,
            Ajki_to_Aijk< ExprBklj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBklj &b, const ExprAjikl &a){

	typedef Ajikl_to_Aijkl < ExprAjikl ,T> Permuted_Obj1;
	typedef Ajki_to_Aijk< ExprBklj,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief A(j,i,k,l)*B(l,j,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajikl_to_Aijkl < ExprAjikl ,T> ,
            Akij_to_Aijk< ExprBljk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjikl &a,const ExprBljk &b){

	typedef Ajikl_to_Aijkl < ExprAjikl ,T> Permuted_Obj1;
	typedef Akij_to_Aijk< ExprBljk,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}

/*!brief B(l,j,k)* A(j,i,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajikl_to_Aijkl < ExprAjikl ,T> ,
            Akij_to_Aijk< ExprBljk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBljk &b, const ExprAjikl &a){

	typedef Ajikl_to_Aijkl < ExprAjikl ,T> Permuted_Obj1;
	typedef Akij_to_Aijk< ExprBljk,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(j,i,k,l)*B(l,k,j) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajikl_to_Aijkl < ExprAjikl ,T> ,
            Akji_to_Aijk< ExprBlkj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjikl &a,const ExprBlkj &b){

	typedef Ajikl_to_Aijkl < ExprAjikl ,T> Permuted_Obj1;
	typedef Akji_to_Aijk< ExprBlkj,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief B(l,k,j)*A(j,i,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajikl_to_Aijkl < ExprAjikl ,T> ,
            Akji_to_Aijk< ExprBlkj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBlkj &b, const ExprAjikl &a){

	typedef Ajikl_to_Aijkl < ExprAjikl ,T> Permuted_Obj1;
	typedef Akji_to_Aijk< ExprBlkj,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(j,k,i,l)*B(j,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkil_to_Aijkl < ExprAjkil ,T> ,
            ExprBjkl ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjkil &a,const ExprBjkl &b){

	typedef Ajkil_to_Aijkl < ExprAjkil ,T> Permuted_Obj1;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , ExprBjkl, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),b));
}




/*!brief B(j,k,l) * A(j,k,i,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkil_to_Aijkl < ExprAjkil ,T> ,
            ExprBjkl ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBjkl &b, const ExprAjkil &a){

	typedef Ajkil_to_Aijkl < ExprAjkil ,T> Permuted_Obj1;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , ExprBjkl, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),b));
}


/*!brief A(j,k,i,l)*B(j,l,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkil_to_Aijkl < ExprAjkil ,T> ,
            Aikj_to_Aijk< ExprBjlk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjkil &a,const ExprBjlk &b){

	typedef Ajkil_to_Aijkl < ExprAjkil ,T> Permuted_Obj1;
	typedef Aikj_to_Aijk< ExprBjlk,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief B(j,l,k) * A(j,k,i,l)*/

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkil_to_Aijkl < ExprAjkil ,T> ,
            Aikj_to_Aijk< ExprBjlk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBjlk &b, const ExprAjkil &a){

	typedef Ajkil_to_Aijkl < ExprAjkil ,T> Permuted_Obj1;
	typedef Aikj_to_Aijk< ExprBjlk,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(j,k,i,l)*B(k,j,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkil_to_Aijkl < ExprAjkil ,T> ,
            Ajik_to_Aijk< ExprBkjl,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjkil &a,const ExprBkjl &b){

	typedef Ajkil_to_Aijkl < ExprAjkil ,T> Permuted_Obj1;
	typedef Ajik_to_Aijk< ExprBkjl,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief B(k,j,l)*A(j,k,i,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkil_to_Aijkl < ExprAjkil ,T> ,
            Ajik_to_Aijk< ExprBkjl,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBkjl &b, const ExprAjkil &a){

	typedef Ajkil_to_Aijkl < ExprAjkil ,T> Permuted_Obj1;
	typedef Ajik_to_Aijk< ExprBkjl,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(j,k,i,l)*B(k,l,j) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkil_to_Aijkl < ExprAjkil ,T> ,
            Ajki_to_Aijk< ExprBklj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjkil &a,const ExprBklj &b){

	typedef Ajkil_to_Aijkl < ExprAjkil ,T> Permuted_Obj1;
	typedef Ajki_to_Aijk< ExprBklj,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief B(k,l,j)*A(j,k,i,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkil_to_Aijkl < ExprAjkil ,T> ,
            Ajki_to_Aijk< ExprBklj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBklj &b, const ExprAjkil &a){

	typedef Ajkil_to_Aijkl < ExprAjkil ,T> Permuted_Obj1;
	typedef Ajki_to_Aijk< ExprBklj,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(j,k,i,l)*B(l,j,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkil_to_Aijkl < ExprAjkil ,T> ,
            Akij_to_Aijk< ExprBljk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjkil &a,const ExprBljk &b){

	typedef Ajkil_to_Aijkl < ExprAjkil ,T> Permuted_Obj1;
	typedef Akij_to_Aijk< ExprBljk,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief B(l,j,k)*A(j,k,i,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkil_to_Aijkl < ExprAjkil ,T> ,
            Akij_to_Aijk< ExprBljk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBljk &b, const ExprAjkil &a){

	typedef Ajkil_to_Aijkl < ExprAjkil ,T> Permuted_Obj1;
	typedef Akij_to_Aijk< ExprBljk,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(j,k,i,l)*B(l,k,j) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkil_to_Aijkl < ExprAjkil ,T> ,
            Akji_to_Aijk< ExprBlkj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjkil &a,const ExprBlkj &b){

	typedef Ajkil_to_Aijkl < ExprAjkil ,T> Permuted_Obj1;
	typedef Akji_to_Aijk< ExprBlkj,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief B(l,k,j)*A(j,k,i,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkil_to_Aijkl < ExprAjkil ,T> ,
            Akji_to_Aijk< ExprBlkj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBlkj &b, const ExprAjkil &a){

	typedef Ajkil_to_Aijkl < ExprAjkil ,T> Permuted_Obj1;
	typedef Akji_to_Aijk< ExprBlkj,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}





/*!brief A(j,k,l,i)*B(j,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkli_to_Aijkl < ExprAjkli ,T> ,
            ExprBjkl ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjkli &a,const ExprBjkl &b){

	typedef Ajkli_to_Aijkl < ExprAjkli ,T> Permuted_Obj1;


	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , ExprBjkl, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),b));
}


/*!brief B(j,k,l)*A(j,k,l,i) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkli_to_Aijkl < ExprAjkli ,T> ,
            ExprBjkl ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBjkl &b, const ExprAjkli &a){

	typedef Ajkli_to_Aijkl < ExprAjkli ,T> Permuted_Obj1;


	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , ExprBjkl, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),b));
}



/*!brief A(j,k,l,i)*B(j,l,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkli_to_Aijkl < ExprAjkli ,T> ,
            Aikj_to_Aijk< ExprBjlk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjkli &a,const ExprBjlk &b){

	typedef Ajkli_to_Aijkl < ExprAjkli ,T> Permuted_Obj1;
	typedef Aikj_to_Aijk< ExprBjlk,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief B(j,l,k)*A(j,k,l,i) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkli_to_Aijkl < ExprAjkli ,T> ,
            Aikj_to_Aijk< ExprBjlk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBjlk &b, const ExprAjkli &a){

	typedef Ajkli_to_Aijkl < ExprAjkli ,T> Permuted_Obj1;
	typedef Aikj_to_Aijk< ExprBjlk,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(j,k,l,i)*B(k,j,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkli_to_Aijkl < ExprAjkli ,T> ,
            Ajik_to_Aijk< ExprBkjl,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjkli &a,const ExprBkjl &b){

	typedef Ajkli_to_Aijkl < ExprAjkli ,T> Permuted_Obj1;
	typedef Ajik_to_Aijk< ExprBkjl,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief B(k,j,l)*A(j,k,l,i) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkli_to_Aijkl < ExprAjkli ,T> ,
            Ajik_to_Aijk< ExprBkjl,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBkjl &b, const ExprAjkli &a){

	typedef Ajkli_to_Aijkl < ExprAjkli ,T> Permuted_Obj1;
	typedef Ajik_to_Aijk< ExprBkjl,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}




/*!brief A(j,k,l,i)*B(k,l,j) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkli_to_Aijkl < ExprAjkli ,T> ,
            Ajki_to_Aijk< ExprBklj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjkli &a,const ExprBklj &b){

	typedef Ajkli_to_Aijkl < ExprAjkli ,T> Permuted_Obj1;
	typedef Ajki_to_Aijk< ExprBklj,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief B(k,l,j)*A(j,k,l,i) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkli_to_Aijkl < ExprAjkli ,T> ,
            Ajki_to_Aijk< ExprBklj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBklj &b,const ExprAjkli &a){

	typedef Ajkli_to_Aijkl < ExprAjkli ,T> Permuted_Obj1;
	typedef Ajki_to_Aijk< ExprBklj,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(j,k,l,i)*B(l,j,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkli_to_Aijkl < ExprAjkli ,T> ,
            Akij_to_Aijk< ExprBljk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjkli &a,const ExprBljk &b){

	typedef Ajkli_to_Aijkl < ExprAjkli ,T> Permuted_Obj1;
	typedef Akij_to_Aijk< ExprBljk,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief B(l,j,k)*A(j,k,l,i) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkli_to_Aijkl < ExprAjkli ,T> ,
            Akij_to_Aijk< ExprBljk,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBljk &b, const ExprAjkli &a){

	typedef Ajkli_to_Aijkl < ExprAjkli ,T> Permuted_Obj1;
	typedef Akij_to_Aijk< ExprBljk,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief A(j,k,l,i)*B(l,k,j) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkli_to_Aijkl < ExprAjkli ,T> ,
            Akji_to_Aijk< ExprBlkj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprAjkli &a,const ExprBlkj &b){

	typedef Ajkli_to_Aijkl < ExprAjkli ,T> Permuted_Obj1;
	typedef Akji_to_Aijk< ExprBlkj,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief B(l,k,j)*A(j,k,l,i) */

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr1<const Aijkl_contracts_Bjkl<
            Ajkli_to_Aijkl < ExprAjkli ,T> ,
            Akji_to_Aijk< ExprBlkj,U> ,
            typename promotedType>,
            typename promotedType,i>
operator*(const ExprBlkj &b, const ExprAjkli &a){

	typedef Ajkli_to_Aijkl < ExprAjkli ,T> Permuted_Obj1;
	typedef Akji_to_Aijk< ExprBlkj,U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bjkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr1 < ExprObj, typename promotedType, i > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}




#undef ExprAijkl
#undef ExprAijlk
#undef ExprAklij
#undef ExprAlkij
#undef ExprAikjl
#undef ExprAjikl
#undef ExprAiklj
#undef ExprAkijl
#undef ExprAkilj
#undef ExprAjkil
#undef ExprAjkli

#undef ExprBklm
#undef ExprBmkl
#undef ExprBmlk
#undef ExprBkml
#undef ExprBlkm
#undef ExprBlmk
#undef ExprBjkl
#undef ExprBjlk
#undef ExprBkjl
#undef ExprBklj
#undef ExprBljk
#undef ExprBlkj


#undef promotedType

#endif
