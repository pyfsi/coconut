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
#ifndef Expr4_contract_Expr2_H
#define Expr4_contract_Expr2_H


////////////////////////////////////////////
/*!brief A(i,j,k,l) * B(l,m) -> C(i,j,k,m) */
///////////////////////////////////////////

#define ExprAijkl	Expr4<A,T,i,j,k,l>
#define ExprAijlk	Expr4<A,T,i,j,l,k>
#define ExprAiljk	Expr4<A,T,i,l,j,k>
#define ExprAlijk	Expr4<A,T,l,i,j,k>
#define ExprAklij	Expr4<A,T,k,l,i,j>
#define ExprAikjl	Expr4<A,T,i,k,j,l>
#define ExprAiklj	Expr4<A,T,i,k,l,j>
#define ExprAkijl	Expr4<A,T,k,i,j,l>
#define ExprAkilj	Expr4<A,T,k,i,l,j>

#define ExprBlm	Expr2<B,U,l,m>
#define ExprBml	Expr2<B,U,m,l>
#define ExprBkl	Expr2<B,U,k,l>
#define ExprBlk	Expr2<B,U,l,k>
#define promotedType	promote<T,U>::V

/////////////////////////////////////////////////////////////////////
// Single Contraction with first index of Second Order Tensor
////////////////////////////////////////////////////////////////////

/*!brief A(i,j,k,l)*B(l,m) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm < ExprAijkl , ExprBlm , typename promotedType>,
				typename promotedType, i , j, k, m >
operator* (const ExprAijkl &a, const ExprBlm &b)
{
	typedef  const Aijkl_contracts_Blm < ExprAijkl , ExprBlm , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (a,b));
}

/*!brief B(l,m)*A(i,j,k,l) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm < ExprAijkl , ExprBlm , typename promotedType>,
				typename promotedType, i , j, k, m >
operator* (const ExprBlm &b, const ExprAijkl &a)
{
	typedef  const Aijkl_contracts_Blm < ExprAijkl , ExprBlm , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (a,b));
}



/*!brief A(i,j,l,k)*B(l,m) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm < Aijlk_to_Aijkl < ExprAijlk ,T > ,
			 ExprBlm , typename promotedType>, typename promotedType, i, j, k, m >
operator* (const ExprAijlk &a, const ExprBlm &b)
{
	typedef Aijlk_to_Aijkl < ExprAijlk ,T > Permuted_Obj;
	typedef const Aijkl_contracts_Blm < Permuted_Obj , ExprBlm , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (Permuted_Obj(a),b));
}

/*!brief B(l,m)*A(i,j,l,k) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm < Aijlk_to_Aijkl < ExprAijlk ,T> ,
			 ExprBlm , typename promotedType>, typename promotedType, i, j, k, m >
operator* (const ExprBlm &b, const ExprAijlk &a)
{
	typedef Aijlk_to_Aijkl < ExprAijlk ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Blm < Permuted_Obj , ExprBlm , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (Permuted_Obj(a),b));
}


/*!brief A(i,l,j,k)*B(l,m) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm < Ailjk_to_Aijkl < ExprAiljk ,T> ,
			 ExprBlm , typename promotedType>, typename promotedType, i, j, k, m >
operator* (const ExprAiljk &a, const ExprBlm &b)
{
	typedef Ailjk_to_Aijkl < ExprAiljk ,T > Permuted_Obj;
	typedef Ailjk_to_Aijkl < ExprAiljk ,T > Permuted_Obj;
	typedef const Aijkl_contracts_Blm < Permuted_Obj , ExprBlm , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (Permuted_Obj(a),b));
}

/*!brief B(l,m)*A(i,l,j,k) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm < Ailjk_to_Aijkl < ExprAiljk ,T> ,
			 ExprBlm , typename promotedType>, typename promotedType, i, j, k, m >
operator* (const ExprBlm &b, const ExprAiljk &a)
{
	typedef Ailjk_to_Aijkl < ExprAiljk ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Blm < Permuted_Obj , ExprBlm , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (Permuted_Obj(a),b));
}

/*!brief A(l,i,j,k)*B(l,m) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm < Alijk_to_Aijkl < ExprAlijk ,T > ,
			 ExprBlm , typename promotedType>, typename promotedType, i, j, k, m >
operator* (const ExprAlijk &a, const ExprBlm &b)
{
	typedef Alijk_to_Aijkl < ExprAlijk ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Blm < Permuted_Obj , ExprBlm , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (Permuted_Obj(a),b));
}

/*!brief B(l,m)*A(l,i,j,k) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm < Alijk_to_Aijkl < ExprAlijk ,T> ,
			 ExprBlm , typename promotedType>, typename promotedType, i, j, k, m >
operator* (const ExprBlm &b, const ExprAlijk &a)
{
	typedef Alijk_to_Aijkl < ExprAlijk ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Blm < Permuted_Obj , ExprBlm , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (Permuted_Obj(a),b));
}

/////////////////////////////////////////////////////////////////////
// Single Contraction with Second index of Second Order Tensor
////////////////////////////////////////////////////////////////////


/*!brief A(i,j,k,l)*B(m,l) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm < ExprAijkl , Aji_to_Aij < ExprBml, U> , typename promotedType>,
				typename promotedType, i , j, k, m >
operator* (const ExprAijkl &a, const ExprBml &b)
{
	typedef Aji_to_Aij < ExprBml ,U> Permuted_Obj1;
	typedef  const Aijkl_contracts_Blm < ExprAijkl , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (a,b));
}

/*!brief B(m,l)*A(i,j,k,l) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm < ExprAijkl , Aji_to_Aij < ExprBml, U > , typename promotedType>,
				typename promotedType, i , j, k, m >
operator* (const ExprBml &b, const ExprAijkl &a)
{
	typedef Aji_to_Aij < ExprBml ,U > Permuted_Obj1;
	typedef  const Aijkl_contracts_Blm < ExprAijkl , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (a,Permuted_Obj1(b)));
}


/*!brief A(i,j,l,k)*B(m,l) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm <
				Aijlk_to_Aijkl < ExprAijlk ,T >,
			 	Aji_to_Aij < ExprBml, U>,
				typename promotedType>,
		 typename promotedType, i, j, k, m >
operator* (const ExprAijlk &a, const ExprBml &b)
{
	typedef Aijlk_to_Aijkl < ExprAijlk ,T > Permuted_Obj;
	typedef Aji_to_Aij < ExprBml ,U > Permuted_Obj1;
	typedef const Aijkl_contracts_Blm < Permuted_Obj , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (Permuted_Obj(a),Permuted_Obj1(b)));
}

/*!brief B(m,l)*A(i,j,l,k) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm <
				Aijlk_to_Aijkl < ExprAijlk ,T >,
			 	Aji_to_Aij < ExprBml, U>,
				typename promotedType>,
		 typename promotedType, i, j, k, m >
operator* (const ExprBml &b, const ExprAijlk &a)
{
	typedef Aijlk_to_Aijkl < ExprAijlk ,T > Permuted_Obj;
	typedef Aji_to_Aij < ExprBml ,U> Permuted_Obj1;
	typedef const Aijkl_contracts_Blm < Permuted_Obj , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (Permuted_Obj(a),Permuted_Obj1(b)));
}

/*!brief A(i,l,j,k)*B(m,l) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm <
				Ailjk_to_Aijkl < ExprAiljk ,T>,
			 	Aji_to_Aij < ExprBml, U>,
				typename promotedType>,
		 typename promotedType, i, j, k, m >
operator* (const ExprAiljk &a, const ExprBml &b)
{
	typedef Ailjk_to_Aijkl < ExprAiljk ,T> Permuted_Obj;
	typedef Aji_to_Aij < ExprBml ,U > Permuted_Obj1;
	typedef const Aijkl_contracts_Blm < Permuted_Obj , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (Permuted_Obj(a),Permuted_Obj1(b)));
}

/*!brief B(m,l)*A(i,l,j,k) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm <
				Ailjk_to_Aijkl < ExprAiljk ,T >,
			 	Aji_to_Aij < ExprBml, U>,
				typename promotedType>,
		 typename promotedType, i, j, k, m >
operator* (const ExprBml &b, const ExprAiljk &a)
{
	typedef Ailjk_to_Aijkl < ExprAiljk ,T> Permuted_Obj;
	typedef Aji_to_Aij < ExprBml ,U > Permuted_Obj1;
	typedef const Aijkl_contracts_Blm < Permuted_Obj , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (Permuted_Obj(a),Permuted_Obj1(b)));
}

/*!brief A(l,i,j,k)*B(m,l) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm <
				Alijk_to_Aijkl < ExprAlijk ,T>,
			 	Aji_to_Aij < ExprBml, U>,
				typename promotedType>,
		 typename promotedType, i, j, k, m >
operator* (const ExprAlijk &a, const ExprBml &b)
{
	typedef Alijk_to_Aijkl < ExprAlijk ,T> Permuted_Obj;
	typedef Aji_to_Aij < ExprBml ,U> Permuted_Obj1;
	typedef const Aijkl_contracts_Blm < Permuted_Obj , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (Permuted_Obj(a),Permuted_Obj1(b)));
}

/*!brief B(m,l)*A(l,i,j,k) ==> C(i,j,k,m) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijkl_contracts_Blm <
				Alijk_to_Aijkl < ExprAlijk ,T>,
			 	Aji_to_Aij < ExprBml, U>,
				typename promotedType>,
		 typename promotedType, i, j, k, m >
operator* (const ExprBml &b, const ExprAlijk &a)
{
	typedef Alijk_to_Aijkl < ExprAlijk ,T> Permuted_Obj;
	typedef Aji_to_Aij < ExprBml ,U > Permuted_Obj1;
	typedef const Aijkl_contracts_Blm < Permuted_Obj , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType , i, j, k, m > (ExprObj (Permuted_Obj(a),Permuted_Obj1(b)));
}

/////////////////////////////////////////////////////////////////////
// Double Index Contraction with Second Order Tensor
////////////////////////////////////////////////////////////////////

/*!brief A(i,j,k,l)*B(k,l) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl < ExprAijkl , ExprBkl , typename promotedType>,
				typename promotedType, i , j>
operator* (const ExprAijkl &a, const ExprBkl &b)
{
	typedef  const Aijkl_contracts_Bkl < ExprAijkl , ExprBkl , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j> (ExprObj (a,b));
}

/*!brief B(k,l)*A(i,j,k,l) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl < ExprAijkl , ExprBkl , typename promotedType>,
				typename promotedType, i , j>
operator* (const ExprBkl &b, const ExprAijkl &a)
{
	typedef  const Aijkl_contracts_Bkl < ExprAijkl , ExprBkl , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j> (ExprObj (a,b));
}

/*!brief A(i,j,k,l)*B(l,k) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl < ExprAijkl,
			 	Aji_to_Aij < ExprBlk, U>,
				typename promotedType>, typename promotedType, i, j>
operator* (const ExprAijkl &a, const ExprBlk &b)
{
	typedef Aji_to_Aij < ExprBlk ,U> Permuted_Obj1;
	typedef const Aijkl_contracts_Bkl < ExprAijkl , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j> (ExprObj (a,Permuted_Obj1(b)));
}

/*!brief B(l,k)*A(i,j,k,l) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl < ExprAijkl,
			 	Aji_to_Aij < ExprBlk, U>,
				typename promotedType>, typename promotedType, i, j>
operator* (const ExprBlk &b, const ExprAijkl &a)
{
	typedef Aji_to_Aij < ExprBlk ,U > Permuted_Obj1;
	typedef const Aijkl_contracts_Bkl < ExprAijkl , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j> (ExprObj (a,Permuted_Obj1(b)));
}

/*!brief A(k,l,i,j)*B(k,l) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Aklij_to_Aijkl < ExprAklij ,T>,
			 	ExprBkl,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprAklij &a, const ExprBkl &b)
{
	typedef Aklij_to_Aijkl < ExprAklij ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Bkl < Permuted_Obj ,  ExprBkl , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj(a),b));
}

/*!brief B(k,l)*A(k,l,i,j) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Aklij_to_Aijkl < ExprAklij ,T>,
			 	ExprBkl,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprBkl &b, const ExprAklij &a)
{
	typedef Aklij_to_Aijkl < ExprAklij ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Bkl < Permuted_Obj ,  ExprBkl , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj(a),b));
}

/*!brief A(k,l,i,j)*B(l,k) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Aklij_to_Aijkl < ExprAklij ,T>,
			 	Aji_to_Aij < ExprBlk, U>,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprAklij &a, const ExprBlk &b)
{
	typedef Aklij_to_Aijkl < ExprAklij ,T> Permuted_Obj;
	typedef Aji_to_Aij < ExprBlk ,U> Permuted_Obj1;
	typedef const Aijkl_contracts_Bkl < Permuted_Obj , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj(a),Permuted_Obj1(b)));
}

/*!brief B(l,k)*A(k,l,i,j) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Aklij_to_Aijkl < ExprAklij ,T>,
			 	Aji_to_Aij < ExprBlk, U>,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprBlk &b, const ExprAklij &a)
{
	typedef Aklij_to_Aijkl < ExprAklij ,T> Permuted_Obj;
	typedef Aji_to_Aij < ExprBlk ,U> Permuted_Obj1;
	typedef const Aijkl_contracts_Bkl < Permuted_Obj , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj(a),Permuted_Obj1(b)));
}


/*!brief A(i,k,j,l)*B(k,l) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Aikjl_to_Aijkl < ExprAikjl ,T>,
			 	ExprBkl ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprAikjl &a, const ExprBkl &b)
{
	typedef Aikjl_to_Aijkl < ExprAikjl ,T> Permuted_Obj;

	typedef const Aijkl_contracts_Bkl < Permuted_Obj , ExprBkl , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj(a),b));
}


/*!brief B(k,l)*A(i,k,j,l) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Aikjl_to_Aijkl < ExprAikjl ,T>,
			 	ExprBkl ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprBkl &b, const ExprAikjl &a)
{
	typedef Aikjl_to_Aijkl < ExprAikjl ,T> Permuted_Obj;

	typedef const Aijkl_contracts_Bkl < Permuted_Obj , ExprBkl , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj(a),b));
}



/*!brief A(i,k,j,l)*B(l,k) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Aikjl_to_Aijkl < ExprAikjl ,T>,
			 	Aji_to_Aij<ExprBlk, U>  ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprAikjl &a, const ExprBlk &b)
{
	typedef Aikjl_to_Aijkl < ExprAikjl ,T> Permuted_Obj1;
	typedef Aji_to_Aij<ExprBlk, U>  Permuted_Obj2;

	typedef const Aijkl_contracts_Bkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief B(l,k)*A(i,k,j,l) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Aikjl_to_Aijkl < ExprAikjl ,T>,
			 	Aji_to_Aij<ExprBlk, U>  ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprBlk &b, const ExprAikjl &a)
{
	typedef Aikjl_to_Aijkl < ExprAikjl ,T> Permuted_Obj1;
	typedef Aji_to_Aij<ExprBlk, U>  Permuted_Obj2;

	typedef const Aijkl_contracts_Bkl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief A(i,k,l,j)*B(k,l) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Aiklj_to_Aijkl < ExprAiklj ,T>,
			 	ExprBkl ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprAiklj &a, const ExprBkl &b)
{
	typedef Aiklj_to_Aijkl < ExprAiklj ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Bkl < Permuted_Obj , ExprBkl , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj(a),b));
}


/*!brief B(k,l)*A(i,k,l,j)==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Aiklj_to_Aijkl < ExprAiklj ,T>,
			 	ExprBkl ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprBkl &b, const ExprAiklj &a)
{
	typedef Aiklj_to_Aijkl < ExprAiklj ,T> Permuted_Obj;

	typedef const Aijkl_contracts_Bkl < Permuted_Obj , ExprBkl , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj(a),b));
}



/*!brief A(i,k,l,j)*B(l,k) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Aiklj_to_Aijkl < ExprAiklj ,T>,
			 	Aji_to_Aij<ExprBlk , U> ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprAiklj &a, const ExprBlk &b)
{
	typedef Aiklj_to_Aijkl < ExprAiklj ,T> Permuted_Obj1;
	typedef Aji_to_Aij<ExprBlk , U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bkl < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!brief B(l,k)*A(i,k,l,j)==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Aiklj_to_Aijkl < ExprAiklj ,T>,
			 	Aji_to_Aij<ExprBlk , U> ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprBlk &b, const ExprAiklj &a)
{
	typedef Aiklj_to_Aijkl < ExprAiklj ,T> Permuted_Obj1;
	typedef Aji_to_Aij<ExprBlk , U> Permuted_Obj2;

	typedef const Aijkl_contracts_Bkl < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief A(k,i,j,l)*B(k,l) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Akijl_to_Aijkl < ExprAkijl ,T>,
			 	ExprBkl  ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprAkijl &a, const ExprBkl &b)
{
	typedef Akijl_to_Aijkl < ExprAkijl ,T> Permuted_Obj1;

	typedef const Aijkl_contracts_Bkl < Permuted_Obj1 , ExprBkl , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj1(a),b));
}


/*!brief B(k,l)*A(k,i,j,l)==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Akijl_to_Aijkl < ExprAkijl ,T>,
			 	ExprBkl  ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprBkl &b, const ExprAkijl &a)
{
	typedef Akijl_to_Aijkl < ExprAkijl ,T> Permuted_Obj1;

	typedef const Aijkl_contracts_Bkl < Permuted_Obj1 , ExprBkl , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj1(a),b));
}


/*!brief A(k,i,j,l)*B(l,k) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Akijl_to_Aijkl < ExprAkijl ,T>,
			 	Aji_to_Aij<ExprBlk ,U>  ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprAkijl &a, const ExprBlk &b)
{
	typedef Akijl_to_Aijkl < ExprAkijl ,T> Permuted_Obj1;
	typedef Aji_to_Aij<ExprBlk ,U> Permuted_Obj2;
	typedef const Aijkl_contracts_Bkl < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief B(l,k)*A(k,i,j,l) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Akijl_to_Aijkl < ExprAkijl ,T>,
			 	Aji_to_Aij<ExprBlk ,U>  ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* ( const ExprBlk &b,const ExprAkijl &a)
{
	typedef Akijl_to_Aijkl < ExprAkijl ,T> Permuted_Obj1;
	typedef Aji_to_Aij<ExprBlk ,U> Permuted_Obj2;
	typedef const Aijkl_contracts_Bkl < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief A(k,i,l,j)*B(l,k) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Akilj_to_Aijkl < ExprAkilj ,T>,
			 	Aji_to_Aij<ExprBlk ,U>  ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprAkilj &a, const ExprBlk &b)
{
	typedef Akilj_to_Aijkl < ExprAkilj ,T> Permuted_Obj1;
	typedef Aji_to_Aij<ExprBlk ,U> Permuted_Obj2;
	typedef const Aijkl_contracts_Bkl < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!brief B(l,k)*A(k,i,l,j) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Akilj_to_Aijkl < ExprAkilj ,T>,
			 	Aji_to_Aij<ExprBlk ,U>  ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprBlk &b, const ExprAkilj &a)
{
	typedef Akilj_to_Aijkl < ExprAkilj ,T> Permuted_Obj1;
	typedef Aji_to_Aij<ExprBlk ,U> Permuted_Obj2;
	typedef const Aijkl_contracts_Bkl < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}




/*!brief A(k,i,l,j)*B(k,l) ==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Akilj_to_Aijkl < ExprAkilj ,T>,
			 	ExprBkl  ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprAkilj &a, const ExprBkl &b)
{
	typedef Akilj_to_Aijkl < ExprAkilj ,T> Permuted_Obj1;

	typedef const Aijkl_contracts_Bkl < Permuted_Obj1 , ExprBkl , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj1(a),b));
}



/*!brief B(k,l)*A(k,i,l,j)==> C(i,j) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijkl_contracts_Bkl <
				Akilj_to_Aijkl < ExprAkilj ,T>,
			 	ExprBkl  ,
				typename promotedType>,
		 typename promotedType, i, j>
operator* (const ExprBkl &b, const ExprAkilj &a)
{
	typedef Akilj_to_Aijkl < ExprAkilj ,T> Permuted_Obj1;

	typedef const Aijkl_contracts_Bkl < Permuted_Obj1 , ExprBkl , typename promotedType> ExprObj;
	return Expr2 < ExprObj, typename promotedType , i, j > (ExprObj (Permuted_Obj1(a),b));
}


#undef ExprAijkl
#undef ExprAijlk
#undef ExprAiljk
#undef ExprAlijk
#undef ExprAklij
#undef ExprAikjl
#undef ExprAiklj
#undef ExprAkijl
#undef ExprAkilj

#undef ExprBlm
#undef ExprBml
#undef ExprBkl
#undef ExprBlk
#undef promotedType

#endif
