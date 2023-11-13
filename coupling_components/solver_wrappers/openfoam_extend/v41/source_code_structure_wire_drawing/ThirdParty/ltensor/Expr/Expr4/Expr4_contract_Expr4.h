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
#ifndef Expr4_contract_Expr4_H
#define Expr4_contract_Expr4_H


////////////////////////////////////////////
/*!brief A(i,j,k,l) * B(k,l,m,n) -> C(i,j,m,n) */
///////////////////////////////////////////

#define ExprAijkl	Expr4<A,T,i,j,k,l>
#define ExprAijlk	Expr4<A,T,i,j,l,k>
#define ExprAklij	Expr4<A,T,k,l,i,j>
#define ExprAlkij	Expr4<A,T,l,k,i,j>
#define ExprAikjl	Expr4<A,T,i,k,j,l>
#define ExprAiklj	Expr4<A,T,i,k,l,j>
#define ExprAkijl	Expr4<A,T,k,i,j,l>
#define ExprAkilj	Expr4<A,T,k,i,l,j>
#define ExprAjikl	Expr4<A,T,j,i,k,l>
#define ExprAjkil   Expr4<A,T,j,k,i,l>
#define ExprAjkli   Expr4<A,T,j,k,l,i>
#define ExprAiljk   Expr4<A,T,i,l,j,k>
#define ExprAilkj   Expr4<A,T,i,l,k,j>
#define ExprAjilk   Expr4<A,T,j,i,l,k>
#define ExprAjlik   Expr4<A,T,j,l,i,k>
#define ExprAjlki   Expr4<A,T,j,l,k,i>
#define ExprAkjil   Expr4<A,T,k,j,i,l>


#define ExprBklmn	Expr4<B,U,k,l,m,n>
#define ExprBmnkl	Expr4<B,U,m,n,k,l>
#define ExprBmnlk	Expr4<B,U,m,n,l,k>
#define ExprBkmln	Expr4<B,U,k,m,l,n>
#define ExprBkmnl	Expr4<B,U,k,m,n,l>
#define ExprBlmkn	Expr4<B,U,l,m,k,n>
#define ExprBlmnk	Expr4<B,U,l,m,n,k>
#define ExprBmkln	Expr4<B,U,m,k,l,n>
#define ExprBmknl	Expr4<B,U,m,k,n,l>
#define ExprBmlkn	Expr4<B,U,m,l,k,n>
#define ExprBmlnk	Expr4<B,U,m,l,n,k>
#define ExprBlkmn	Expr4<B,U,l,k,m,n>
#define ExprBjklm	Expr4<B,U,j,k,l,m>
#define ExprBjkml	Expr4<B,U,j,k,m,l>
#define ExprBjlkm	Expr4<B,U,j,l,k,m>
#define ExprBjlmk	Expr4<B,U,j,l,m,k>
#define ExprBjmkl	Expr4<B,U,j,m,k,l>
#define ExprBjmlk	Expr4<B,U,j,m,l,k>
#define ExprBkjlm	Expr4<B,U,k,j,l,m>
#define ExprBkjml   Expr4<B,U,k,j,m,l>
#define ExprBkljm   Expr4<B,U,k,l,j,m>
#define ExprBklmj   Expr4<B,U,k,l,m,j>
#define ExprBkmjl   Expr4<B,U,k,m,j,l>
#define ExprBkmlj   Expr4<B,U,k,m,l,j>
#define ExprBljkm   Expr4<B,U,l,j,k,m>
#define ExprBljmk   Expr4<B,U,l,j,m,k>
#define ExprBlkjm   Expr4<B,U,l,k,j,m>
#define ExprBlkmj   Expr4<B,U,l,k,m,j>
#define ExprBlmjk   Expr4<B,U,l,m,j,k>
#define ExprBlmkj   Expr4<B,U,l,m,k,j>
#define ExprBmjkl   Expr4<B,U,m,j,k,l>
#define ExprBmjlk   Expr4<B,U,m,j,l,k>
#define ExprBmkjl   Expr4<B,U,m,k,j,l>
#define ExprBmklj   Expr4<B,U,m,k,l,j>
#define ExprBmljk   Expr4<B,U,m,l,j,k>
#define ExprBmlkj   Expr4<B,U,m,l,k,j>
#define ExprBijkl   Expr4<B,U,i,j,k,l>
#define ExprBijlk   Expr4<B,U,i,j,l,k>
#define ExprBikjl   Expr4<B,U,i,k,j,l>
#define ExprBiklj   Expr4<B,U,i,k,l,j>
#define ExprBiljk   Expr4<B,U,i,l,j,k>
#define ExprBilkj   Expr4<B,U,i,l,k,j>
#define ExprBjikl   Expr4<B,U,j,i,k,l>
#define ExprBjilk   Expr4<B,U,j,i,l,k>
#define ExprBjkil   Expr4<B,U,j,k,i,l>
#define ExprBjkli   Expr4<B,U,j,k,l,i>
#define ExprBjlik   Expr4<B,U,j,l,i,k>
#define ExprBjlki   Expr4<B,U,j,l,k,i>
#define ExprBkijl   Expr4<B,U,k,i,j,l>
#define ExprBkilj   Expr4<B,U,k,i,l,j>
#define ExprBkjil   Expr4<B,U,k,j,i,l>
#define ExprBkjli   Expr4<B,U,k,j,l,i>
#define ExprBklij   Expr4<B,U,k,l,i,j>
#define ExprBklji   Expr4<B,U,k,l,j,i>
#define ExprBlijk   Expr4<B,U,l,i,j,k>
#define ExprBlikj   Expr4<B,U,l,i,k,j>
#define ExprBljik   Expr4<B,U,l,j,i,k>
#define ExprBljki   Expr4<B,U,l,j,k,i>
#define ExprBlkij   Expr4<B,U,l,k,i,j>
#define ExprBlkji   Expr4<B,U,l,k,j,i>

#define promotedType	promote<T,U>::V

/////////////////////////////////////////////////////////////////////
// Double Contraction with first 2 indices of Fourth Order Tensor
////////////////////////////////////////////////////////////////////

/*!brief A(i,j,k,l)*B(k,l,m,n) ==> C(i,j,m,n) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4 < const Aijkl_contracts_Bklmn < ExprAijkl , ExprBklmn , typename promotedType>,
				typename promotedType, i , j, m, n>
operator* (const ExprAijkl &a, const ExprBklmn &b)
{
	typedef  const Aijkl_contracts_Bklmn < ExprAijkl , ExprBklmn , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType, i, j, m, n > (ExprObj (a,b));
}//CHECKED 1


/*!brief A(i,j,l,k)*B(k,l,m,n) ==> C(i,j,m,n) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4 < const Aijkl_contracts_Bklmn < Aijlk_to_Aijkl < ExprAijlk ,T> ,
			 ExprBklmn , typename promotedType>, typename promotedType, i, j, m, n >
operator* (const ExprAijlk &a, const ExprBklmn &b)
{
	typedef Aijlk_to_Aijkl < ExprAijlk ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Bklmn < Permuted_Obj , ExprBklmn , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType, i, j, m, n > (ExprObj (Permuted_Obj(a),b));
}//CHECKED 2



/////////////////////////////////////////////////////////////////////
// Double Contraction of first 2 indices with first 2 indices of Fourth Order Tensors
////////////////////////////////////////////////////////////////////

/*!brief A(k,l,i,j)*B(k,l,m,n) ==> C(i,j,m,n) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4 < const Aijkl_contracts_Bklmn <
				Aklij_to_Aijkl < ExprAklij ,T>,
			 	ExprBklmn,
				typename promotedType>,
		 typename promotedType, i, j, m, n>
operator* (const ExprAklij &a, const ExprBklmn &b)
{
	typedef Aklij_to_Aijkl < ExprAklij ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Bklmn < Permuted_Obj ,  ExprBklmn , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType, i, j, m, n > (ExprObj (Permuted_Obj(a),b));
}//CHECKED 3

/*!brief A(l,k,i,j)*B(k,l,m,n) ==> C(i,j,m,n) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4 < const Aijkl_contracts_Bklmn <
				Alkij_to_Aijkl < ExprAlkij ,T>,
			 	ExprBklmn,
				typename promotedType>,
		 typename promotedType, i, j, m, n>
operator* (const ExprAlkij &a, const ExprBklmn &b)
{
	typedef Alkij_to_Aijkl < ExprAlkij ,T> Permuted_Obj;
	typedef const Aijkl_contracts_Bklmn < Permuted_Obj ,  ExprBklmn , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType, i, j, m, n > (ExprObj (Permuted_Obj(a),b));
}//CHECKED 4


/////////////////////////////////////////////////////////////////////
// Double Contraction of last 2 indices with last 2 indices of Fourth Order Tensor
////////////////////////////////////////////////////////////////////

/*!brief A(i,j,k,l)*B(m,n,k,l) ==> C(i,j,m,n) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4 < const Aijkl_contracts_Bklmn < ExprAijkl,
			 	Aklij_to_Aijkl < ExprBmnkl, U>,
				typename promotedType>, typename promotedType, i, j, m, n>
operator* (const ExprAijkl &a, const ExprBmnkl &b)
{
	typedef Aklij_to_Aijkl < ExprBmnkl ,U> Permuted_Obj1;
	typedef const Aijkl_contracts_Bklmn < ExprAijkl , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType, i, j, m, n> (ExprObj (a,Permuted_Obj1(b)));
}//CHECKED 5

/*!brief A(i,j,k,l)*B(m,n,l,k) ==> C(i,j,m,n) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4 < const Aijkl_contracts_Bklmn < ExprAijkl,
			 	Aklji_to_Aijkl < ExprBmnlk, U>,
				typename promotedType>, typename promotedType, i, j, m, n>
operator* (const ExprAijkl &a, const ExprBmnlk &b)
{
	typedef Aklji_to_Aijkl < ExprBmnlk ,U> Permuted_Obj1;
	typedef const Aijkl_contracts_Bklmn < ExprAijkl , Permuted_Obj1 , typename promotedType> ExprObj;
	return Expr4 < ExprObj, typename promotedType, i, j, m, n> (ExprObj (a,Permuted_Obj1(b)));
}//CHECKED 6


/*!brief A(i,j,k,l)*B(k,m,l,n)=C(i,j,m,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4<const Aijkl_contracts_Bklmn < ExprAijkl, Aikjl_to_Aijkl<ExprBkmln,U>,typename promotedType>
                    ,typename promotedType,i,j,m,n>
operator*(const ExprAijkl &a, const ExprBkmln &b){
    typedef Aikjl_to_Aijkl<ExprBkmln,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklmn< ExprAijkl, Permuted_Obj1, typename promotedType> ExprObj;
    return Expr4 < ExprObj,typename promotedType,i,j,m,n> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 7


/*!brief A(i,j,k,l)*B(k,m,n,l)=C(i,j,m,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4<const Aijkl_contracts_Bklmn<ExprAijkl,Aiklj_to_Aijkl<ExprBkmnl,U>,typename promotedType>
                    , typename promotedType,i,j,m,n>
operator*(const ExprAijkl &a, const ExprBkmnl &b){
    typedef Aiklj_to_Aijkl<ExprBkmnl,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklmn< ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr4< ExprObj,typename promotedType,i,j,m,n> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 8


/*!brief A(i,j,k,l)*B(l,m,k,n)=C(i,j,m,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4<const Aijkl_contracts_Bklmn<ExprAijkl,Ajkil_to_Aijkl<ExprBlmkn,U>,typename promotedType>
                    ,typename promotedType,i,j,m,n>
operator*(const ExprAijkl &a, const ExprBlmkn &b){

    typedef Ajkil_to_Aijkl<ExprBlmkn,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklmn< ExprAijkl, Permuted_Obj1,typename promotedType> ExprObj;
    return Expr4<ExprObj, typename promotedType,i,j,m,n> (ExprObj(a,Permuted_Obj1(b)));

}//CHECKED 9


/*!brief A(i,j,k,l)*B(l,m,n,k)=C(i,j,m,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4 <const Aijkl_contracts_Bklmn<ExprAijkl, Ajkli_to_Aijkl<ExprBlmnk,U>, typename promotedType>
                    , typename promotedType,i,j,m,n>
operator*(const ExprAijkl &a, const ExprBlmnk &b){

    typedef Ajkli_to_Aijkl<ExprBlmnk,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklmn< ExprAijkl, Permuted_Obj1,typename promotedType> ExprObj;
    return Expr4<ExprObj, typename promotedType,i,j,m,n> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 10

/*!brief A(i,j,k,l)*B(m,k,l,n)=C(i,j,m,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4< const Aijkl_contracts_Bklmn<ExprAijkl,Akijl_to_Aijkl<ExprBmkln,U>, typename promotedType>
                    , typename promotedType,i,j,m,n>
operator*(const ExprAijkl &a, const ExprBmkln &b){

    typedef Akijl_to_Aijkl<ExprBmkln,U> Permuted_Obj1;
     typedef const Aijkl_contracts_Bklmn< ExprAijkl, Permuted_Obj1,typename promotedType> ExprObj;
    return Expr4<ExprObj, typename promotedType,i,j,m,n> (ExprObj(a,Permuted_Obj1(b)));


}//CHECKED 11


/*!brief A(i,j,k,l)*B(m,k,n,l)=C(i,j,m,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4<const Aijkl_contracts_Bklmn<ExprAijkl,Akilj_to_Aijkl<ExprBmknl,U>, typename promotedType>
                    , typename promotedType,i,j,m,n>
operator*(const ExprAijkl &a, const ExprBmknl &b){
    typedef Akilj_to_Aijkl<ExprBmknl,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklmn< ExprAijkl, Permuted_Obj1,typename promotedType> ExprObj;
    return Expr4<ExprObj, typename promotedType,i,j,m,n> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 12

/*!brief A(i,j,k,l)*B(m,l,k,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4<const Aijkl_contracts_Bklmn<ExprAijkl,Akjil_to_Aijkl<ExprBmlkn,U>, typename promotedType>
                    , typename promotedType,i,j,m,n>
operator*(const ExprAijkl &a, const ExprBmlkn &b){

    typedef Akjil_to_Aijkl<ExprBmlkn,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklmn< ExprAijkl, Permuted_Obj1,typename promotedType> ExprObj;
    return Expr4<ExprObj, typename promotedType,i,j,m,n> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 13

/*!brief A(i,j,k,l)*B(m,l,n,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4<const Aijkl_contracts_Bklmn<ExprAijkl,Akjli_to_Aijkl<ExprBmlnk,U>, typename promotedType>
                    , typename promotedType,i,j,m,n>
operator*(const ExprAijkl &a, const ExprBmlnk &b){

 typedef Akjli_to_Aijkl<ExprBmlnk,U> Permuted_Obj1;
 typedef const Aijkl_contracts_Bklmn< ExprAijkl, Permuted_Obj1,typename promotedType> ExprObj;
 return Expr4<ExprObj, typename promotedType,i,j,m,n> (ExprObj(a,Permuted_Obj1(b)));

}//CHECKED 14



/*!brief A(i,k,j,l)*B(k,l,m,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4<const Aijkl_contracts_Bklmn<Aikjl_to_Aijkl<ExprAikjl,T>,ExprBklmn, typename promotedType>
                    , typename promotedType,i,j,m,n>
operator*( const ExprAikjl &a, const ExprBklmn &b){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,ExprBklmn, typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),b));

}//CHECKED 15


/*!brief A(i,k,j,l)*B(k,m,l,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4<const Aijkl_contracts_Bklmn<Aikjl_to_Aijkl<ExprAikjl,T>,Aikjl_to_Aijkl<ExprBkmln,U>, typename promotedType>
                    , typename promotedType,i,j,m,n>
operator*(const ExprAikjl &a, const ExprBkmln &b){

    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Aikjl_to_Aijkl<ExprBkmln,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

    } //CHECKED 16

/*!brief A(i,k,j,l)*B(k,m,n,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4<const Aijkl_contracts_Bklmn<Aikjl_to_Aijkl<ExprAikjl,T>,Aiklj_to_Aijkl<ExprBkmnl,U>, typename promotedType>
                    , typename promotedType, i,j,m,n>
operator*(const ExprAikjl &a, const ExprBkmnl &b){

    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Aiklj_to_Aijkl<ExprBkmnl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 17

/*!brief A(i,k,j,l)*B(l,k,m,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline const Expr4<const Aijkl_contracts_Bklmn<Aikjl_to_Aijkl<ExprAikjl,T>,Ajikl_to_Aijkl<ExprBlkmn,U>,typename promotedType>
                    , typename promotedType,i,j,m,n>
operator*(const ExprAikjl &a, const ExprBlkmn &b){

    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Ajikl_to_Aijkl<ExprBlkmn,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 18

/*!brief A(i,k,j,l)*B(l,m,k,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aikjl_to_Aijkl<ExprAikjl,T>,Ajkil_to_Aijkl<ExprBlmkn,U>, typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAikjl &a,const ExprBlmkn &b){
 typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
 typedef Ajkil_to_Aijkl<ExprBlmkn,U> Permuted_Obj2;
 typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
 return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
 } //CHECKED 19

 /*!brief A(i,k,j,l)*B(l,m,n,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aikjl_to_Aijkl<ExprAikjl,T>,Ajkli_to_Aijkl<ExprBlmnk,U>, typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAikjl &a, const ExprBlmnk &b){
     typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
     typedef Ajkli_to_Aijkl<ExprBlmnk,U> Permuted_Obj2;
     typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
     return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 20

/*!brief A(i,k,j,l)*B(m,k,l,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aikjl_to_Aijkl<ExprAikjl,T>,Akijl_to_Aijkl<ExprBmkln,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAikjl &a, const ExprBmkln &b){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Akijl_to_Aijkl<ExprBmkln,U> Permuted_Obj2;
     typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
     return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 21

//A(i,k,j,l)*B(m,k,n,l)

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aikjl_to_Aijkl<ExprAikjl,T>,Akilj_to_Aijkl<ExprBmknl,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAikjl &a, const ExprBmknl &b){

    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Akilj_to_Aijkl<ExprBmknl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 22

/*!brief A(i,k,j,l)*B(m,l,k,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aikjl_to_Aijkl<ExprAikjl,T>,Akjil_to_Aijkl<ExprBmlkn,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAikjl &a, const ExprBmlkn &b){
     typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
     typedef Akjil_to_Aijkl<ExprBmlkn,U> Permuted_Obj2;
     typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
     return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}// CHECKED 23


/*!brief A(i,k,j,l)*B(m,l,n,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aikjl_to_Aijkl<ExprAikjl,T>,Akjli_to_Aijkl<ExprBmlnk,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAikjl &a, const ExprBmlnk &b){
     typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
     typedef Akjli_to_Aijkl<ExprBmlnk,U> Permuted_Obj2;
     typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
     return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 24

/*!brief A(i,k,j,l)*B(m,n,k,l) */


template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aikjl_to_Aijkl<ExprAikjl,T>,Aklij_to_Aijkl<ExprBmnkl,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAikjl &a, const ExprBmnkl &b){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Aklij_to_Aijkl<ExprBmnkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 25


/*!brief A(i,k,j,l)*B(m,n,l,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aikjl_to_Aijkl<ExprAikjl,T>,Aklji_to_Aijkl<ExprBmnlk,U>,typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAikjl &a, const ExprBmnlk &b){
    typedef Aikjl_to_Aijkl<ExprAikjl,T> Permuted_Obj1;
    typedef Aklji_to_Aijkl<ExprBmnlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 26



/*!brief A(i,k,l,j)*B(k,l,m,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aiklj_to_Aijkl<ExprAiklj,T>, ExprBklmn, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAiklj &a, const ExprBklmn &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,ExprBklmn, typename promotedType> ExprObj;
    return Expr4<ExprObj, typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),b));
}//CHECKED 27


/*!brief A(i,k,l,j)*B(k,m,l,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aiklj_to_Aijkl<ExprAiklj,T>,Aikjl_to_Aijkl<ExprBkmln,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAiklj &a, const ExprBkmln &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Aikjl_to_Aijkl<ExprBkmln,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 28


/*!brief A(i,k,l,j)*B(k,m,n,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aiklj_to_Aijkl<ExprAiklj,T>,Aiklj_to_Aijkl<ExprBkmnl,U>,typename promotedType>
                , typename promotedType, i,j,m,n>
operator*( const ExprAiklj &a, const ExprBkmnl &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Aiklj_to_Aijkl<ExprBkmnl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}// CHECKED 29


/*!brief A(i,k,l,j)*B(l,k,m,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aiklj_to_Aijkl<ExprAiklj,T>,Ajikl_to_Aijkl<ExprBlkmn,U>, typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAiklj &a, const ExprBlkmn &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Ajikl_to_Aijkl<ExprBlkmn,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 30


/*!brief A(i,k,l,j)*B(l,m,k,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aiklj_to_Aijkl<ExprAiklj,T>,Ajkil_to_Aijkl<ExprBlmkn,U>,typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAiklj &a, const ExprBlmkn &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Ajkil_to_Aijkl<ExprBlmkn,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 31


/*!brief A(i,k,l,j)*B(l,m,n,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aiklj_to_Aijkl<ExprAiklj,T>,Ajkli_to_Aijkl<ExprBlmnk,U>, typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAiklj &a, const ExprBlmnk &b){
     typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
     typedef Ajkli_to_Aijkl<ExprBlmnk,U> Permuted_Obj2;
     typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
     return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 32


/*!brief A(i,k,l,j)*B(m,k,l,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aiklj_to_Aijkl<ExprAiklj,T>,Akijl_to_Aijkl<ExprBmkln,U>,typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAiklj &a, const ExprBmkln &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Akijl_to_Aijkl<ExprBmkln,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 33

/*!brief A(i,k,l,j)*B(m,k,n,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aiklj_to_Aijkl<ExprAiklj,T>,Akilj_to_Aijkl<ExprBmknl,U>,typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAiklj &a,const ExprBmknl &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Akilj_to_Aijkl<ExprBmknl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
} //CHECKED 34


/*!brief A(i,k,l,j)*B(m,l,k,n) */
template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aiklj_to_Aijkl<ExprAiklj,T>,Akjil_to_Aijkl<ExprBmlkn,U>,typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAiklj &a, const ExprBmlkn &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Akjil_to_Aijkl<ExprBmlkn,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 35

/*!brief A(i,k,l,j)*B(m,l,n,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aiklj_to_Aijkl<ExprAiklj,T>,Akjli_to_Aijkl<ExprBmlnk,U>,typename promotedType>
                ,   typename promotedType,i,j,m,n>
operator*(const ExprAiklj &a, const ExprBmlnk &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Akjli_to_Aijkl<ExprBmlnk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 36

/*!brief A(i,k,l,j)*B(m,n,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aiklj_to_Aijkl<ExprAiklj,T>,Aklij_to_Aijkl<ExprBmnkl,U>,typename promotedType>
                ,   typename promotedType,i,j,m,n>
operator*(const ExprAiklj &a, const ExprBmnkl &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Aklij_to_Aijkl<ExprBmnkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 37


/*!brief A(i,k,l,j)*B(m,n,l,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aiklj_to_Aijkl<ExprAiklj,T>,Aklji_to_Aijkl<ExprBmnlk,U>,typename promotedType>
                ,   typename promotedType,i,j,m,n>
operator*(const ExprAiklj &a, const ExprBmnlk &b){
    typedef Aiklj_to_Aijkl<ExprAiklj,T> Permuted_Obj1;
    typedef Aklji_to_Aijkl<ExprBmnlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 38


/*!brief A(k,i,j,l)*B(k,l,m,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akijl_to_Aijkl<ExprAkijl,T>, ExprBklmn,typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAkijl &a, const ExprBklmn &b){

    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,ExprBklmn,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),b));
} //CHECKED 39

/*!brief A(k,i,j,l)*B(k,m,l,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akijl_to_Aijkl<ExprAkijl,T>,Aikjl_to_Aijkl<ExprBkmln,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAkijl &a, const ExprBkmln &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Aikjl_to_Aijkl<ExprBkmln,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

} //CHECKED 40

/*!brief A(k,i,j,l)*B(k,m,n,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akijl_to_Aijkl<ExprAkijl,T>,Aiklj_to_Aijkl<ExprBkmnl,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*( const ExprAkijl &a, const ExprBkmnl &b){

    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Aiklj_to_Aijkl<ExprBkmnl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 41

/*!brief A(k,i,j,l)*B(l,k,m,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akijl_to_Aijkl<ExprAkijl,T>,Ajikl_to_Aijkl<ExprBlkmn,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAkijl &a,const ExprBlkmn &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Ajikl_to_Aijkl<ExprBlkmn,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
} //CHECKED 42

/*!brief A(k,i,j,l)*B(l,m,k,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akijl_to_Aijkl<ExprAkijl,T>,Ajkil_to_Aijkl<ExprBlmkn,U>, typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAkijl &a, const ExprBlmkn &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Ajkil_to_Aijkl<ExprBlmkn,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 43


/*!brief A(k,i,j,l)*B(l,m,n,k) */


template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akijl_to_Aijkl<ExprAkijl,T>,Ajkli_to_Aijkl<ExprBlmnk,U>, typename promotedType>
                ,   typename promotedType, i,j,m,n>
operator*(const ExprAkijl &a, const ExprBlmnk &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Ajkli_to_Aijkl<ExprBlmnk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 44

/*!brief A(k,i,j,l)*B(m,k,l,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akijl_to_Aijkl<ExprAkijl,T>,Akijl_to_Aijkl<ExprBmkln,U>, typename promotedType>
                ,   typename promotedType, i,j,m,n>
operator*(const ExprAkijl &a, const ExprBmkln &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Akijl_to_Aijkl<ExprBmkln,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 45

/*!brief A(k,i,j,l)*B(m,k,n,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akijl_to_Aijkl<ExprAkijl,T>, Akilj_to_Aijkl<ExprBmknl,U>,typename promotedType>
                ,   typename promotedType, i,j,m,n>
operator*(const ExprAkijl &a, const ExprBmknl &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Akilj_to_Aijkl<ExprBmknl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 46


/*!brief A(k,i,j,l)*B(m,l,k,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akijl_to_Aijkl<ExprAkijl,T>,Akjil_to_Aijkl<ExprBmlkn,U>, typename promotedType>
                , typename promotedType, i ,j ,m ,n>
operator*(const ExprAkijl &a, const ExprBmlkn &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Akjil_to_Aijkl<ExprBmlkn,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
} //CHECKED 47

/*!brief A(k,i,j,l)*B(m,l,n,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akijl_to_Aijkl<ExprAkijl,T>,Akjli_to_Aijkl<ExprBmlnk,U>, typename promotedType>
                ,   typename promotedType,i,j,m,n>
operator*(const ExprAkijl &a, const ExprBmlnk &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Akjli_to_Aijkl<ExprBmlnk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 48

/*!brief A(k,i,j,l)*B(m,n,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akijl_to_Aijkl<ExprAkijl,T>,Aklij_to_Aijkl<ExprBmnkl,U>, typename promotedType>
                ,   typename promotedType,i,j,m,n>
operator*(const ExprAkijl &a, const ExprBmnkl &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Aklij_to_Aijkl<ExprBmnkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
} //CHECKED 49

/*!brief A(k,i,j,l)*B(m,n,l,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akijl_to_Aijkl<ExprAkijl,T>, Aklji_to_Aijkl<ExprBmnlk,U>, typename promotedType>
                , typename promotedType,i,j,m,n>
operator *(const ExprAkijl &a, const ExprBmnlk &b){
    typedef Akijl_to_Aijkl<ExprAkijl,T> Permuted_Obj1;
    typedef Aklji_to_Aijkl<ExprBmnlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 50

/*!brief A(k,i,l,j)*B(k,l,m,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akilj_to_Aijkl<ExprAkilj,T>, ExprBklmn, typename promotedType>
                ,   typename promotedType, i,j,m,n>
operator*(const ExprAkilj &a, const ExprBklmn &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,ExprBklmn,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),b));

}//CHECKED 51

/*!brief A(k,i,l,j)*B(k,m,l,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akilj_to_Aijkl<ExprAkilj,T>, Aikjl_to_Aijkl<ExprBkmln,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*( const ExprAkilj &a, const ExprBkmln &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Aikjl_to_Aijkl<ExprBkmln,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 52


/*!brief A(k,i,l,j)*B(k,m,n,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akilj_to_Aijkl<ExprAkilj,T>, Aiklj_to_Aijkl<ExprBkmnl,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAkilj &a, const ExprBkmnl &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Aiklj_to_Aijkl<ExprBkmnl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
} //CHECKED 53


/*!brief A(k,i,l,j)*B(l,k,m,n) */


template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akilj_to_Aijkl<ExprAkilj,T>,Ajikl_to_Aijkl<ExprBlkmn,U>, typename promotedType>
                ,   typename promotedType,i,j,m,n>
operator*(const ExprAkilj &a, const ExprBlkmn &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Ajikl_to_Aijkl<ExprBlkmn,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 54


/*!brief A(k,i,l,j)*B(l,m,k,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akilj_to_Aijkl<ExprAkilj,T>, Ajkil_to_Aijkl<ExprBlmkn,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAkilj &a, const ExprBlmkn &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Ajkil_to_Aijkl<ExprBlmkn,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 55

/*!brief A(k,i,l,j)*B(l,m,n,k) */


template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akilj_to_Aijkl<ExprAkilj,T>,Ajkli_to_Aijkl<ExprBlmnk,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAkilj &a, const ExprBlmnk &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Ajkli_to_Aijkl<ExprBlmnk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 56


/*!brief A(k,i,l,j)*B(m,k,l,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akilj_to_Aijkl<ExprAkilj,T>, Akijl_to_Aijkl<ExprBmkln,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*( const ExprAkilj &a, const ExprBmkln &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Akijl_to_Aijkl<ExprBmkln,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

} //CHECKED 57


/*!brief A(k,i,l,j)*B(m,k,n,l) */


template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akilj_to_Aijkl<ExprAkilj,T>, Akilj_to_Aijkl<ExprBmknl,U>, typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAkilj &a, const ExprBmknl &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Akilj_to_Aijkl<ExprBmknl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 58

/*!brief A(k,i,l,j)*B(m,l,k,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akilj_to_Aijkl<ExprAkilj,T>, Akjil_to_Aijkl<ExprBmlkn,U>, typename promotedType>
                , typename promotedType,i,j,m,n>
operator*( const ExprAkilj &a, const ExprBmlkn &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Akjil_to_Aijkl<ExprBmlkn,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 59

/*!brief A(k,i,l,j)*B(m,l,n,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akilj_to_Aijkl<ExprAkilj,T>, Akjli_to_Aijkl<ExprBmlnk,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAkilj &a, const ExprBmlnk &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Akjli_to_Aijkl<ExprBmlnk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 60

/*!brief A(k,i,l,j)*B(m,n,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akilj_to_Aijkl<ExprAkilj,T>, Aklij_to_Aijkl<ExprBmnkl,U>, typename promotedType>
                ,   typename promotedType, i,j,m,n>

operator*(const ExprAkilj &a, const ExprBmnkl &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Aklij_to_Aijkl<ExprBmnkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 61


/*!brief A(k,i,l,j)*B(m,n,l,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Akilj_to_Aijkl<ExprAkilj,T>, Aklji_to_Aijkl<ExprBmnlk,U>, typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAkilj &a, const ExprBmnlk &b){
    typedef Akilj_to_Aijkl<ExprAkilj,T> Permuted_Obj1;
    typedef Aklji_to_Aijkl<ExprBmnlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
} //CHECKED 62



/*!brief A(k,l,i,j)*B(k,m,l,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aklij_to_Aijkl<ExprAklij,T>,Aikjl_to_Aijkl<ExprBkmln,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAklij &a, const ExprBkmln &b){
    typedef Aklij_to_Aijkl<ExprAklij,T> Permuted_Obj1;
    typedef Aikjl_to_Aijkl<ExprBkmln,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 63

/*!brief A(k,l,i,j)*B(k,m,n,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aklij_to_Aijkl<ExprAklij,T>,Aiklj_to_Aijkl<ExprBkmnl,U>,typename promotedType>
                ,typename promotedType,i,j,m,n>
operator*(const ExprAklij &a, const ExprBkmnl &b){

    typedef Aklij_to_Aijkl<ExprAklij,T> Permuted_Obj1;
    typedef Aiklj_to_Aijkl<ExprBkmnl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
} // CHECKED 64


/*!brief A(k,l,i,j)*B(l,m,k,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aklij_to_Aijkl<ExprAklij,T>,Ajkil_to_Aijkl<ExprBlmkn,U>, typename promotedType>
                ,typename promotedType, i,j,m,n>
operator*(const ExprAklij &a, const ExprBlmkn &b){
    typedef Aklij_to_Aijkl<ExprAklij,T> Permuted_Obj1;
    typedef Ajkil_to_Aijkl<ExprBlmkn,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 65


/*!brief A(k,l,i,j)*B(l,m,n,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aklij_to_Aijkl<ExprAklij,T>, Ajkli_to_Aijkl<ExprBlmnk,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAklij &a, const ExprBlmnk &b){
    typedef Aklij_to_Aijkl<ExprAklij,T> Permuted_Obj1;
    typedef Ajkli_to_Aijkl<ExprBlmnk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 66


/*!brief A(k,l,i,j)*B(m,k,l,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aklij_to_Aijkl<ExprAklij,T>, Akijl_to_Aijkl<ExprBmkln,U>, typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAklij &a, const ExprBmkln &b){
    typedef Aklij_to_Aijkl<ExprAklij,T> Permuted_Obj1;
    typedef Akijl_to_Aijkl<ExprBmkln,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 67

/*!brief A(k,l,i,j)*B(m,k,n,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aklij_to_Aijkl<ExprAklij,T>, Akilj_to_Aijkl<ExprBmknl,U>, typename promotedType>
                , typename promotedType, i,j,m,n>
operator*(const ExprAklij &a, const ExprBmknl &b){
    typedef Aklij_to_Aijkl<ExprAklij,T> Permuted_Obj1;
    typedef Akilj_to_Aijkl<ExprBmknl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 68

/*!brief A(k,l,i,j)*B(m,l,k,n) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aklij_to_Aijkl<ExprAklij,T>, Akjil_to_Aijkl<ExprBmlkn,U>, typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAklij &a, const ExprBmlkn &b){
    typedef Aklij_to_Aijkl<ExprAklij,T> Permuted_Obj1;
    typedef Akjil_to_Aijkl<ExprBmlkn,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 69


/*!brief A(k,l,i,j)*B(m,l,n,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aklij_to_Aijkl<ExprAklij,T>, Akjli_to_Aijkl<ExprBmlnk,U>, typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAklij &a, const ExprBmlnk &b){
    typedef Aklij_to_Aijkl<ExprAklij,T> Permuted_Obj1;
    typedef Akjli_to_Aijkl<ExprBmlnk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 70


/*!brief A(k,l,i,j)*B(m,n,k,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aklij_to_Aijkl<ExprAklij,T>, Aklij_to_Aijkl<ExprBmnkl,U>, typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAklij &a, const ExprBmnkl &b){
    typedef Aklij_to_Aijkl<ExprAklij,T> Permuted_Obj1;
    typedef Aklij_to_Aijkl<ExprBmnkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 71



/*!brief A(k,l,i,j)*B(m,n,l,k) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m, char n>
inline Expr4<const Aijkl_contracts_Bklmn<Aklij_to_Aijkl<ExprAklij,T>, Aklji_to_Aijkl<ExprBmnlk,U>, typename promotedType>
                , typename promotedType,i,j,m,n>
operator*(const ExprAklij &a, const ExprBmnlk &b){
    typedef Aklij_to_Aijkl<ExprAklij,T> Permuted_Obj1;
    typedef Aklji_to_Aijkl<ExprBmnlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bklmn<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr4<ExprObj,typename promotedType,i,j,m,n> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));

}//CHECKED 72


//FINISHED






//3 Indexes Contractions

//A(i,j,k,l)B*(j,k,l,m) = C(i,m)

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            ExprBjklm,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBjklm &b){
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,ExprBjklm,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,b));
}//CHECKED 74

/*!brief A(i,j,k,l)*B(j,k,m,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Aijlk_to_Aijkl<ExprBjkml,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBjkml &b){
    typedef Aijlk_to_Aijkl<ExprBjkml,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 75

/*!brief A(i,j,k,l)*B(j,l,k,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Aikjl_to_Aijkl<ExprBjlkm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBjlkm &b){
    typedef Aikjl_to_Aijkl<ExprBjlkm,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 76


/*!brief A(i,j,k,l)*B(j,l,m,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Aiklj_to_Aijkl<ExprBjlmk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBjlmk &b){
    typedef Aiklj_to_Aijkl<ExprBjlmk,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 77


/*!brief A(i,j,k,l)*B(j,m,k,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Ailjk_to_Aijkl<ExprBjmkl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBjmkl &b){
    typedef Ailjk_to_Aijkl<ExprBjmkl,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 78


/*!brief A(i,j,k,l)*B(j,m,l,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Ailkj_to_Aijkl<ExprBjmlk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBjmlk &b){
    typedef Ailkj_to_Aijkl<ExprBjmlk,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 79



/*!brief A(i,j,k,l)*B(k,j,l,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Ajikl_to_Aijkl<ExprBkjlm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBkjlm &b){
    typedef Ajikl_to_Aijkl<ExprBkjlm,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 80

/*!brief A(i,j,k,l)*B(k,j,m,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Ajilk_to_Aijkl<ExprBkjml,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBkjml &b){
    typedef Ajilk_to_Aijkl<ExprBkjml,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 81


/*!brief A(i,j,k,l)*B(k,l,j,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Ajkil_to_Aijkl<ExprBkljm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBkljm &b){
    typedef Ajkil_to_Aijkl<ExprBkljm,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 82


/*!brief A(i,j,k,l)*B(k,l,m,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Ajkli_to_Aijkl<ExprBklmj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBklmj &b){
    typedef Ajkli_to_Aijkl<ExprBklmj,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 83



/*!brief A(i,j,k,l)*B(k,m,j,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Ajlik_to_Aijkl<ExprBkmjl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBkmjl &b){
    typedef Ajlik_to_Aijkl<ExprBkmjl,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 84


/*!brief A(i,j,k,l)*B(k,m,l,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Ajlki_to_Aijkl<ExprBkmlj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBkmlj &b){
    typedef Ajlki_to_Aijkl<ExprBkmlj,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 85




/*!brief A(i,j,k,l)*B(l,j,k,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Akijl_to_Aijkl<ExprBljkm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBljkm &b){
    typedef Akijl_to_Aijkl<ExprBljkm,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 86


/*!brief A(i,j,k,l)*B(l,j,m,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Akilj_to_Aijkl<ExprBljmk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBljmk &b){
    typedef Akilj_to_Aijkl<ExprBljmk,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 87



/*!brief A(i,j,k,l)*B(l,k,j,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Akjil_to_Aijkl<ExprBlkjm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBlkjm &b){
    typedef Akjil_to_Aijkl<ExprBlkjm,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 88



/*!brief A(i,j,k,l)*B(l,k,m,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Akjli_to_Aijkl<ExprBlkmj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBlkmj &b){
    typedef Akjli_to_Aijkl<ExprBlkmj,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 89



/*!brief A(i,j,k,l)*B(l,m,j,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Aklij_to_Aijkl<ExprBlmjk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBlmjk &b){
    typedef Aklij_to_Aijkl<ExprBlmjk,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 90



/*!brief A(i,j,k,l)*B(l,m,k,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Aklji_to_Aijkl<ExprBlmkj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBlmkj &b){
    typedef Aklji_to_Aijkl<ExprBlmkj,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 91


/*!brief A(i,j,k,l)*B(m,j,k,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Alijk_to_Aijkl<ExprBmjkl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBmjkl &b){
    typedef Alijk_to_Aijkl<ExprBmjkl,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 92


/*!brief A(i,j,k,l)*B(m,j,l,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Alikj_to_Aijkl<ExprBmjlk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBmjlk &b){
    typedef Alikj_to_Aijkl<ExprBmjlk,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 93


/*!brief A(i,j,k,l)*B(m,k,j,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Aljik_to_Aijkl<ExprBmkjl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBmkjl &b){
    typedef Aljik_to_Aijkl<ExprBmkjl,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 94



/*!brief A(i,j,k,l)*B(m,k,l,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Aljki_to_Aijkl<ExprBmklj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBmklj &b){
    typedef Aljki_to_Aijkl<ExprBmklj,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 95

/*!brief A(i,j,k,l)*B(m,l,j,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Alkij_to_Aijkl<ExprBmljk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBmljk &b){
    typedef Alkij_to_Aijkl<ExprBmljk,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 96


/*!brief A(i,j,k,l)*B(m,l,k,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            ExprAijkl,
            Alkji_to_Aijkl<ExprBmlkj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAijkl &a, const ExprBmlkj &b){
    typedef Alkji_to_Aijkl<ExprBmlkj,U> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<ExprAijkl,Permuted_Obj1,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(a,Permuted_Obj1(b)));
}//CHECKED 97





/*!brief A(j,i,k,l)*B(j,k,l,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            ExprBjklm,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBjklm &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,ExprBjklm,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),b));
}//CHECKED 99



/*!brief A(j,i,k,l)*B(j,k,m,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Aijlk_to_Aijkl<ExprBjkml,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBjkml &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Aijlk_to_Aijkl<ExprBjkml,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 100


/*!brief A(j,i,k,l)*B(j,l,k,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Aikjl_to_Aijkl<ExprBjlkm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBjlkm &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Aikjl_to_Aijkl<ExprBjlkm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 101


/*!brief A(j,i,k,l)*B(j,l,m,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Aiklj_to_Aijkl<ExprBjlmk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBjlmk &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Aiklj_to_Aijkl<ExprBjlmk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 102


/*!brief A(j,i,k,l)*B(j,m,k,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Ailjk_to_Aijkl<ExprBjmkl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBjmkl &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Ailjk_to_Aijkl<ExprBjmkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 103



/*!brief A(j,i,k,l)*B(j,m,l,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Ailkj_to_Aijkl<ExprBjmlk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBjmlk &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Ailkj_to_Aijkl<ExprBjmlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 104



/*!brief A(j,i,k,l)*B(k,j,l,m) */


template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Ajikl_to_Aijkl<ExprBkjlm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBkjlm &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Ajikl_to_Aijkl<ExprBkjlm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 105


/*!brief A(j,i,k,l)*B(k,j,m,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Ajilk_to_Aijkl<ExprBkjml,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBkjml &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Ajilk_to_Aijkl<ExprBkjml,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 106


/*!brief A(j,i,k,l)*B(k,l,j,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Ajkil_to_Aijkl<ExprBkljm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBkljm &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Ajkil_to_Aijkl<ExprBkljm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 107


/*!brief A(j,i,k,l)*B(k,l,m,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Ajkli_to_Aijkl<ExprBklmj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBklmj &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Ajkli_to_Aijkl<ExprBklmj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 108


/*!brief A(j,i,k,l)*B(k,m,j,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Ajlik_to_Aijkl<ExprBkmjl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBkmjl &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Ajlik_to_Aijkl<ExprBkmjl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 109



/*!brief A(j,i,k,l)*B(k,m,l,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Ajlki_to_Aijkl<ExprBkmlj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBkmlj &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Ajlki_to_Aijkl<ExprBkmlj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 110


/*!brief A(j,i,k,l)*B(l,j,k,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Akijl_to_Aijkl<ExprBljkm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBljkm &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Akijl_to_Aijkl<ExprBljkm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 111


/*!brief A(j,i,k,l)*B(l,j,m,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Akilj_to_Aijkl<ExprBljmk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBljmk &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Akilj_to_Aijkl<ExprBljmk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 112


/*!brief A(j,i,k,l)*B(l,k,j,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Akjil_to_Aijkl<ExprBlkjm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBlkjm &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Akjil_to_Aijkl<ExprBlkjm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 113


/*!brief A(j,i,k,l)*B(l,k,m,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Akjli_to_Aijkl<ExprBlkmj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBlkmj &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Akjli_to_Aijkl<ExprBlkmj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 114



/*!brief A(j,i,k,l)*B(l,m,j,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Aklij_to_Aijkl<ExprBlmjk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBlmjk &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Aklij_to_Aijkl<ExprBlmjk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 115


/*!brief A(j,i,k,l)*B(l,m,k,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Aklji_to_Aijkl<ExprBlmkj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBlmkj &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Aklji_to_Aijkl<ExprBlmkj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 116



/*!brief A(j,i,k,l)*B(m,j,k,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Alijk_to_Aijkl<ExprBmjkl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBmjkl &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Alijk_to_Aijkl<ExprBmjkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 117


/*!brief A(j,i,k,l)*B(m,j,l,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Alikj_to_Aijkl<ExprBmjlk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBmjlk &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Alikj_to_Aijkl<ExprBmjlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 118



/*!brief A(j,i,k,l)*B(m,k,j,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Aljik_to_Aijkl<ExprBmkjl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBmkjl &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Aljik_to_Aijkl<ExprBmkjl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 119



/*!brief A(j,i,k,l)*B(m,k,l,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Aljki_to_Aijkl<ExprBmklj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBmklj &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Aljki_to_Aijkl<ExprBmklj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 119 bis



/*!brief A(j,i,k,l)*B(m,l,j,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Alkij_to_Aijkl<ExprBmljk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBmljk &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Alkij_to_Aijkl<ExprBmljk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 120




/*!brief A(j,i,k,l)*B(m,l,k,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajikl_to_Aijkl<ExprAjikl,T>,
            Alkji_to_Aijkl<ExprBmlkj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjikl &a, const ExprBmlkj &b){
    typedef Ajikl_to_Aijkl<ExprAjikl,T> Permuted_Obj1;
    typedef Alkji_to_Aijkl<ExprBmlkj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 122


/*!brief A(j,k,i,l)*B(j,k,l,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            ExprBjklm,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBjklm &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,ExprBjklm,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),b));
}//CHECKED 123


/*!brief A(j,k,i,l)*B(j,k,m,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Aijlk_to_Aijkl<ExprBjkml,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBjkml &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Aijlk_to_Aijkl<ExprBjkml,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 124



/*!brief A(j,k,i,l)*B(j,l,k,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Aikjl_to_Aijkl<ExprBjlkm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBjlkm &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Aikjl_to_Aijkl<ExprBjlkm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 125


/*!brief A(j,k,i,l)*B(j,l,m,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Aiklj_to_Aijkl<ExprBjlmk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBjlmk &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Aiklj_to_Aijkl<ExprBjlmk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 126


/*!brief A(j,k,i,l)*B(j,m,k,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Ailjk_to_Aijkl<ExprBjmkl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBjmkl &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Ailjk_to_Aijkl<ExprBjmkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 127


/*!brief A(j,k,i,l)*B(j,m,l,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Ailkj_to_Aijkl<ExprBjmlk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBjmlk &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Ailkj_to_Aijkl<ExprBjmlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 128



/*!brief A(j,k,i,l)*B(k,j,l,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Ajikl_to_Aijkl<ExprBkjlm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBkjlm &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Ajikl_to_Aijkl<ExprBkjlm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 129


/*!brief A(j,k,i,l)*B(k,j,m,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Ajilk_to_Aijkl<ExprBkjml,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBkjml &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Ajilk_to_Aijkl<ExprBkjml,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 130


/*!brief A(j,k,i,l)*B(k,l,j,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Ajkil_to_Aijkl<ExprBkljm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBkljm &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Ajkil_to_Aijkl<ExprBkljm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 131


/*!brief A(j,k,i,l)*B(k,l,m,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Ajkli_to_Aijkl<ExprBklmj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBklmj &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Ajkli_to_Aijkl<ExprBklmj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 132



/*!brief A(j,k,i,l)*B(k,m,j,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Ajlik_to_Aijkl<ExprBkmjl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBkmjl &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Ajlik_to_Aijkl<ExprBkmjl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 133


/*!brief A(j,k,i,l)*B(k,m,l,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Ajlki_to_Aijkl<ExprBkmlj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBkmlj &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Ajlki_to_Aijkl<ExprBkmlj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 134


/*!brief A(j,k,i,l)*B(l,j,k,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Akijl_to_Aijkl<ExprBljkm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBljkm &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Akijl_to_Aijkl<ExprBljkm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 135


/*!brief A(j,k,i,l)*B(l,j,m,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Akilj_to_Aijkl<ExprBljmk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBljmk &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Akilj_to_Aijkl<ExprBljmk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 136


/*!brief A(j,k,i,l)*B(l,k,j,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Akjil_to_Aijkl<ExprBlkjm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBlkjm &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Akjil_to_Aijkl<ExprBlkjm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 137



/*!brief A(j,k,i,l)*B(l,k,m,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Akjli_to_Aijkl<ExprBlkmj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBlkmj &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Akjli_to_Aijkl<ExprBlkmj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 138




/*!brief A(j,k,i,l)*B(l,m,j,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Aklij_to_Aijkl<ExprBlmjk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBlmjk &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Aklij_to_Aijkl<ExprBlmjk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 139



/*!brief A(j,k,i,l)*B(l,m,k,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Aklji_to_Aijkl<ExprBlmkj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBlmkj &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Aklji_to_Aijkl<ExprBlmkj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 140




/*!brief A(j,k,i,l)*B(m,j,k,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Alijk_to_Aijkl<ExprBmjkl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBmjkl &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Alijk_to_Aijkl<ExprBmjkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 141



/*!brief A(j,k,i,l)*B(m,j,l,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Alikj_to_Aijkl<ExprBmjlk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBmjlk &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Alikj_to_Aijkl<ExprBmjlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 142



/*!brief A(j,k,i,l)*B(m,k,j,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Aljik_to_Aijkl<ExprBmkjl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBmkjl &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Aljik_to_Aijkl<ExprBmkjl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 143


/*!brief A(j,k,i,l)*B(m,k,l,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Aljki_to_Aijkl<ExprBmklj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBmklj &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Aljki_to_Aijkl<ExprBmklj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 144


/*!brief A(j,k,i,l)*B(m,l,j,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Alkij_to_Aijkl<ExprBmljk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBmljk &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Alkij_to_Aijkl<ExprBmljk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 145



/*!brief A(j,k,i,l)*B(m,l,k,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkil_to_Aijkl<ExprAjkil,T>,
            Alkji_to_Aijkl<ExprBmlkj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkil &a, const ExprBmlkj &b){
    typedef Ajkil_to_Aijkl<ExprAjkil,T> Permuted_Obj1;
    typedef Alkji_to_Aijkl<ExprBmlkj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2 ,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 146



/*!brief A(j,k,l,i)*B(j,k,l,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            ExprBjklm,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBjklm &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,ExprBjklm,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),b));
}//CHECKED 147


/*!brief A(j,k,l,i)*B(j,k,m,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Aijlk_to_Aijkl<ExprBjkml,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBjkml &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Aijlk_to_Aijkl<ExprBjkml,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 148

/*!brief A(j,k,l,i)*B(j,l,k,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Aikjl_to_Aijkl<ExprBjlkm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBjlkm &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Aikjl_to_Aijkl<ExprBjlkm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 149


/*!brief A(j,k,l,i)*B(j,l,m,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Aiklj_to_Aijkl<ExprBjlmk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBjlmk &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Aiklj_to_Aijkl<ExprBjlmk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 150



/*!brief A(j,k,l,i)*B(j,m,k,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Ailjk_to_Aijkl<ExprBjmkl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBjmkl &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Ailjk_to_Aijkl<ExprBjmkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 151


/*!brief A(j,k,l,i)*B(j,m,l,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Ailkj_to_Aijkl<ExprBjmlk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBjmlk &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Ailkj_to_Aijkl<ExprBjmlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 152



/*!brief A(j,k,l,i)*B(k,j,l,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Ajikl_to_Aijkl<ExprBkjlm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBkjlm &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Ajikl_to_Aijkl<ExprBkjlm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 153



/*!brief A(j,k,l,i)*B(k,j,m,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Ajilk_to_Aijkl<ExprBkjml,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBkjml &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Ajilk_to_Aijkl<ExprBkjml,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 154




/*!brief A(j,k,l,i)*B(k,l,j,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Ajkil_to_Aijkl<ExprBkljm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBkljm &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Ajkil_to_Aijkl<ExprBkljm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 155



/*!brief A(j,k,l,i)*B(k,l,m,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Ajkli_to_Aijkl<ExprBklmj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBklmj &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Ajkli_to_Aijkl<ExprBklmj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 156




/*!brief A(j,k,l,i)*B(k,m,j,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Ajlik_to_Aijkl<ExprBkmjl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBkmjl &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Ajlik_to_Aijkl<ExprBkmjl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 157



/*!brief A(j,k,l,i)*B(k,m,l,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Ajlki_to_Aijkl<ExprBkmlj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBkmlj &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Ajlki_to_Aijkl<ExprBkmlj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 158


/*!brief A(j,k,l,i)*B(l,j,k,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Akijl_to_Aijkl<ExprBljkm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBljkm &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Akijl_to_Aijkl<ExprBljkm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 159



/*!brief A(j,k,l,i)*B(l,j,m,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Akilj_to_Aijkl<ExprBljmk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBljmk &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Akilj_to_Aijkl<ExprBljmk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 160


/*!brief A(j,k,l,i)*B(l,k,j,m) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Akjil_to_Aijkl<ExprBlkjm,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBlkjm &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Akjil_to_Aijkl<ExprBlkjm,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 161



/*!brief A(j,k,l,i)*B(l,k,m,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Akjli_to_Aijkl<ExprBlkmj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBlkmj &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Akjli_to_Aijkl<ExprBlkmj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 162



/*!brief A(j,k,l,i)*B(l,m,j,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Aklij_to_Aijkl<ExprBlmjk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBlmjk &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Aklij_to_Aijkl<ExprBlmjk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 163



/*!brief A(j,k,l,i)*B(l,m,k,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Aklji_to_Aijkl<ExprBlmkj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBlmkj &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Aklji_to_Aijkl<ExprBlmkj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 164



/*!brief A(j,k,l,i)*B(m,j,k,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Alijk_to_Aijkl<ExprBmjkl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBmjkl &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Alijk_to_Aijkl<ExprBmjkl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 165


/*!brief A(j,k,l,i)*B(m,j,l,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Alikj_to_Aijkl<ExprBmjlk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBmjlk &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Alikj_to_Aijkl<ExprBmjlk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 166

/*!brief A(j,k,l,i)*B(m,k,j,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Aljik_to_Aijkl<ExprBmkjl,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBmkjl &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Aljik_to_Aijkl<ExprBmkjl,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 167



/*!brief A(j,k,l,i)*B(m,k,l,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Aljki_to_Aijkl<ExprBmklj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBmklj &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Aljki_to_Aijkl<ExprBmklj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 168




/*!brief A(j,k,l,i)*B(m,l,j,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Alkij_to_Aijkl<ExprBmljk,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBmljk &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Alkij_to_Aijkl<ExprBmljk,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 169




/*!brief A(j,k,l,i)*B(m,l,k,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l, char m>
inline Expr2<const Aijkl_contracts_Bjklm<
            Ajkli_to_Aijkl<ExprAjkli,T>,
            Alkji_to_Aijkl<ExprBmlkj,U>,
            typename promotedType>,
            typename promotedType, i , m>
operator*(const ExprAjkli &a, const ExprBmlkj &b){
    typedef Ajkli_to_Aijkl<ExprAjkli,T> Permuted_Obj1;
    typedef Alkji_to_Aijkl<ExprBmlkj,U> Permuted_Obj2;
    typedef const Aijkl_contracts_Bjklm<Permuted_Obj1,Permuted_Obj2,typename promotedType> ExprObj;
    return Expr2<ExprObj,typename promotedType,i,m> (ExprObj(Permuted_Obj1(a),Permuted_Obj2(b)));
}//CHECKED 170



/////////////////////////////////////
//FULL INDEX CONTRACTIONS
////////////////////////////////////



/*!brief A(i,j,k,l)*B(i,j,k,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline const typename promotedType
operator*(const ExprAijkl &a, const ExprBijkl &b){

    return Aijkl_contracts_Bijkl<ExprAijkl ,ExprBijkl,typename promotedType> (a,b);

}//CHECKED 171


/*!brief A(i,j,k,l)*B(i,j,l,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBijlk &b){
    typedef Aijlk_to_Aijkl<ExprBijlk,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 172


/*!brief A(i,j,k,l)*B(i,k,j,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBikjl &b){
    typedef Aikjl_to_Aijkl<ExprBikjl,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 173


/*!brief A(i,j,k,l)*B(i,k,l,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBiklj &b){
    typedef Aiklj_to_Aijkl<ExprBiklj,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 174




/*!brief A(i,j,k,l)*B(i,l,j,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBiljk &b){
    typedef Ailjk_to_Aijkl<ExprBiljk,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 175



/*!brief A(i,j,k,l)*B(i,l,k,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBilkj &b){
    typedef Ailkj_to_Aijkl<ExprBilkj,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 176


/*!brief A(i,j,k,l)*B(j,i,k,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBjikl &b){
    typedef Ajikl_to_Aijkl<ExprBjikl,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 177


/*!brief A(i,j,k,l)*B(j,i,l,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBjilk &b){
    typedef Ajilk_to_Aijkl<ExprBjilk,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 178


/*!brief A(i,j,k,l)*B(j,k,i,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBjkil &b){
    typedef Ajkil_to_Aijkl<ExprBjkil,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 179





/*!brief A(i,j,k,l)*B(j,k,l,i) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBjkli &b){
    typedef Ajkli_to_Aijkl<ExprBjkli,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 180




/*!brief A(i,j,k,l)*B(j,l,i,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBjlik &b){
    typedef Ajlik_to_Aijkl<ExprBjlik,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 181



/*!brief A(i,j,k,l)*B(j,l,k,i) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBjlki &b){
    typedef Ajlki_to_Aijkl<ExprBjlki,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 182



/*!brief A(i,j,k,l)*B(k,i,j,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBkijl &b){
    typedef Akijl_to_Aijkl<ExprBkijl,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 183


/*!brief A(i,j,k,l)*B(k,i,l,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBkilj &b){
    typedef Akilj_to_Aijkl<ExprBkilj,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 184



/*!brief A(i,j,k,l)*B(k,j,i,l) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBkjil &b){
    typedef Akjil_to_Aijkl<ExprBkjil,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 185


/*!brief A(i,j,k,l)*B(k,j,l,i) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBkjli &b){
    typedef Akjli_to_Aijkl<ExprBkjli,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 186


/*!brief A(i,j,k,l)*B(k,l,i,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBklij &b){
    typedef Aklij_to_Aijkl<ExprBklij,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 187


/*!brief A(i,j,k,l)*B(k,l,j,i) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBklji &b){
    typedef Aklji_to_Aijkl<ExprBklji,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 188



/*!brief A(i,j,k,l)*B(l,i,j,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBlijk &b){
    typedef Alijk_to_Aijkl<ExprBlijk,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 189


/*!brief A(i,j,k,l)*B(l,i,k,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBlikj &b){
    typedef Alikj_to_Aijkl<ExprBlikj,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 190



/*!brief A(i,j,k,l)*B(l,j,i,k) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBljik &b){
    typedef Aljik_to_Aijkl<ExprBljik,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 191



/*!brief A(i,j,k,l)*B(l,j,k,i) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBljki &b){
    typedef Aljki_to_Aijkl<ExprBljki,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 192


/*!brief A(i,j,k,l)*B(l,k,i,j) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBlkij &b){
    typedef Alkij_to_Aijkl<ExprBlkij,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 193


/*!brief A(i,j,k,l)*B(l,k,j,i) */

template <class A, class B, class T, class U, char i, char j, char k, char l>
inline typename promotedType
operator*(const ExprAijkl &a, const ExprBlkji &b){
    typedef Alkji_to_Aijkl<ExprBlkji,U> Permuted_Obj1;
    return Aijkl_contracts_Bijkl<ExprAijkl,Permuted_Obj1,typename promotedType> (a,Permuted_Obj1(b));

}//CHECKED 194


#undef ExprAijkl
#undef ExprAijlk
#undef ExprAklij
#undef ExprAlkij
#undef ExprAikjl
#undef ExprAiklj
#undef ExprAkijl
#undef ExprAkilj
#undef ExprAjikl
#undef ExprAjkil
#undef ExprAjkli
#undef ExprAiljk
#undef ExprAilkj
#undef ExprAjilk
#undef ExprAjlik
#undef ExprAjlki
#undef ExprAkjil


#undef ExprBklmn
#undef ExprBmnkl
#undef ExprBmnlk
#undef ExprBkmln
#undef ExprBkmnl
#undef ExprBlmkn
#undef ExprBlmnk
#undef ExprBmkln
#undef ExprBmknl
#undef ExprBmlkn
#undef ExprBmlnk
#undef ExprBlkmn
#undef ExprBjklm
#undef ExprBjkml
#undef ExprBjlkm
#undef ExprBjlmk
#undef ExprBjmkl
#undef ExprBjmlk
#undef ExprBkjlm
#undef ExprBkjml
#undef ExprBkljm
#undef ExprBklmj
#undef ExprBkmjl
#undef ExprBkmlj
#undef ExprBljkm
#undef ExprBljmk
#undef ExprBlkjm
#undef ExprBlkmj
#undef ExprBlmjk
#undef ExprBlmkj
#undef ExprBmjkl
#undef ExprBmjlk
#undef ExprBmkjl
#undef ExprBmklj
#undef ExprBmljk
#undef ExprBmlkj
#undef ExprBijkl
#undef ExprBijlk
#undef ExprBikjl
#undef ExprBiklj
#undef ExprBiljk
#undef ExprBilkj
#undef ExprBjikl
#undef ExprBjilk
#undef ExprBjkil
#undef ExprBjkli
#undef ExprBjlik
#undef ExprBjlki
#undef ExprBkijl
#undef ExprBkilj
#undef ExprBkjil
#undef ExprBkjli
#undef ExprBklij
#undef ExprBklji
#undef ExprBlijk
#undef ExprBlikj
#undef ExprBljik
#undef ExprBljki
#undef ExprBlkij
#undef ExprBlkji
#undef promotedType
#endif
