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
#ifndef Expr3_contract_Expr3_H
#define Expr3_contract_Expr3_H


////////////////////////////////////////////////////////////////////
// FULL CONTRACTIONS OF INDiCES
//////////////////////////////////////////////////////////////////


#define ExprAijk   Expr3<A,T,i,j,k>
#define ExprAikj   Expr3<A,T,i,k,j>
#define ExprAjik   Expr3<A,T,j,i,k>
#define ExprAkij   Expr3<A,T,k,i,j>
#define ExprAjki   Expr3<A,T,j,k,i>
#define ExprAkji   Expr3<A,T,k,j,i>

#define ExprBijk   Expr3<B,U,i,j,k>
#define ExprBilj   Expr3<B,U,i,l,j>

#define promotedType promote<T,U>::V


//
/*!\ brief A(i,j,k)*B(i,j,k) */
//
template < class A, class B, class T, class U, char i, char j, char k >
inline const typename promotedType
operator*  (const ExprAijk &ExprL, const ExprBijk &ExprR)
{
	return Aijk_contracts_Bijk<ExprAijk,ExprBijk,typename promotedType >(ExprL,ExprR);
}

//
/*!\ brief A(i,k,j)*B(i,j,k) */
//

template < class A, class B, class T, class U, char i, char j, char k >
inline const typename promotedType
operator*  (const ExprAikj &ExprL, const ExprBijk &ExprR)
{
	typedef Aikj_to_Aijk < ExprAikj ,T > Permuted_Obj;
	return Aijk_contracts_Bijk<Permuted_Obj,ExprBijk,typename promotedType >(Permuted_Obj(ExprL),ExprR);
}

//
/*!\ brief A(j,i,k)*B(i,j,k) */
//
template < class A, class B, class T, class U, char i, char j, char k >
inline const typename promotedType
operator*  (const ExprAjik &ExprL, const ExprBijk &ExprR)
{
	typedef Ajik_to_Aijk < ExprAjik ,T > Permuted_Obj;
	return Aijk_contracts_Bijk<Permuted_Obj,ExprBijk,typename promotedType >(Permuted_Obj(ExprL),ExprR);
}

//
/*!\ brief A(k,i,j)*B(i,j,k) */
//
template < class A, class B, class T, class U, char i, char j, char k >
inline const typename promotedType
operator*  (const ExprAkij &ExprL, const ExprBijk &ExprR)
{
	typedef Akij_to_Aijk < ExprAkij ,T > Permuted_Obj;
	return Aijk_contracts_Bijk<Permuted_Obj,ExprBijk,typename promotedType >(Permuted_Obj(ExprL),ExprR);
}

//
/*!\ brief A(j,k,i)*B(i,j,k) */
//
template < class A, class B, class T, class U, char i, char j, char k >
inline const typename promotedType
operator*  (const ExprAjki &ExprL, const ExprBijk &ExprR)
{
	typedef Ajki_to_Aijk < ExprAjki ,T > Permuted_Obj;
	return Aijk_contracts_Bijk<Permuted_Obj,ExprBijk,typename promotedType >(Permuted_Obj(ExprL),ExprR);
}

//
/*!\ brief A(k,j,i)*B(i,j,k) */
//
template < class A, class B, class T, class U, char i, char j, char k >
inline const typename promotedType
operator*  (const ExprAkji &ExprL, const ExprBijk &ExprR)
{
	typedef Akji_to_Aijk < ExprAkji ,T > Permuted_Obj;
	return Aijk_contracts_Bijk<Permuted_Obj,ExprBijk,typename promotedType >(Permuted_Obj(ExprL),ExprR);
}


#undef ExprAijk
#undef ExprAikj
#undef ExprAjik
#undef ExprAkij
#undef ExprAjki
#undef ExprAkji

#undef ExprBijk



////////////////////////////////////////////////////////////////////
// SINGLE CONTRACTIONS OF INDiCES
//

//
////////////////////////////////////////////////////////////////////



#define ExprAijk   Expr3<A,T,i,j,k>
#define ExprAikj   Expr3<A,T,i,k,j>
#define ExprAjik   Expr3<A,T,j,i,k>
#define ExprAkij   Expr3<A,T,k,i,j>
#define ExprAjki   Expr3<A,T,j,k,i>
#define ExprAkji   Expr3<A,T,k,j,i>

#define ExprBilm   Expr3<B,U,i,l,m>
#define ExprBlim   Expr3<B,U,l,i,m>
#define ExprBlmi   Expr3<B,U,l,m,i>

/*!\ brief A(i,j,k)*B(i,l,m) ----> C(j,k,m,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijk_contracts_Biml < ExprAijk , Aikj_to_Aijk<ExprBilm,U> , typename promotedType>,
					typename promotedType, j,k,m,l >
operator* (const ExprAijk &a, const ExprBilm &b)
{

typedef Aikj_to_Aijk<ExprBilm,U> Permuted_Obj1;
typedef const Aijk_contracts_Biml < ExprAijk , Permuted_Obj1 , typename promotedType> ExprObj;

return Expr4 < const ExprObj,typename promotedType , j,k,m,l> (ExprObj (a,Permuted_Obj1(b)));
}


/*!\ brief A(i,j,k)*B(l,i,m) ----> C(j,k,m,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijk_contracts_Biml < ExprAijk , Akij_to_Aijk<ExprBlim,U> , typename promotedType>,
					typename promotedType, j,k,m,l >
operator* (const ExprAijk &a, const ExprBlim &b)
{
typedef Akij_to_Aijk<ExprBlim,U> Permuted_Obj1;
typedef const Aijk_contracts_Biml < ExprAijk , Permuted_Obj1, typename promotedType> ExprObj;
return Expr4 < ExprObj,typename promotedType , j,k,m,l> (ExprObj (a,Permuted_Obj1(b)));
}


/*!\ brief A(i,j,k)*B(l,m,i) ----> C(j,k,m,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijk_contracts_Biml < ExprAijk , Akji_to_Aijk<ExprBlmi,U> , typename promotedType>,
					typename promotedType, j,k,m,l >
operator* (const ExprAijk &a, const ExprBlmi &b)
{
typedef Akji_to_Aijk<ExprBlmi,U> Permuted_Obj1;
typedef const Aijk_contracts_Biml < ExprAijk , Permuted_Obj1, typename promotedType> ExprObj;
return Expr4 < ExprObj,typename promotedType , j,k,m,l> (ExprObj (a,Permuted_Obj1(b)));
}



/*!\ brief A(j,i,k)*B(i,l,m) ----> C(j,k,m,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijk_contracts_Biml < Ajik_to_Aijk<ExprAjik,T> , Aikj_to_Aijk<ExprBilm,U> , typename promotedType>,
					typename promotedType, j,k,m,l>
operator* (const ExprAjik &a, const ExprBilm &b)
{
typedef  Ajik_to_Aijk<ExprAjik,T> Permuted_Obj1;
typedef  Aikj_to_Aijk<ExprBilm,U> Permuted_Obj2;
typedef const Aijk_contracts_Biml < Permuted_Obj1, Permuted_Obj2, typename promotedType> ExprObj;
return Expr4 < ExprObj,typename promotedType , j,k,m,l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!\ brief A(j,i,k)*B(l,i,m) ----> C(j,k,m,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijk_contracts_Biml < Ajik_to_Aijk<ExprAjik,T> , Akij_to_Aijk<ExprBlim,U> , typename promotedType>,
					typename promotedType, j,k,m,l>
operator* (const ExprAjik &a, const ExprBlim &b)
{
typedef Ajik_to_Aijk<ExprAjik,T> Permuted_Obj1;
typedef Akij_to_Aijk<ExprBlim,U> Permuted_Obj2;
typedef const Aijk_contracts_Biml < Permuted_Obj1, Permuted_Obj2, typename promotedType> ExprObj;
return Expr4 < ExprObj,typename promotedType , j,k,m,l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


/*!\ brief A(j,i,k)*B(l,m,i) ----> C(j,k,m,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijk_contracts_Biml < Ajik_to_Aijk<ExprAjik,T> , Akji_to_Aijk<ExprBlmi,U> , typename promotedType>,
					typename promotedType, j,k,m,l >
operator* (const ExprAjik &a, const ExprBlmi &b)
{
typedef Ajik_to_Aijk<ExprAjik,T> Permuted_Obj1;
typedef Akji_to_Aijk<ExprBlmi,U> Permuted_Obj2;
typedef const Aijk_contracts_Biml < Permuted_Obj1, Permuted_Obj2, typename promotedType> ExprObj;
return Expr4 < ExprObj,typename promotedType , j,k,m,l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!\ brief A(j,k,i)*B(i,l,m) ----> C(j,k,m,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijk_contracts_Biml < Ajki_to_Aijk<ExprAjki,T> , Aikj_to_Aijk<ExprBilm,U> , typename promotedType>,
					typename promotedType, j,k,m,l >
operator* (const ExprAjki &a, const ExprBilm &b)
{
typedef Ajki_to_Aijk<ExprAjki,T> Permuted_Obj1;
typedef Aikj_to_Aijk<ExprBilm,U> Permuted_Obj2;
typedef const Aijk_contracts_Biml < Permuted_Obj1, Permuted_Obj2, typename promotedType> ExprObj;
return Expr4 < ExprObj,typename promotedType , j,k,m,l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}




/*!\ brief A(j,k,i)*B(l,i,m) ----> C(j,k,m,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijk_contracts_Biml < Ajki_to_Aijk<ExprAjki,T> , Akij_to_Aijk<ExprBlim,U> , typename promotedType>,
					typename promotedType, j,k,m,l >
operator* (const ExprAjki &a, const ExprBlim &b)
{
typedef Ajki_to_Aijk<ExprAjki,T> Permuted_Obj1;
typedef Akij_to_Aijk<ExprBlim,U> Permuted_Obj2;
typedef const Aijk_contracts_Biml < Permuted_Obj1, Permuted_Obj2, typename promotedType> ExprObj;
return Expr4 < ExprObj,typename promotedType , j,k,m,l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



/*!\ brief A(j,k,i)*B(l,m,i) ----> C(j,k,m,l) */

template < class A, class B, class T, class U, char i , char j, char k, char l, char m>
inline const Expr4 < const Aijk_contracts_Biml < Ajki_to_Aijk<ExprAjki,T> , Akji_to_Aijk<ExprBlmi,U> , typename promotedType>,
					typename promotedType, j,k,m,l >
operator* (const ExprAjki &a, const ExprBlmi &b)
{
typedef Ajki_to_Aijk<ExprAjki,T> Permuted_Obj1;
typedef Akji_to_Aijk<ExprBlmi,U> Permuted_Obj2;
typedef const Aijk_contracts_Biml < Permuted_Obj1, Permuted_Obj2, typename promotedType> ExprObj;
return Expr4 < ExprObj,typename promotedType , j,k,m,l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


#undef ExprBilm
#undef ExprBlim
#undef ExprBlmi



#undef ExprAijk
#undef ExprAikj
#undef ExprAjik
#undef ExprAkij
#undef ExprAjki
#undef ExprAkji



////////////////////////////////////////////////////////////////////
// DOUBLE CONTRACTIONS OF INDiCES
//
/*!\ brief A(i,j,k)*B(i,j,l)   ----> C(k,l) */
//
////////////////////////////////////////////////////////////////////

#define ExprAijk   Expr3<A,T,i,j,k>
#define ExprAjik   Expr3<A,T,j,i,k>
#define ExprAkij   Expr3<A,T,k,i,j>
#define ExprAkji   Expr3<A,T,k,j,i>
#define ExprAikj   Expr3<A,T,i,k,j>

#define ExprBijl   Expr3<B,U,i,j,l>
#define ExprBlij   Expr3<B,U,l,i,j>
#define ExprBjli   Expr3<B,U,j,l,i>
#define ExprBjil   Expr3<B,U,j,i,l>
#define ExprBlji   Expr3<B,U,l,j,i>

#define promotedType promote<T,U>::V

//
/*!\ brief A(i,j,k)*B(i,j,l)   ----> C(k,l) */
//

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl < ExprAijk , ExprBijl , typename promotedType>,
				typename promotedType, k, l >
operator* (const ExprAijk &a, const ExprBijl &b)
{
	typedef const Aijk_contracts_Bijl < ExprAijk , ExprBijl , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (a,b));
}
//checked

//
/*!\ brief A(j,i,k)*B(i,j,l)   ----> C(k,l) */
//

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl< Ajik_to_Aijk < ExprAjik ,T > ,
					ExprBijl , typename promotedType>, typename promotedType, k, l >
operator* (const ExprAjik &a, const ExprBijl &b)
{
	typedef Ajik_to_Aijk < ExprAjik ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bijl < Permuted_Obj , ExprBijl , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj(a),b));
}

//
/*!\ brief A(k,i,j)*B(i,j,l)   ----> C(k,l) */
//

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl< Akij_to_Aijk < ExprAkij ,T > ,
					ExprBijl , typename promotedType>, typename promotedType, k, l >
operator* (const ExprAkij &a, const ExprBijl &b)
{
	typedef Akij_to_Aijk < ExprAkij ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bijl < Permuted_Obj , ExprBijl , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj(a),b));
}

//
/*!\ brief A(k,j,i)*B(i,j,l)   ----> C(k,l) */
//

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl< Akji_to_Aijk < ExprAkji ,T > ,
					ExprBijl , typename promotedType>, typename promotedType, k, l >
operator* (const ExprAkji &a, const ExprBijl &b)
{
	typedef Akji_to_Aijk < ExprAkji ,T > Permuted_Obj;
	typedef const Aijk_contracts_Bijl < Permuted_Obj , ExprBijl , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj(a),b));
}

//
/*!\ brief A(k,i,j)*B(l,i,j)   ----> C(k,l) */
//

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl<
				Akij_to_Aijk < ExprAkij ,T > ,
				Akij_to_Aijk < ExprBlij ,U >	,
				typename promotedType>, typename promotedType, k, l >
operator* (const ExprAkij &a, const ExprBlij &b)
{
	typedef Akij_to_Aijk < ExprAkij ,T > Permuted_Obj1;
	typedef Akij_to_Aijk < ExprBlij ,U > Permuted_Obj2;

	typedef const Aijk_contracts_Bijl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}

//
/*!\ brief A(k,j,i)*B(l,i,j)   ----> C(k,l) */
//

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl<
				Akji_to_Aijk < ExprAkji ,T > ,
				Akij_to_Aijk < ExprBlij ,U >	,
				typename promotedType>, typename promotedType, k, l >
operator* (const ExprAkji &a, const ExprBlij &b)
{
	typedef Akji_to_Aijk < ExprAkji ,T > Permuted_Obj1;
	typedef Akij_to_Aijk < ExprBlij ,U > Permuted_Obj2;

	typedef const Aijk_contracts_Bijl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}





//
/*!\ brief A(k,j,i)*B(i,l,j)   ----> C(k,l) */
//

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl<
				Akji_to_Aijk < ExprAkji ,T > ,
				Aikj_to_Aijk < ExprBilj ,U >	,
				typename promotedType>, typename promotedType, k, l >
operator* (const ExprAkji &a, const ExprBilj &b)
{
	typedef Akji_to_Aijk < ExprAkji ,T > Permuted_Obj1;
	typedef Aikj_to_Aijk < ExprBilj ,U > Permuted_Obj2;

	typedef const Aijk_contracts_Bijl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}




//
/*!\ brief A(i,j,k)*B(l,i,j)   ----> C(k,l) */
//

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl<
				Aijk_to_Aijk < ExprAijk ,T > ,
				Akij_to_Aijk < ExprBlij ,U >	,
				typename promotedType>, typename promotedType, k, l >
operator* (const ExprAijk &a, const ExprBlij &b)
{
	typedef Aijk_to_Aijk < ExprAijk ,T > Permuted_Obj1;
	typedef Akij_to_Aijk < ExprBlij ,U > Permuted_Obj2;

	typedef const Aijk_contracts_Bijl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}

//
/*!\ brief A(j,i,k)*B(l,i,j)   ----> C(k,l) */
//

template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl<
				Ajik_to_Aijk < ExprAjik ,T > ,
				Akij_to_Aijk < ExprBlij ,U >	,
				typename promotedType>, typename promotedType, k, l >
operator* (const ExprAjik &a, const ExprBlij &b)
{
	typedef Ajik_to_Aijk < ExprAjik ,T > Permuted_Obj1;
	typedef Akij_to_Aijk < ExprBlij ,U > Permuted_Obj2;

	typedef const Aijk_contracts_Bijl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}

//added 4/9/08

//
/*!\ brief A(i,j,k)*B(i,l,j)   ----> C(k,l) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl<
				ExprAijk ,
				Aikj_to_Aijk < ExprBilj ,U >	,
				typename promotedType>, typename promotedType, k, l >
operator* (const ExprAijk &a, const ExprBilj &b)
{
	typedef Aikj_to_Aijk < ExprBilj ,U > Permuted_Obj2;

	typedef const Aijk_contracts_Bijl < ExprAijk , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (a,Permuted_Obj2(b)));
}



//
/*!\ brief A(i,j,k)*B(j,l,i)   ----> C(k,l) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl<
				ExprAijk ,
				Ajki_to_Aijk < ExprBjli ,U >	,
				typename promotedType>, typename promotedType, k, l >
operator* (const ExprAijk &a, const ExprBjli &b)
{
	typedef Ajki_to_Aijk < ExprBjli ,U > Permuted_Obj2;

	typedef const Aijk_contracts_Bijl < ExprAijk , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (a,Permuted_Obj2(b)));
}



//
/*!\ brief A(i,k,j)*B(i,j,l)   ----> C(k,l) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl<
				Aikj_to_Aijk< ExprAikj , T>  ,
				ExprBijl ,
				typename promotedType>, typename promotedType, k, l >
operator* (const ExprAikj &a, const ExprBijl &b)
{
	typedef Aikj_to_Aijk < ExprAikj ,T > Permuted_Obj1;

	typedef const Aijk_contracts_Bijl < Permuted_Obj1 , ExprBijl, typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj1(a),b));
}


//
/*!\ brief A(i,k,j)*B(i,l,j)   ----> C(k,l) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl<
				Aikj_to_Aijk< ExprAikj , T>  ,
				Aikj_to_Aijk< ExprBilj , U>  ,
				typename promotedType>, typename promotedType, k, l >
operator* (const ExprAikj &a, const ExprBilj &b)
{
	typedef Aikj_to_Aijk < ExprAikj ,T > Permuted_Obj1;
	typedef Aikj_to_Aijk < ExprBilj ,U > Permuted_Obj2;

	typedef const Aijk_contracts_Bijl < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



//
/*!\ brief A(i,k,j)*B(j,i,l)   ----> C(k,l) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl<
				Aikj_to_Aijk< ExprAikj , T>  ,
				Ajik_to_Aijk< ExprBjil , U>  ,
				typename promotedType>, typename promotedType, k, l >
operator* (const ExprAikj &a, const ExprBjil &b)
{
	typedef Aikj_to_Aijk < ExprAikj ,T > Permuted_Obj1;
	typedef Ajik_to_Aijk< ExprBjil , U> Permuted_Obj2;

	typedef const Aijk_contracts_Bijl < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


//
/*!\ brief A(i,k,j)*B(j,l,i)   ----> C(k,l) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl<
				Aikj_to_Aijk< ExprAikj , T>  ,
				Ajki_to_Aijk< ExprBjli , U>  ,
				typename promotedType>, typename promotedType, k, l >
operator* (const ExprAikj &a, const ExprBjli &b)
{
	typedef Aikj_to_Aijk < ExprAikj ,T > Permuted_Obj1;
	typedef Ajki_to_Aijk< ExprBjli , U> Permuted_Obj2;

	typedef const Aijk_contracts_Bijl < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



//
/*!\ brief A(i,k,j)*B(l,i,j)   ----> C(k,l) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl<
				Aikj_to_Aijk< ExprAikj , T>  ,
				Akij_to_Aijk< ExprBlij , U>  ,
				typename promotedType>, typename promotedType, k, l >
operator* (const ExprAikj &a, const ExprBlij &b)
{
	typedef Aikj_to_Aijk < ExprAikj ,T > Permuted_Obj1;
	typedef Akij_to_Aijk< ExprBlij , U> Permuted_Obj2;

	typedef const Aijk_contracts_Bijl < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



//
/*!\ brief A(i,k,j)*B(l,j,i)   ----> C(k,l) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl<
				Aikj_to_Aijk< ExprAikj , T>  ,
				Akji_to_Aijk< ExprBlji , U>  ,
				typename promotedType>, typename promotedType, k, l >
operator* (const ExprAikj &a, const ExprBlji &b)
{
	typedef Aikj_to_Aijk < ExprAikj ,T > Permuted_Obj1;
	typedef Akji_to_Aijk< ExprBlji , U> Permuted_Obj2;

	typedef const Aijk_contracts_Bijl < Permuted_Obj1 , Permuted_Obj2 , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}



//
/*!\ brief A(k,i,j)*B(i,l,j)   ----> C(k,l) */
//
template < class A, class B, class T, class U, char i , char j, char k, char l>
inline const Expr2 < const Aijk_contracts_Bijl<
				Akij_to_Aijk< ExprAkij , T>  ,
				Aikj_to_Aijk< ExprBilj , U>  ,
				typename promotedType>, typename promotedType, k, l >
operator* (const ExprAkij &a, const ExprBilj &b)
{
	typedef Akij_to_Aijk < ExprAkij ,T > Permuted_Obj1;
	typedef Aikj_to_Aijk< ExprBilj , U> Permuted_Obj2;


	typedef const Aijk_contracts_Bijl < Permuted_Obj1 , Permuted_Obj2, typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , k, l> (ExprObj (Permuted_Obj1(a),Permuted_Obj2(b)));
}


#undef ExprAijk
#undef ExprAjik
#undef ExprAkij
#undef ExprAkji
#undef ExprAikj
#undef ExprBjil


#undef ExprBijl
#undef ExprBlij
#undef ExprBilj
#undef ExprBjli
#undef ExprBlji
#undef promotedType

#endif
