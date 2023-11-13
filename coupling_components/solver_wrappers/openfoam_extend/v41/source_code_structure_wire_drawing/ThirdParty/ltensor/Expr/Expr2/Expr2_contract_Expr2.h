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
#ifndef Expr2_contract_Expr2_H
#define Expr2_contract_Expr2_H


////////////////////////////////////////////////////////////////////
// FULL CONTRACTIONS OF INDiCES
//////////////////////////////////////////////////////////////////



#define ExprAij   Expr2<A,T,i,j>
#define ExprBij   Expr2<B,U,i,j>
#define ExprBji   Expr2<B,U,j,i>
#define promotedType promote<T,U>::V

/*! \brief The operator that defines a contraction
*
*	This operator performs the contraction along indexes (i,j)and returns the value
*	res=A(i,j)*B(i,j)
*/

template < class A, class B, class T, class U, char i, char j >
inline const typename promotedType
operator*  (const ExprAij &ExprL, const ExprBij &ExprR)
{
	return Aij_contracts_Bij<ExprAij,ExprBij,typename promotedType >(ExprL,ExprR);
}

/*! \brief The operator that defines a contraction
*
*	This operator performs the contraction along indexes (i,j)and returns the value
*	res=A(i,j)*B(j,i)
*	\remarks NOTE the permutted indexes
*/

template < class A, class B, class T, class U, char i, char j >
inline const typename promotedType
operator*  (const ExprAij &ExprL, const ExprBji &ExprR)
{
	return Aij_contracts_Bji<ExprAij,ExprBji,typename promotedType >(ExprL,ExprR);
}


#undef ExprAij
#undef ExprBij
#undef ExprBji


////////////////////////////////////////////////////////////////////
// PARTIAL CONTRACTIONS OF INDiCES
////////////////////////////////////////////////////////////////////

#define ExprAik   Expr2<A,T,i,k>
#define ExprAki   Expr2<A,T,k,i>
#define ExprBjk   Expr2<B,U,j,k>
#define ExprBkj   Expr2<B,U,k,j>


//
//
//



/*! \brief The operator that defines a contraction
*
*	This operator performs the contraction along index k and returns the value as Expr of rank 2
*	A(i,k)*B(k,j)==>C(i,j)
*	
*/

template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr2 < const Aik_contracts_Bkj < ExprAik , ExprBkj , typename promotedType>,
				typename promotedType, i, j >
operator* (const ExprAik &a, const ExprBkj &b)
{
	typedef const Aik_contracts_Bkj < ExprAik , ExprBkj , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , i, j> (ExprObj (a,b));
}

//
//
//

/*! \brief The operator that defines a contraction
*
*	This operator performs the contraction along index k and returns the value as Expr of rank 2
*	A(k,i)*B(j,k)==>C(i,j)	
*/



template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr2 < const Aki_contracts_Bjk < ExprAki , ExprBjk , typename promotedType>,
				typename promotedType, i, j >
operator* (const ExprAki &a, const ExprBjk &b)
{
	typedef const Aki_contracts_Bjk < ExprAki , ExprBjk , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , i, j> (ExprObj (a,b));
}



/*! \brief The operator that defines a contraction
*
*	This operator performs the contraction along index k and returns the value as Expr of rank 2
*	A(i,k)*B(j,k)==>C(i,j)	
*/

template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr2 < const Aik_contracts_Bjk < ExprAik , ExprBjk , typename promotedType>,
				typename promotedType, i, j >
operator* (const ExprAik &a, const ExprBjk &b)
{
	typedef const Aik_contracts_Bjk < ExprAik , ExprBjk , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , i, j> (ExprObj (a,b));
}


/*! \brief The operator that defines a contraction
*
*	This operator performs the contraction along index k and returns the value as Expr of rank 2
*	A(k,i)*B(k,j)==>C(i,j)	
*/


template < class A, class B, class T, class U, char i , char j, char k>
inline const Expr2 < const Aki_contracts_Bkj < ExprAki , ExprBkj , typename promotedType>,
				typename promotedType, i, j >
operator* (const ExprAki &a, const ExprBkj &b)
{
	typedef const Aki_contracts_Bkj < ExprAki , ExprBkj , typename promotedType> ExprObj;
	return Expr2 < ExprObj,typename promotedType , i, j> (ExprObj (a,b));
}


#undef ExprAik
#undef ExprAki
#undef ExprBjk
#undef ExprBkj
#undef pType

#endif
