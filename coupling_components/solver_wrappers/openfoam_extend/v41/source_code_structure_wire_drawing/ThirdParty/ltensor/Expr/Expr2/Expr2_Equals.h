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
#ifndef Expr2_Equals_H
#define Expr2_Equals_H

//#define CHECK_Expr2_Definitions
//#define USE_ASSERT_Expr2
#ifdef USE_ASSERT_Expr2
#include <assert.h>
#endif

//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// T(i,j) = C(i,j)
//////////////////////////////////////////////////////////////////////

#define promotedType promote<T,U>::V

template<class A, class T, char i, char j>
template<class B,class U>
inline const Expr2<A,T,i,j> &
Expr2<A,T,i,j>::operator=(const Expr2<B,U,i,j> &rhs)
{
	Lij_equals_Rij((*this),rhs);
	return *this;
}




template<class A, class T, char i, char j>
inline const Expr2<A,T,i,j> &
Expr2<A,T,i,j>::operator=(const Expr2<A,T,i,j> &rhs)
{
	Lij_equals_Rij((*this),rhs);
	return *this;
}


template<class A, class T, char i, char j>
template<class B,class U>
inline const Expr2<A,T,i,j> &
Expr2<A,T,i,j>::operator+=(const Expr2<B,U,i,j> &rhs)
{
	Lij_plusequals_Rij((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j>
template<class B,class U>
inline const Expr2<A,T,i,j> &
Expr2<A,T,i,j>::operator-=(const Expr2<B,U,i,j> &rhs)
{
	Lij_minusequals_Rij((*this),rhs);
	return *this;
}

///////////////////////////////////////////////////////
// Complement of above operations...
// T(i,j) = B(j,i)....
///////////////////////////////////////////////////////


template<class A, class T, char i, char j>
template<class B,class U>
inline const Expr2<A,T,i,j> &
Expr2<A,T,i,j>::operator=(const Expr2<B,U,j,i> &rhs)
{
	Lij_equals_Rji((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j>
inline const Expr2<A,T,i,j> &
Expr2<A,T,i,j>::operator=(const Expr2<A,T,j,i> &rhs)
{
	Lij_equals_Rji((*this),rhs);
	return *this;
}


template<class A, class T, char i, char j>
template<class B,class U>
inline const Expr2<A,T,i,j> &
Expr2<A,T,i,j>::operator+=(const Expr2<B,U,j,i> &rhs)
{
	Lij_plusequals_Rji((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j>
template<class B,class U>
inline const Expr2<A,T,i,j> &
Expr2<A,T,i,j>::operator-=(const Expr2<B,U,j,i> &rhs)
{
	Lij_minusequals_Rji((*this),rhs);
	return *this;
}



#endif
