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
#ifndef Expr3_Equals_H
#define Expr3_Equals_H

//#define CHECK_Expr3_Definitions
//#define USE_ASSERT_Expr3
#ifdef USE_ASSERT_Expr3
#include <assert.h>
#endif



// T(i,j,k) = B(i,j,k)....
//

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator=(const Expr3<B,T,i,j,k> &rhs)
{
	Lijk_equals_Rijk((*this),rhs);
	return *this;
}



template<class A, class T, char i, char j, char k>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator=(const Expr3<A,T,i,j,k> &rhs)
{
	Lijk_equals_Rijk((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator+=(const Expr3<B,T,i,j,k> &rhs)
{
	Lijk_plusequals_Rijk((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator-=(const Expr3<B,T,i,j,k> &rhs)
{
	Lijk_minusequals_Rijk((*this),rhs);
	return *this;
}

//
// T(i,j,k) = B(i,k,j)....
//

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator=(const Expr3<B,T,i,k,j> &rhs)
{
	Lijk_equals_Rikj((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator=(const Expr3<A,T,i,k,j> &rhs)
{
	Lijk_equals_Rikj((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator+=(const Expr3<B,T,i,k,j> &rhs)
{
	Lijk_plusequals_Rikj((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator-=(const Expr3<B,T,i,k,j> &rhs)
{
	Lijk_minusequals_Rikj((*this),rhs);
	return *this;
}

//
// T(i,j,k) = B(j,i,k)....
//

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator=(const Expr3<B,T,j,i,k> &rhs)
{
	Lijk_equals_Rjik((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator=(const Expr3<A,T,j,i,k> &rhs)
{
	Lijk_equals_Rjik((*this),rhs);
	return *this;
}


template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator+=(const Expr3<B,T,j,i,k> &rhs)
{
	Lijk_plusequals_Rjik((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator-=(const Expr3<B,T,j,i,k> &rhs)
{
	Lijk_minusequals_Rjik((*this),rhs);
	return *this;
}

//
// T(i,j,k) = B(k,i,j)....
//

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator=(const Expr3<B,T,k,i,j> &rhs)
{
	Lijk_equals_Rkij((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator=(const Expr3<A,T,k,i,j> &rhs)
{
	Lijk_equals_Rkij((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator+=(const Expr3<B,T,k,i,j> &rhs)
{
	Lijk_plusequals_Rkij((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator-=(const Expr3<B,T,k,i,j> &rhs)
{
	Lijk_minusequals_Rkij((*this),rhs);
	return *this;
}

//
// T(i,j,k) = B(j,k,i)....
//

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator=(const Expr3<B,T,j,k,i> &rhs)
{
	Lijk_equals_Rjki((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator=(const Expr3<A,T,j,k,i> &rhs)
{
	Lijk_equals_Rjki((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator+=(const Expr3<B,T,j,k,i> &rhs)
{
	Lijk_plusequals_Rjki((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator-=(const Expr3<B,T,j,k,i> &rhs)
{
	Lijk_minusequals_Rjki((*this),rhs);
	return *this;
}

//
// T(i,j,k) = B(k,j,i)....
//

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator=(const Expr3<B,T,k,j,i> &rhs)
{
	Lijk_equals_Rkji((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator=(const Expr3<A,T,k,j,i> &rhs)
{
	Lijk_equals_Rkji((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator+=(const Expr3<B,T,k,j,i> &rhs)
{
	Lijk_plusequals_Rkji((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k>
template<class B>
inline const Expr3<A,T,i,j,k> &
Expr3<A,T,i,j,k>::operator-=(const Expr3<B,T,k,j,i> &rhs)
{
	Lijk_minusequals_Rkji((*this),rhs);
	return *this;
}



#endif
