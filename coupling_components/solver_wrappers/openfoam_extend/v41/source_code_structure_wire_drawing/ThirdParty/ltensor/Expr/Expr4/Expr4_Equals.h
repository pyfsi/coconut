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
#ifndef Expr4_Equals_H
#define Expr4_Equals_H



// L(i,j,k,l) = R(i,j,k,l)....00

template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,i,j,k,l> &rhs)
{
	Lijkl_equals_Rijkl((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,i,j,k,l> &rhs)
{
	Lijkl_equals_Rijkl((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,i,j,k,l> &rhs)
{
	Lijkl_plusequals_Rijkl((*this),rhs);
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,i,j,k,l> &rhs)
{
	Lijkl_minusequals_Rijkl((*this),rhs);
	return *this;
}

// L(i,j,k,l) = R(j,i,k,l)....01
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,j,i,k,l> &rhs)
{
	typedef Ajikl_to_Aijkl<Expr4<A,T,j,i,k,l>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,j,i,k,l> &rhs)
{
	typedef Ajikl_to_Aijkl<Expr4<B,U,j,i,k,l>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,j,i,k,l> &rhs)
{
	typedef Ajikl_to_Aijkl<Expr4<B,U,j,i,k,l>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,j,i,k,l> &rhs)
{
	typedef Ajikl_to_Aijkl<Expr4<B,U,j,i,k,l>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

// L(i,j,k,l) = R(i,j,l,k)....03
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,i,j,l,k> &rhs)
{
	typedef Aijlk_to_Aijkl<Expr4<A,T,i,j,l,k>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,i,j,l,k> &rhs)
{
	typedef Aijlk_to_Aijkl<Expr4<B,U,i,j,l,k>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,i,j,l,k> &rhs)
{
	typedef Aijlk_to_Aijkl<Expr4<B,U,i,j,l,k>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,i,j,l,k> &rhs)
{
	typedef Aijlk_to_Aijkl<Expr4<B,U,i,j,l,k>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

// L(i,j,k,l) = R(j,i,l,k)....04
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,j,i,l,k> &rhs)
{
	typedef Ajilk_to_Aijkl<Expr4<A,T,j,i,l,k>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,j,i,l,k> &rhs)
{
	typedef Ajilk_to_Aijkl<Expr4<B,U,j,i,l,k>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,j,i,l,k> &rhs)
{
	typedef Ajilk_to_Aijkl<Expr4<B,U,j,i,l,k>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,j,i,l,k> &rhs)
{
	typedef Ajilk_to_Aijkl<Expr4<B,U,j,i,l,k>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}


// L(i,j,k,l) = R(k,l,i,j)....05
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,k,l,i,j> &rhs)
{
	typedef Aklij_to_Aijkl<Expr4<A,T,k,l,i,j>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,k,l,i,j> &rhs)
{
	typedef Aklij_to_Aijkl<Expr4<B,U,k,l,i,j>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,k,l,i,j> &rhs)
{
	typedef Aklij_to_Aijkl<Expr4<B,U,k,l,i,j>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,k,l,i,j> &rhs)
{
	typedef Aklij_to_Aijkl<Expr4<B,U,k,l,i,j>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

// L(i,j,k,l) = R(l,k,i,j)....06
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,l,k,i,j> &rhs)
{
	typedef Alkij_to_Aijkl<Expr4<A,T,l,k,i,j>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,l,k,i,j> &rhs)
{
	typedef Alkij_to_Aijkl<Expr4<B,U,l,k,i,j>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,l,k,i,j> &rhs)
{
	typedef Alkij_to_Aijkl<Expr4<B,U,l,k,i,j>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,l,k,i,j> &rhs)
{
	typedef Alkij_to_Aijkl<Expr4<B,U,l,k,i,j>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

// L(i,j,k,l) = R(l,k,j,i)....07
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,l,k,j,i> &rhs)
{
	typedef Alkji_to_Aijkl<Expr4<A,T,l,k,j,i>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,l,k,j,i> &rhs)
{
	typedef Alkji_to_Aijkl<Expr4<B,U,l,k,j,i>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,l,k,j,i> &rhs)
{
	typedef Alkji_to_Aijkl<Expr4<B,U,l,k,j,i>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,l,k,j,i> &rhs)
{
	typedef Alkji_to_Aijkl<Expr4<B,U,l,k,j,i>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

// L(i,j,k,l) = R(k,l,j,i)....08
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,k,l,j,i> &rhs)
{
	typedef Aklji_to_Aijkl<Expr4<A,T,k,l,j,i>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,k,l,j,i> &rhs)
{
	typedef Aklji_to_Aijkl<Expr4<B,U,k,l,j,i>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,k,l,j,i> &rhs)
{
	typedef Aklji_to_Aijkl<Expr4<B,U,k,l,j,i>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,k,l,j,i> &rhs)
{
	typedef Aklji_to_Aijkl<Expr4<B,U,k,l,j,i>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

// L(i,j,k,l) = R(l,j,k,i)....10
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,l,j,k,i> &rhs)
{
	typedef Aljki_to_Aijkl<Expr4<A,T,l,j,k,i>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,l,j,k,i> &rhs)
{
	typedef Aljki_to_Aijkl<Expr4<B,U,l,j,k,i>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,l,j,k,i> &rhs)
{
	typedef Aljki_to_Aijkl<Expr4<B,U,l,j,k,i>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,l,j,k,i> &rhs)
{
	typedef Aljki_to_Aijkl<Expr4<B,U,l,j,k,i>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}



// L(i,j,k,l) = R(i,k,j,l)....11
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,i,k,j,l> &rhs)
{
	typedef Aikjl_to_Aijkl<Expr4<A,T,i,k,j,l>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,i,k,j,l> &rhs)
{
	typedef Aikjl_to_Aijkl<Expr4<B,U,i,k,j,l>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,i,k,j,l> &rhs)
{
	typedef Aikjl_to_Aijkl<Expr4<B,U,i,k,j,l>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,i,k,j,l> &rhs)
{
	typedef Aikjl_to_Aijkl<Expr4<B,U,i,k,j,l>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}


// L(i,j,k,l) = R(i,k,l,j)....12
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,i,k,l,j> &rhs)
{
	typedef Aiklj_to_Aijkl<Expr4<A,T,i,k,l,j>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,i,k,l,j> &rhs)
{
	typedef Aiklj_to_Aijkl<Expr4<B,U,i,k,l,j>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,i,k,l,j> &rhs)
{
	typedef Aiklj_to_Aijkl<Expr4<B,U,i,k,l,j>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,i,k,l,j> &rhs)
{
	typedef Aiklj_to_Aijkl<Expr4<B,U,i,k,l,j>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}


// L(i,j,k,l) = R(i,l,j,k)....13
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,i,l,j,k> &rhs)
{
	typedef Ailjk_to_Aijkl<Expr4<A,T,i,l,j,k>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,i,l,j,k> &rhs)
{
	typedef Ailjk_to_Aijkl<Expr4<B,U,i,l,j,k>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,i,l,j,k> &rhs)
{
	typedef Ailjk_to_Aijkl<Expr4<B,U,i,l,j,k>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,i,l,j,k> &rhs)
{
	typedef Ailjk_to_Aijkl<Expr4<B,U,i,l,j,k>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

// L(i,j,k,l) = R(i,l,k,j)....14
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,i,l,k,j> &rhs)
{
	typedef Ailkj_to_Aijkl<Expr4<A,T,i,l,k,j>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,i,l,k,j> &rhs)
{
	typedef Ailkj_to_Aijkl<Expr4<B,U,i,l,k,j>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,i,l,k,j> &rhs)
{
	typedef Ailkj_to_Aijkl<Expr4<B,U,i,l,k,j>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,i,l,k,j> &rhs)
{
	typedef Ailkj_to_Aijkl<Expr4<B,U,i,l,k,j>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

// L(i,j,k,l) = R(j,k,i,l)....15
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,j,k,i,l> &rhs)
{
	typedef Ajkil_to_Aijkl<Expr4<A,T,j,k,i,l>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,j,k,i,l> &rhs)
{
	typedef Ajkil_to_Aijkl<Expr4<B,U,j,k,i,l>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,j,k,i,l> &rhs)
{
	typedef Ajkil_to_Aijkl<Expr4<B,U,j,k,i,l>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,j,k,i,l> &rhs)
{
	typedef Ajkil_to_Aijkl<Expr4<B,U,j,k,i,l>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}



// L(i,j,k,l) = R(j,k,l,i)....16
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,j,k,l,i> &rhs)
{
	typedef Ajkli_to_Aijkl<Expr4<A,T,j,k,l,i>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,j,k,l,i> &rhs)
{
	typedef Ajkli_to_Aijkl<Expr4<B,U,j,k,l,i>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,j,k,l,i> &rhs)
{
	typedef Ajkli_to_Aijkl<Expr4<B,U,j,k,l,i>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,j,k,l,i> &rhs)
{
	typedef Ajkli_to_Aijkl<Expr4<B,U,j,k,l,i>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}



// L(i,j,k,l) = R(j,l,i,k)....17
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,j,l,i,k> &rhs)
{
	typedef Ajlik_to_Aijkl<Expr4<A,T,j,l,i,k>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,j,l,i,k> &rhs)
{
	typedef Ajlik_to_Aijkl<Expr4<B,U,j,l,i,k>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,j,l,i,k> &rhs)
{
	typedef Ajlik_to_Aijkl<Expr4<B,U,j,l,i,k>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,j,l,i,k> &rhs)
{
	typedef Ajlik_to_Aijkl<Expr4<B,U,j,l,i,k>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}


// L(i,j,k,l) = R(j,l,k,i)....18
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,j,l,k,i> &rhs)
{
	typedef Ajlki_to_Aijkl<Expr4<A,T,j,l,k,i>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,j,l,k,i> &rhs)
{
	typedef Ajlki_to_Aijkl<Expr4<B,U,j,l,k,i>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,j,l,k,i> &rhs)
{
	typedef Ajlki_to_Aijkl<Expr4<B,U,j,l,k,i>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,j,l,k,i> &rhs)
{
	typedef Ajlki_to_Aijkl<Expr4<B,U,j,l,k,i>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}


// L(i,j,k,l) = R(k,i,j,l)....19
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,k,i,j,l> &rhs)
{
	typedef Akijl_to_Aijkl<Expr4<A,T,k,i,j,l>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,k,i,j,l> &rhs)
{
	typedef Akijl_to_Aijkl<Expr4<B,U,k,i,j,l>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,k,i,j,l> &rhs)
{
	typedef Akijl_to_Aijkl<Expr4<B,U,k,i,j,l>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,k,i,j,l> &rhs)
{
	typedef Akijl_to_Aijkl<Expr4<B,U,k,i,j,l>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}



// L(i,j,k,l) = R(k,i,l,j)....20
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,k,i,l,j> &rhs)
{
	typedef Akilj_to_Aijkl<Expr4<A,T,k,i,l,j>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,k,i,l,j> &rhs)
{
	typedef Akilj_to_Aijkl<Expr4<B,U,k,i,l,j>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,k,i,l,j> &rhs)
{
	typedef Akilj_to_Aijkl<Expr4<B,U,k,i,l,j>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,k,i,l,j> &rhs)
{
	typedef Akilj_to_Aijkl<Expr4<B,U,k,i,l,j>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}




// L(i,j,k,l) = R(k,j,i,l)....21
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,k,j,i,l> &rhs)
{
	typedef Akjil_to_Aijkl<Expr4<A,T,k,j,i,l>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,k,j,i,l> &rhs)
{
	typedef Akjil_to_Aijkl<Expr4<B,U,k,j,i,l>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,k,j,i,l> &rhs)
{
	typedef Akjil_to_Aijkl<Expr4<B,U,k,j,i,l>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,k,j,i,l> &rhs)
{
	typedef Akjil_to_Aijkl<Expr4<B,U,k,j,i,l>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}



// L(i,j,k,l) = R(k,j,l,i)....22
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,k,j,l,i> &rhs)
{
	typedef Akjli_to_Aijkl<Expr4<A,T,k,j,l,i>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,k,j,l,i> &rhs)
{
	typedef Akjli_to_Aijkl<Expr4<B,U,k,j,l,i>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,k,j,l,i> &rhs)
{
	typedef Akjli_to_Aijkl<Expr4<B,U,k,j,l,i>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,k,j,l,i> &rhs)
{
	typedef Akjli_to_Aijkl<Expr4<B,U,k,j,l,i>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}


// L(i,j,k,l) = R(l,i,j,k)....23
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,l,i,j,k> &rhs)
{
	typedef Alijk_to_Aijkl<Expr4<A,T,l,i,j,k>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,l,i,j,k> &rhs)
{
	typedef Alijk_to_Aijkl<Expr4<B,U,l,i,j,k>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,l,i,j,k> &rhs)
{
	typedef Alijk_to_Aijkl<Expr4<B,U,l,i,j,k>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,l,i,j,k> &rhs)
{
	typedef Alijk_to_Aijkl<Expr4<B,U,l,i,j,k>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}


// L(i,j,k,l) = R(l,i,k,j)....24
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,l,i,k,j> &rhs)
{
	typedef Alikj_to_Aijkl<Expr4<A,T,l,i,k,j>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,l,i,k,j> &rhs)
{
	typedef Alikj_to_Aijkl<Expr4<B,U,l,i,k,j>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,l,i,k,j> &rhs)
{
	typedef Alikj_to_Aijkl<Expr4<B,U,l,i,k,j>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,l,i,k,j> &rhs)
{
	typedef Alikj_to_Aijkl<Expr4<B,U,l,i,k,j>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}



// L(i,j,k,l) = R(l,j,i,k)....24
//
template<class A, class T, char i, char j, char k, char l>	
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<A,T,l,j,i,k> &rhs)
{
	typedef Aljik_to_Aijkl<Expr4<A,T,l,j,i,k>,T> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator=(const Expr4<B,U,l,j,i,k> &rhs)
{
	typedef Aljik_to_Aijkl<Expr4<B,U,l,j,i,k>,U> Expr_Obj;
	Lijkl_equals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator+=(const Expr4<B,U,l,j,i,k> &rhs)
{
	typedef Aljik_to_Aijkl<Expr4<B,U,l,j,i,k>,U> Expr_Obj;
	Lijkl_plusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}

template<class A, class T, char i, char j, char k, char l>	
template<class B, class U>
inline const Expr4<A,T,i,j,k,l> & 
Expr4<A,T,i,j,k,l>::operator-=(const Expr4<B,U,l,j,i,k> &rhs)
{
	typedef Aljik_to_Aijkl<Expr4<B,U,l,j,i,k>,U> Expr_Obj;
	Lijkl_minusequals_Rijkl((*this),Expr_Obj(rhs));
	return *this;
}


#endif
