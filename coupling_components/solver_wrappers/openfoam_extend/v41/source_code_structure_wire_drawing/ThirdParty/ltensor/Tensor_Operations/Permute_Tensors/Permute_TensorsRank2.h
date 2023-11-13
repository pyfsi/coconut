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
/////////////////////////////////////////////////////
// Given Expr2=A(j,i) define a new permuted Expr2 P object
// so to define: B(i,j)=A(j,i)
/////////////////////////////////////////////////////

#ifndef Permute_TensorsRank2_H
#define Permute_TensorsRank2_H

////////////////////////////////////////////
// Define B(i,j)=A(i,j)  ==> No permutation
///////////////////////////////////////////

template < class A, class T> 
class Aij_to_Aij
{
	A  TA;
      public:

	Aij_to_Aij (const A &a): TA(a)
	{
	}

	int get_dim1 () const
	{
		return TA.get_dim1();
	}
	
	int get_dim2 () const
	{
		return TA.get_dim2();
	}
	
	T operator () (const int N1,const int N2) const
	{
		return TA(N1,N2);
	}
};

////////////////////////////////////////////
// Define B(i,j)=A(j,i)  
///////////////////////////////////////////


template < class A, class T> 
class Aji_to_Aij
{
	A TA;
      public:

	Aji_to_Aij (const A &a): TA(a)
	{
	}

	int get_dim1 () const
	{
		return TA.get_dim2();
	}
	
	int get_dim2 () const
	{
		return TA.get_dim1();
	}
	
	T operator () (const int N1,const int N2) const
	{
		return TA(N2,N1);
	}
};

#endif
