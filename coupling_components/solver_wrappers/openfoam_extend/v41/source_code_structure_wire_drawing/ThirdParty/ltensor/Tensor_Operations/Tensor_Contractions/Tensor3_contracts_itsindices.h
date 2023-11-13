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
/////////////////////////////////////////////////////////
// Define FUNCTIONS to perform operations of type:
// Tensor3 = A(i,j,k)   ----> Tensor1 = A(i,j,j)
//////////////////////////////////////////////////////////

#ifndef Tensor3_contracts_itsindices_H
#define Tensor3_contracts_itsindices_H


template < class A, class T> 
class Aijj_contraction
{
	const A ObjA;
      public:

	Aijj_contraction (const A &a): ObjA(a)
	{
		#ifdef USE_ASSERT_Tensor_Operations
		assert (a.get_dim2() == a.get_dim3());
		#endif    
	}

	T operator () (const int N1) const
	{
		T res=0;
		int n2;
		const int dim2 = ObjA.get_dim2();
		for(n2 = 0; n2< dim2;++n2)
		{	
			res += ObjA(N1,n2,n2);
		}
		return res;
	}

	int get_dim1 () const
	{
		return ObjA.get_dim1();
	}

};

template < class A, class T> 
class Ajij_contraction
{
	const A ObjA;
      public:

	Ajij_contraction (const A &a): ObjA(a)
	{
		#ifdef USE_ASSERT_Tensor_Operations
		assert (a.get_dim1() == a.get_dim3());
		#endif    
	}

	T operator () (const int N1) const
	{
		T res=0;
		int n2;
		const int dim2 = ObjA.get_dim2();
		for(n2 = 0; n2< dim2;++n2)
		{	
			res += ObjA(n2,N1,n2);
		}
		return res;
	}

	int get_dim1 () const
	{
		return ObjA.get_dim2();
	}

};

template < class A, class T> 
class Ajji_contraction
{
	const A ObjA;
      public:

	Ajji_contraction (const A &a): ObjA(a)
	{
		#ifdef USE_ASSERT_Tensor_Operations
		assert (a.get_dim1() == a.get_dim3());
		#endif    
	}

	T operator () (const int N1) const
	{
		T res=0;
		int n2;
		const int dim2 = ObjA.get_dim2();
		for(n2 = 0; n2< dim2;++n2)
		{	
			res += ObjA(n2,n2,N1);
		}
		return res;
	}

	int get_dim1 () const
	{
		return ObjA.get_dim3();
	}

};


#endif
