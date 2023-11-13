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
/*!brief A(i,j,k)*B(j,k) */
//////////////////////////////////////////////////////////

#ifndef Tensor3_contracts_Tensor2_H
#define Tensor3_contracts_Tensor2_H

//////////////////////////////////////////////////////////////////////////
// Contraction of two indices
//////////////////////////////////////////////////////////////////////////

template < class A, class B, class T> 
class Aijk_contracts_Bjk
{
	const A ObjA;
	const B ObjB;
      public:

	Aijk_contracts_Bjk (const A &a, const B &b): ObjA(a), ObjB(b)
	{
		#ifdef USE_ASSERT_Tensor_Operations
		assert (a.get_dim2() == b.get_dim1());
		assert (a.get_dim3() == b.get_dim2());
		#endif    
	}

	T operator () (const int N1) const
	{
		T res=0;
		int n1;
		int n2;
		const int dim1 = ObjA.get_dim2();
		const int dim2 = ObjA.get_dim3();
		for(n1 = 0; n1< dim1;++n1)
		{	
			for(n2 = 0; n2< dim2;++n2)
			{
				res += ObjA(N1,n1,n2)*ObjB(n1,n2);
			}
		}
		return res;
	}

	int get_dim1 () const
	{
		return ObjA.get_dim1();
	}

};

//
// Contraction of 1 index
//

template < class A, class B, class T> 
class Aijk_contracts_Bkl
{
	const A ObjA;
	const B ObjB;
      public:

	Aijk_contracts_Bkl (const A &a, const B &b): ObjA(a), ObjB(b)
	{
		#ifdef USE_ASSERT_Tensor_Operations
		assert (a.get_dim3() == b.get_dim1());
		#endif    
	}

	T operator () (const int N1,const int N2,const int N3) const
	{
		T res=0;
		int n1;
		const int dim1 = ObjA.get_dim3();
		for(n1 = 0; n1< dim1;++n1)
		{	
				res += ObjA(N1,N2,n1)*ObjB(n1,N3);
		}
		return res;
	}

	int get_dim1() const
	{
		return ObjA.get_dim1();
	}
	int get_dim2() const
	{
		return ObjA.get_dim2();
	}
	int get_dim3() const
	{
		return ObjB.get_dim2();
	}

};


#endif
