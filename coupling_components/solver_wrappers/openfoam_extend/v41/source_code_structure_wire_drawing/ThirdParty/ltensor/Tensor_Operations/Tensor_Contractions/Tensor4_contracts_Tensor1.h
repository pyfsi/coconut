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
// A(i,j,k,l)*B(l), etc
//////////////////////////////////////////////////////////

#ifndef Tensor4_contracts_Tensor1_H
#define Tensor4_contracts_Tensor1_H


template < class A, class B, class T> 
class Aijkl_contracts_Bl
{
	const A ObjA;
	const B ObjB;
      public:

	Aijkl_contracts_Bl (const A &a, const B &b): ObjA(a), ObjB(b)
	{
		#ifdef USE_ASSERT_Tensor_Operations
		assert (a.get_dim4() == b.get_dim1());
		#endif    
	}

	T operator () (const int N1, const int N2, const int N3) const
	{
		T res=0;
		int n1;
		const int dim1 = ObjA.get_dim4();
		for(n1 = 0; n1< dim1;++n1)
		{	
			res += ObjA(N1,N2,N3,n1)*ObjB(n1);
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
		return ObjA.get_dim3();
	}

};

#endif
