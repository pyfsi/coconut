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
#ifndef TensorRank4_times_Scalar_H
#define TensorRank4_times_Scalar_H


///////////////////////////////////////////////
// Define multiplication by scalar = A(i,j,k,l)*u
///////////////////////////////////////////////

template < class A, class T, class U> class Aijkl_times_u
{
	const A objA;
	const U u;
      public:
	

	Aijkl_times_u (const A &a, const U & u0): objA(a), u(u0) {}

	T operator () (const int N1,const int N2,const int N3,const int N4) const
	{
		return objA (N1,N2,N3,N4) * u;
	}

	int get_dim1() const
	{
		return objA.get_dim1();
	}

	int get_dim2() const
	{
		return objA.get_dim2();
	}
	int get_dim3() const
	{
		return objA.get_dim3();
	}
	int get_dim4() const
	{
		return objA.get_dim4();
	}
};



#endif
