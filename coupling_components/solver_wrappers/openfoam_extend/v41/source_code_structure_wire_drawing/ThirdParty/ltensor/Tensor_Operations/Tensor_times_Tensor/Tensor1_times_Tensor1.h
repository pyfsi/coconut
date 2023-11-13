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
#ifndef Tensor1_times_Tensor1_H
#define Tensor1_times_Tensor1_H


///////////////////////////////////////////////
/*!brief A(i)*B(j) */
///////////////////////////////////////////////

template < class A, class T, class B> class Ai_times_Bj
{
	const A objA;
	const B objB;
      public:
	

	Ai_times_Bj (const A &a, const B & b):
	objA(a), objB(b)
	{
	}

	T operator () (const int N1,const int N2) const
	{
		return objA (N1) * objB (N2);
	}

	int get_dim1() const
	{
		return objA.get_dim1();
	}

	int get_dim2() const
	{
		return objB.get_dim1();
	}
};



#endif
