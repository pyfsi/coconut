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
/*!brief A(i)*B(i) */
//////////////////////////////////////////////////////////

#ifndef Tensor1_contracts_Tensor1_H
#define Tensor1_contracts_Tensor1_H


template<class A, class B, class T >
inline T Ai_contracts_Bi(const A & ExprL, const B & ExprR)
{
	#ifdef USE_ASSERT_Tensor_Operations
		assert (ExprL.get_dim1() == ExprR.get_dim1());
	#endif      
	T res = 0;
	int n1;
	const int dim1 = ExprL.get_dim1();
	for(n1 = 0 ; n1 < dim1; ++n1)
	{
		res += ExprL(n1)*ExprR(n1);
	}
	return res;
}

#endif
