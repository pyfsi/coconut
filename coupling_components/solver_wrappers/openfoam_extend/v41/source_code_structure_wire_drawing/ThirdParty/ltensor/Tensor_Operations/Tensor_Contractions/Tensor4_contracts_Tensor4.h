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
// A(i,j,k,l)*B(i,j,k,l), etc
//////////////////////////////////////////////////////////

#ifndef Tensor4_contracts_Tensor4_H
#define Tensor4_contracts_Tensor4_H

////////////////////////////////////////////////////////////////////
// FULL CONTRACTIONS OF INDiCES 
////////////////////////////////////////////////////////////////////

//
//A(i,j,k,l)*B(i,j,k,l)
//
template<class A, class B, class T >
inline T Aijkl_contracts_Bijkl(const A & ExprL, const B & ExprR)
{
	#ifdef USE_ASSERT_Tensor_Operations
		assert (ExprL.get_dim1() == ExprR.get_dim1());
		assert (ExprL.get_dim2() == ExprR.get_dim2());
		assert (ExprL.get_dim3() == ExprR.get_dim3());
		assert (ExprL.get_dim4() == ExprR.get_dim4());
	#endif      
	T res = 0;
	int n1;
	int n2;
	int n3;
	int n4;
	const int dim1 = ExprL.get_dim1();
	const int dim2 = ExprL.get_dim2();
	const int dim3 = ExprL.get_dim3();
	const int dim4 = ExprL.get_dim4();
	for(n1 = 0 ; n1 < dim1; ++n1)
	{
		for(n2 = 0 ; n2 < dim2; ++n2)
		{
			for(n3 = 0 ; n3 < dim3; ++n3)
			{
				for(n4 = 0 ; n4 < dim4; ++n4)
				{
					res += ExprL(n1,n2,n3,n4)*ExprR(n1,n2,n3,n4);
				}
			}
		}
	}
	return res;
}



////////////////////////////////////////////////////////////////////
// TRIPLE CONTRACTION OF INDiCES 
////////////////////////////////////////////////////////////////////

//
//A(i,j,k,l)*B(j,k,l,m)----> C(i,m)
//

template < class A, class B, class T> 
class Aijkl_contracts_Bjklm
{
	const A ObjA;
	const B ObjB;
      public:

	Aijkl_contracts_Bjklm (const A &a, const B &b): ObjA(a), ObjB(b)
	{
		#ifdef USE_ASSERT_Tensor_Operations
		assert (ObjA.get_dim2() == ObjB.get_dim1());
		assert (ObjA.get_dim3() == ObjB.get_dim2());
		assert (ObjA.get_dim4() == ObjB.get_dim3());
		#endif    
	}

	T operator () (const int N1,const int N2) const
	{
		T res=0;
		int n1;
		int n2;
		int n3;
		const int dim1 = ObjA.get_dim2();
		const int dim2 = ObjA.get_dim3();
		const int dim3 = ObjA.get_dim4();
		for(n1 = 0; n1< dim1;++n1)
		{	
			for(n2 = 0; n2< dim2;++n2)
			{	
				for(n3 = 0; n3< dim3;++n3)
				{	
					res += ObjA(N1,n1,n2,n3)*ObjB(n1,n2,n3,N2);
				}
			}
		}
		return res;
	}

	int get_dim1() const
	{
		return ObjA.get_dim1();
	}

	int get_dim2() const
	{
		return ObjB.get_dim4();
	}

};

////////////////////////////////////////////////////////////////////
// DOUBLE CONTRACTION OF INDiCES 
////////////////////////////////////////////////////////////////////

//
//A(i,j,k,l)*B(k,l,m,n)----> C(i,j,m,n)
//

template < class A, class B, class T> 
class Aijkl_contracts_Bklmn
{
	const A ObjA;
	const B ObjB;
      public:

	Aijkl_contracts_Bklmn (const A &a, const B &b): ObjA(a), ObjB(b)
	{
		#ifdef USE_ASSERT_Tensor_Operations
		assert (ObjA.get_dim3() == ObjB.get_dim1());
		assert (ObjA.get_dim4() == ObjB.get_dim2());
		#endif    
	}

	T operator () (const int N1,const int N2,const int N3,const int N4) const
	{
		T res=0;
		int n1;
		int n2;
		const int dim1 = ObjA.get_dim3();
		const int dim2 = ObjA.get_dim4();
		for(n1 = 0; n1< dim1;++n1)
		{	
			for(n2 = 0; n2< dim2;++n2)
			{	
				res += ObjA(N1,N2,n1,n2)*ObjB(n1,n2,N3,N4);
			}
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
		return ObjB.get_dim3();
	}

	int get_dim4() const
	{
		return ObjB.get_dim4();
	}
};

#endif
