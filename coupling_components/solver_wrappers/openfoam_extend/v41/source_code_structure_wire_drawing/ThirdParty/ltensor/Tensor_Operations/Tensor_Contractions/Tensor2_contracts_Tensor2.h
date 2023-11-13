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
/*!brief A(i,j)*B(j) */
//////////////////////////////////////////////////////////

#ifndef Tensor2_contracts_Tensor2_H
#define Tensor2_contracts_Tensor2_H

////////////////////////////////////////////////////////////////////
// FULL CONTRACTIONS OF INDiCES 
////////////////////////////////////////////////////////////////////

//
/*!brief A(i,j)*B(i,j) */
//
template<class A, class B, class T >
inline T Aij_contracts_Bij(const A & ExprL, const B & ExprR)
{
	#ifdef USE_ASSERT_Tensor_Operations
		assert (ExprL.get_dim1() == ExprR.get_dim1());
		assert (ExprL.get_dim2() == ExprR.get_dim2());
	#endif      
	T res = 0;
	int n1;
	int n2;
	const int dim1 = ExprL.get_dim1();
	const int dim2 = ExprL.get_dim2();
	for(n1 = 0 ; n1 < dim1; ++n1)
	{
		for(n2 = 0 ; n2 < dim2; ++n2)
		{
			res += ExprL(n1,n2)*ExprR(n1,n2);
		}
	}
	return res;
}

//
/*!brief A(i,j)*B(j,i) */
//
template<class A, class B, class T >
inline T Aij_contracts_Bji(const A & ExprL, const B & ExprR)
{
	#ifdef USE_ASSERT_Tensor_Operations
		assert (ExprL.get_dim1() == ExprR.get_dim2());
		assert (ExprL.get_dim2() == ExprR.get_dim1());
	#endif      
	T res = 0;
	int n1;
	int n2;
	const int dim1 = ExprL.get_dim1();
	const int dim2 = ExprL.get_dim2();
	for(n1 = 0 ; n1 < dim1; ++n1)
	{
		for(n2 = 0 ; n2 < dim2; ++n2)
		{
			res += ExprL(n1,n2)*ExprR(n2,n1);
		}
	}
	return res;
}

////////////////////////////////////////////////////////////////////
// PARTIAL CONTRACTIONS OF INDiCES 
////////////////////////////////////////////////////////////////////

//
/*!brief A(i,k)*B(k,j) */
//

template < class A, class B, class T> 
class Aik_contracts_Bkj
{
	const A ObjA;
	const B ObjB;
      public:

	Aik_contracts_Bkj (const A &a, const B &b):
	ObjA(a), ObjB(b)
	{
		#ifdef USE_ASSERT_Tensor_Operations
		assert (a.get_dim2() == b.get_dim1());
		#endif    
	}

	T operator () (const int N1,const int N2) const
	{
		T res=0;
		int n1;
		const int dim1 = ObjA.get_dim2();
		for(n1 = 0; n1< dim1;++n1)
		{	
			res += ObjA(N1,n1)*ObjB(n1,N2);
		}
		return res;
	}

	int get_dim1() const
	{
		return ObjA.get_dim1();
	}

	int get_dim2() const
	{
		return ObjB.get_dim2();
	}

};

//
/*!brief A(i,k)*B(j,k) */
//

template < class A, class B, class T> 
class Aik_contracts_Bjk
{
	const A ObjA;
	const B ObjB;
      public:

	Aik_contracts_Bjk (const A &a, const B &b):
	ObjA(a), ObjB(b)
	{
		#ifdef USE_ASSERT_Tensor_Operations
		assert (ObjA.get_dim2() == ObjB.get_dim2());
		#endif    
	}

	T operator () (const int N1,const int N2) const
	{
		T res=0;
		int n1;
		const int dim1 = ObjA.get_dim2();
		for(n1 = 0; n1< dim1;++n1)
		{	
			res += ObjA(N1,n1)*ObjB(N2,n1);
		}
		return res;
	}


	int get_dim1() const
	{
		return ObjA.get_dim1();
	}

	int get_dim2() const
	{
		return ObjB.get_dim1();
	}

};

//
/*!brief A(k,i)*B(j,k) = B(j,k)*A(k,i) */
//

template < class A, class B, class T> 
class Aki_contracts_Bjk
{
	const A ObjA;
	const B ObjB;
      public:

	Aki_contracts_Bjk (const A &a, const B &b):
	ObjA(a), ObjB(b)
	{
		#ifdef USE_ASSERT_Tensor_Operations
		assert (ObjA.get_dim1() == ObjB.get_dim2());
		#endif    
	}

	T operator () (const int N1,const int N2) const
	{
		T res=0;
		int n1;
		const int dim1 = ObjA.get_dim1();
		for(n1 = 0; n1< dim1;++n1)
		{	
			res += ObjA(n1,N1)*ObjB(N2,n1);
		}
		return res;
	}

	int get_dim1() const
	{
		return ObjA.get_dim2();
	}

	int get_dim2() const
	{
		return ObjB.get_dim1();
	}

};


//
/*!brief A(k,i)*B(k,j) */
//

template < class A, class B, class T> 
class Aki_contracts_Bkj
{
	const A ObjA;
	const B ObjB;
      public:

	Aki_contracts_Bkj (const A &a, const B &b):
	ObjA(a), ObjB(b)
	{
		#ifdef USE_ASSERT_Tensor_Operations
		assert (ObjA.get_dim1() == ObjB.get_dim1());
		#endif    
	}

	T operator () (const int N1,const int N2) const
	{
		T res=0;
		int n1;
		const int dim1 = ObjA.get_dim1();
		for(n1 = 0; n1< dim1;++n1)
		{	
			res += ObjA(n1,N1)*ObjB(n1,N2);
		}
		return res;
	}

	int get_dim1() const
	{
		return ObjA.get_dim2();
	}

	int get_dim2() const
	{
		return ObjB.get_dim2();
	}

};



#endif
