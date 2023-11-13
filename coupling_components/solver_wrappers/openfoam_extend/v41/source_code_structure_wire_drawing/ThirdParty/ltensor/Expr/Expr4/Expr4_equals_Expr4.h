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
// Expr2 = Expr2
//////////////////////////////////////////////////////////
#ifndef Expr2_equals_Expr2_H
#define Expr2_equals_Expr2_H


template<class A, class B, class T, class U, char i, char j>
inline void Expr2_equals_Expr2(Expr2<A,T,i,j> & TLeft, const Expr2<B,U,i,j> TRight)
{
	#ifdef USE_ASSERT_Expr2
		assert ( TLeft.get_dim1() == TRight.get_dim1() );
		assert ( TLeft.get_dim2() == TRight.get_dim2() );
	#endif
	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	int n1;
	int n2;
 	for(n1 = 0; n1<dim1;++n1)
	{ 
		for(n2 = 0; n2<dim2;++n2)
		{
			TLeft(n1,n2) = TRight(n1,n2);
		}
	}
}

template<class A, class B, class T, class U, char i, char j>
inline void Expr2_plusequals_Expr2(Expr2<A,T,i,j> & TLeft, const Expr2<B,U,i,j> TRight)
{
	#ifdef USE_ASSERT_Expr2
		assert ( TLeft.get_dim1() == TRight.get_dim1() );
		assert ( TLeft.get_dim2() == TRight.get_dim2() );
	#endif
	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	int n1;
	int n2;
 	for(n1 = 0; n1<dim1;++n1)
	{ 
		for(n2 = 0; n2<dim2;++n2)
		{
			TLeft(n1,n2) += TRight(n1,n2);
		}
	}
}

template<class A, class B, class T, class U, char i, char j>
inline void Expr2_minusequals_Expr2(Expr2<A,T,i,j> & TLeft, const Expr2<B,U,i,j> TRight)
{
	#ifdef USE_ASSERT_Expr2
		assert ( TLeft.get_dim1() == TRight.get_dim1() );
		assert ( TLeft.get_dim2() == TRight.get_dim2() );
	#endif
	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	int n1;
	int n2;
 	for(n1 = 0; n1<dim1;++n1)
	{ 
		for(n2 = 0; n2<dim2;++n2)
		{
			TLeft(n1,n2) -= TRight(n1,n2);
		}
	}
}


#endif
