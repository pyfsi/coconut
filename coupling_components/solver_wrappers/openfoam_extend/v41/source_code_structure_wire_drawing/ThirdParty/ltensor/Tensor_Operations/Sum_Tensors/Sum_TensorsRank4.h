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
#ifndef Sum_TensorsRank4_H
#define Sum_TensorsRank4_H


///////////////////////////////////////////////
// Define SUM = A(i,j,k,l)+B(i,j,k,l)
///////////////////////////////////////////////

template < class A, class B, class T>
class Aijkl_plus_Bijkl
{     

      private:

	const A  objA;
	const B  objB;

      public:

	~Aijkl_plus_Bijkl(){ }

	Aijkl_plus_Bijkl (const A &a, const B &b) : objA(a), objB(b)
	{
		#ifdef USE_ASSERT_Tensor_Operations
			assert (objA.get_dim1() == objB.get_dim1());
			assert (objA.get_dim2() == objB.get_dim2());
			assert (objA.get_dim3() == objB.get_dim3());
			assert (objA.get_dim4() == objB.get_dim4());
		#endif
	}

	T operator () (const int N1,const int N2,const int N3,const int N4) const
	{
		return objA (N1,N2,N3,N4) + objB (N1,N2,N3,N4);
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

	friend std::ostream & operator<< (std::ostream & os, const Aijkl_plus_Bijkl & v)
	{
		os << std::endl << "BEGIN:Aijkl_plus_Bijkl =" << &v << std::endl
		<< ", L = "<< (v.objA) 
		<< ", R = " << (v.objB) 
		<< "END.Aijkl_plus_Bijkl " << std::endl;
		return os;
	}
};


#endif
