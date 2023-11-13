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
#ifndef Sum_TensorsRank3_H
#define Sum_TensorsRank3_H


///////////////////////////////////////////////
// Define SUM = A(i,j)+B(i,j)
///////////////////////////////////////////////

template < class A, class B, class T>
class Aijk_plus_Bijk
{     

      private:

	const A  objA;
	const B  objB;

      public:

	~Aijk_plus_Bijk(){ }

	Aijk_plus_Bijk (const A &a, const B &b) : objA(a), objB(b)
	{
		#ifdef USE_ASSERT_Tensor_Operations
			assert (objA.get_dim1() == objB.get_dim1());
			assert (objA.get_dim2() == objB.get_dim2());
			assert (objA.get_dim3() == objB.get_dim3());
		#endif
	}

	T operator () (const int N1,const int N2,const int N3) const
	{
		return objA (N1,N2,N3) + objB (N1,N2,N3);
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

	friend std::ostream & operator<< (std::ostream & os, const Aijk_plus_Bijk & v)
	{
		os << std::endl << "BEGIN:Aijk_plus_Bijk =" << &v << std::endl
		<< ", L = "<< (v.objA) 
		<< ", R = " << (v.objB) 
		<< "END.Aij_plus_Bij " << std::endl;
		return os;
	}
};


#endif
