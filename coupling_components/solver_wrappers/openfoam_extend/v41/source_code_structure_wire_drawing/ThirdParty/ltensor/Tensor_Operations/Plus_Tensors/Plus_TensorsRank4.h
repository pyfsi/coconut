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
#ifndef Plus_TensorsRank4_H
#define Plus_TensorsRank4_H


///////////////////////////////////////////////
// Define Plus = +A(i,j,k,l)
///////////////////////////////////////////////

template < class A, class T>
class plus_Aijkl
{     

      private:

	const A  objA;

      public:

	~plus_Aijkl(){ }

	plus_Aijkl (const A &a) : objA(a)	{ }

	T operator () (const int N1,const int N2,const int N3,const int N4) const
	{
		return objA (N1,N2,N3,N4);
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
 
	friend std::ostream & operator<< (std::ostream & os, const plus_Aijkl & v)
	{
		os << std::endl << "BEGIN:plus_Aijkl =" << &v << std::endl
		<< ", A = "<< (v.objA)  
		<< "END.plus_Aijkl " << std::endl;
		return os;
	}
};


#endif
