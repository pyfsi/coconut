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
#ifndef Plus_TensorsRank1_H
#define Plus_TensorsRank1_H


///////////////////////////////////////////////
// Define PLUS = +A(i)
///////////////////////////////////////////////

template < class A, class T>
class plus_Ai
{     

      private:

	const A  objA;

      public:

	~plus_Ai(){ }

	plus_Ai (const A &a) : objA(a){ }

	T operator () (const int N1) const
	{
		return objA (N1);
	}

	int get_dim1() const
	{
		return objA.get_dim1();
	}

	friend std::ostream & operator<< (std::ostream & os, const plus_Ai & v)
	{
		os << std::endl << "BEGIN:plus_Ai =" << &v << std::endl
		<< ", A = "<< (v.objA) 
		<< "END.plus_Ai " << std::endl;
		return os;
	}
};

#endif
