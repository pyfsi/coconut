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
#ifndef IndexH_H
#define IndexH_H

#include <iostream>
#include "../Marray/Marray_rank1.h"

template <char i>
class Index<i,2>
{


    public:

    Marray<int,1,TinyArray_base<int,1> > indexes;

    Index();


	Index(const Index &idx);


    Index(Marray<int,1> idxs);
    Index(int* idxs,int n);
    Index(int init,int end,int stride=1);

	 Index(int size);
    ~Index();
    int get_dim1() const;
    int get_size()const;

//weird index formulation

	void rebuild(int init,int end, int stride=1);
	inline int operator() (const int N) const;

	inline int & operator()(const int N);

    template <class base2>
	Index<i,2> & operator=(const Marray<int,1,base2 > &idxs);

	void showIndexes();


	friend std::ostream & operator<< (std::ostream & os, const Index<i,2> & p)
	{
	std::cout << "IndexH = " << p.indexes;
	return os;
	}
	template <class A>
	inline int get_dim1(const A &contenedor) const;
	template <class A>
	inline int get_dim2(const A &contenedor) const;
	template <class A>
	inline int get_dim3(const A &contenedor) const;
	template <class A>
	inline int get_dim4(const A &contenedor) const;


};


#include "IndexH_bis.h"

#endif

