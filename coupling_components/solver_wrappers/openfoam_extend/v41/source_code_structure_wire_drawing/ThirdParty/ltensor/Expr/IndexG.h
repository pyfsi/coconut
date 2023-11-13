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
#ifndef IndexG_H
#define IndexG_H

//specialization for indexG

template < char i >
class Index <i,1>
{
	

public:
	int size;
	inline int operator()(int n){
		return n;
	}

	inline int operator()(int n)const {
		return n;
	}

	template <class A>
	inline int get_dim1(const A &contenedor)const{
		return contenedor.get_dim1();
	}

	template <class A>
	inline int get_dim2(const A &contenedor)const{
		return contenedor.get_dim2();
	}

	template <class A>
	inline int get_dim3(const A &contenedor)const{
		return contenedor.get_dim3();
	}

	template <class A>
	inline int get_dim4(const A &contenedor)const{
		return contenedor.get_dim4();
	}

	

};






#endif
