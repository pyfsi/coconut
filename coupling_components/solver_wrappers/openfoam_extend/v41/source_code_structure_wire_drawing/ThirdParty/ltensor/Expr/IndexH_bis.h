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
#ifndef Index_BIS_H
#define Index_BIS_H

#include <iostream>
#include "../Marray/Marray_rank1.h"
#include "Index.h"

	template <char i>
    Index<i,2>::Index()
    {
    }

template <char i>
    Index<i,2>::Index(const Index<i,2> &idx){
        indexes.resize(idx.indexes.get_dim1());
		for(int k=0;k<idx.indexes.get_dim1();k++)
			indexes(k)=idx.indexes(k);


        }


template <char i>
    Index<i,2>::Index(Marray<int,1> idxs){
		indexes.resize(idxs.get_dim1());
		for(int k=0;k<idxs.get_dim1();k++)
			indexes(k)=idxs(k);
     }
	template <char i>
    Index<i,2>::Index(int* idxs,int n){

        indexes.resize(n);
        for(int k=0;k<n;k++)
            indexes(k)=idxs[k];

     }
	template <char i>
    Index<i,2>::Index(int init,int end,int stride){
		rebuild(init,end,stride);


    }
template <char i>
	 Index<i,2>::Index(int size){
	    #ifdef CHECK_INDEXH
        assert(size>0);
	    #endif
		 indexes.resize(size);
        for(int k=0;k<size;k++)
            indexes(k)=k;


	}
	 template <char i>
    Index<i,2>::~Index(){

    }
template <char i>
    int Index<i,2>::get_dim1() const
    {
        return indexes.get_dim1();
    }
	template <char i>
    int Index<i,2>::get_size()const {
        return indexes.get_dim1();
    }

//weird index formulation
	//Ojo q end es n+1
template <char i>
	void Index<i,2>::rebuild(int init,int end, int stride){
		 #ifdef CHECK_INDEXH
        assert(stride!=0)
        #endif


        int count=0;
        if(stride>0){
            indexes.resize( ( (int)ceil ( ( (end-init)+1) /(double)stride )  ) -1 );
        for(int k=init;k<end;k+=stride)

            indexes(count++)=k;
        }else{
            indexes.resize( ((int)ceil ( ( (init-end)+1) /(double)(-stride) ))-1 );

            for(int k=init;k>end;k+=stride)
            indexes(count++)=k;
        }
	}
	template <char i>
	inline int Index<i,2>::operator() (const int N) const
	{
		return indexes(N);
	}

	template <char i>
	inline int & Index<i,2>::operator()(const int N)
	{
		return indexes(N);
	}

	template <char i>
    template <class base2>
	Index<i,2> & Index<i,2>::operator=(const Marray<int,1,base2 > &idxs)
	{

	    indexes.resize(idxs.get_dim1());
		for(int k=0;k<idxs.get_dim1();k++)
			indexes(k)=idxs(k);
	   //make sure sizes are compatible
		//*this).indexes = rterm;
		//el cambio es atribuido a que podemos igualar un Marray de cualquier base, no necesariamente debe ser el mismo tipo
		return *this;
	}

template <char i>
    void Index<i,2>::showIndexes(){
        std::cout<<indexes;
    }

template < char i>
template < class A>
	inline int Index<i,2>::get_dim1(const A &contenedor) const {
		return indexes.get_dim1();
	}

	template < char i>
template < class A>
	inline int Index<i,2>::get_dim2(const A &contenedor) const {
		return indexes.get_dim1();
	}

	template < char i>
template < class A>
	inline int Index<i,2>::get_dim3(const A &contenedor) const {
		return indexes.get_dim1();
	}

	template < char i>
template < class A>
	inline int Index<i,2>::get_dim4(const A &contenedor) const {
		return indexes.get_dim1();
	}

template <char i>
std::ostream & operator<< (std::ostream & os, const Index<i,2> & p)
{
	std::cout << "Index = " << p.indexes;
	return os;
}
#endif

