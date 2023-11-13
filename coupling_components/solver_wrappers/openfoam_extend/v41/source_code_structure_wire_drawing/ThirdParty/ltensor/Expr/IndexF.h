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
#ifndef IndexF_H
#define IndexF_H

#include <iostream>
#include "../Marray/Marray_rank1.h"


class IndexF
{


    public:

    Marray<int,1,TinyArray_base<int,1> > indexes;

	IndexF(){


	}


	IndexF(const IndexF &idx){
		indexes.resize(idx.indexes.get_dim1());
		for(int k=0;k<idx.indexes.get_dim1();k++)
			indexes(k)=idx.indexes(k);
	}


    IndexF(Marray<int,1> idxs)	{
		indexes.resize(idxs.get_dim1());
		for(int k=0;k<idxs.get_dim1();k++)
			indexes(k)=idxs(k);

	}
	IndexF(int* idxs,int n){
		  indexes.resize(n);
        for(int k=0;k<n;k++)
            indexes(k)=idxs[k];

	}

	IndexF(int init,int end,int stride=1){
		rebuild(init,end,stride);
	}


	IndexF(int size){
		 #ifdef CHECK_INDEXF
        assert(size>0);
	    #endif
		 indexes.resize(size);
        for(int k=0;k<size;k++)
            indexes(k)=k;
	}
	~IndexF(){

	}
	template <class A>
	int get_dim1(const A &contenedor) const{
		return indexes.get_dim1();
	}


	int get_dim1() const{
		return indexes.get_dim1();
	}

	template <class A>
	int get_dim2(const A &contenedor) const{
		return indexes.get_dim1();
	}

	template <class A>
	int get_dim3(const A &contenedor) const{
		return indexes.get_dim1();
	}

	template <class A>
	int get_dim4(const A &contenedor) const{
		return indexes.get_dim1();
	}
	int get_size()const{
		 return indexes.get_dim1();
	}

//weird index formulation
    void resize(long dimension){

		indexes.resize(dimension);
	}

	void rebuild(int init,int end, int stride=1){
		 #ifdef CHECK_INDEXF
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


	inline int operator() (const int N) const{
return indexes(N);

	}

	inline int & operator()(const int N){
return indexes(N);
	}

    template <class base2>
	IndexF & operator=(const Marray<int,1,base2 > &idxs){
 indexes.resize(idxs.get_dim1());
		for(int k=0;k<idxs.get_dim1();k++)
			indexes(k)=idxs(k);
	   //make sure sizes are compatible
		//*this).indexes = rterm;
		//el cambio es atribuido a que podemos igualar un Marray de cualquier base, no necesariamente debe ser el mismo tipo
		return *this;

	}

	void showIndexes(){
		std::cout<<indexes;
	}

	IndexF ciclic(const int ncic)
	{
		const int dim1 = get_size();
		IndexF ciclictensor(dim1);
		div_t result;
		for(int i=0;i<dim1;i++)
		{
			result = div(i+ncic+dim1,dim1);
			ciclictensor(i)=(*this)(result.rem);
		}
		return ciclictensor;
	}

    	void sucesion(int init,int stride=1)
	{
		int count=init;
		for(int i=0;i<get_size();i++)
		{
			(*this)(i)=count;
			count = count + stride;
		}
	}


 void set_sucesion(int init,int end,int stride=1)
     {

	#ifdef CHECK_INDEXF
	     assert( (init<end) && (stride>0) )
	#endif
	indexes.resize( (int)ceil ( ( (end-init)) /(double)stride ) );
	int count=0;
	for(int k=init;k<end;k+=stride)
		  indexes(count++)=k;
     }

	int is_number_in(const int n)
	{
		int i=0;
		int idpos = -1;
		while((idpos == -1)&&(i<get_size()))
		{
			if ((*this)(i)==n)
			{ idpos = i;
			};
			++i;
		}
		return idpos;
	}
	bool insert_if_isnot(const int n)
	{
		int idpos = (*this).is_number_in(n);
		std::cout << "idpos " << idpos;
		bool inserted = false;
		IndexF iFaux(*this);
		if (idpos==-1)
		{
			inserted = true;
			iFaux = (*this);
			(*this).indexes.resize(get_size()+1);
			for(int i=0;i<get_size()-1;i++)
			{
				(*this)(i)=iFaux(i);
			}
			(*this)(get_size()-1)=n;
		}
		return inserted;
	}
	bool insert_value_when_value(const int ivalue, const int atvalue)
	{
		int idpos = (*this).is_number_in(atvalue);
		bool inserted = false;
		if (idpos>=0)
		{
			inserted = true;
			(*this)(idpos)=ivalue;
		}
		return inserted;
	}


	template <class U>
    inline IndexF & operator+= (const U &u){
        for(int i=0;i<get_dim1();i++)
            (*this)(i)+=u;
        return *this;

    }


    template <class U>
    inline IndexF & operator-= (const U &u){
        for(int i=0;i<get_dim1();i++)
            (*this)(i)-=u;
        return *this;

    }



	friend std::ostream & operator<< (std::ostream & os, const IndexF & p)
	{
	std::cout << "IndexF = " << p.indexes;
	return os;
	}

};


#endif

