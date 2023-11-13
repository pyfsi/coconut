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
// -*- c++ -*-
#ifndef Marray_rank1_H
#define Marray_rank1_H

#define CHECK				\
	assert(n<get_dim1());	\
	assert(n>=0);

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <vector>

#include "../base/Array_base.h"
#include "../meta/Metaprograms.h"
#include  "../Expr/Expr1/Expr1.h"
#include  "../Expr/Expr1/F_Expr1.h"
#include "../Tensor_Operations/Encapsulate_Tensors/Encapsulate_Tensors.h"

template < class A, class T, char i > class Expr1;

#define rank 1

//
// FUNTIONS DECLARATIONS in order to FRIEND Function "ostream <<"  to work
//
template < class T,class base> class Marray <T,rank,base>;

template<class T,class base>
std::ostream & operator << (std::ostream & os, const Marray<T,rank,base> & v);

/*! \brief Class Representing a Tensor of rank 1
	*
	*	This class represents a Tensor of rank 1. It has no intrinsic funtionality instead it inherits from the template
	*	parameter base. The main purpose of this class is to implement all the tensor expression's functionality.
	*	This way, any class Base can have tensor expression compatibility without any code rewritting. This class is nothing
	*	but a wrapper that contains tensor expression comptatibility and some generic functions.
	*	\tparam T The type of the data the Marray will hold.
	*	\tparam rank The rank of the Tensor to instantiate. This parameter is used to choose between the specializations of this class.
	*	\tparam base The base class that provides all the functionality of a Tensor of Rank rank.
	*/
template < class T , class base >
class Marray <T,rank,base> : public base
{

    public:
    ///////////////
    //Constructors
    ///////////////

    Marray (long dimension):base(dimension){
	    (*this)=0;

        }

     Marray (long dimension, T valor):base(dimension){
        *this=valor;
        }

    Marray(){

        };


	//ordering makes no sense here
    Marray(long* dimensions):base(dimensions){
	    (*this)=0;
        }


    //copy operator>
    Marray(const Marray<T, rank,base> &R){

            *this=R;

        }

    // norm =  1.. N < Inf

    T getNormN(int norm) const {

	#ifdef CHECK_OOB
	    assert(norm!=0);
	#endif
        double ret=pow((double)absTrait<T>::function( (*this)(0)) ,norm);
        for(int i=1;i<get_dim1();i++)
            ret+=pow((double)absTrait<T>::function( (*this)(i)),norm);
        return pow(ret,1.0/norm);

        }

    T getNormZero() const {
        T ret=0;
		for(int i=0;i<get_dim1();i++)
            ret+= ((*this)(i)==0?0:1);
        return ret;
    }

    T getNormInf() const {
        T ret=absTrait<T>::function((*this)(0));
		for(int i=1;i<get_dim1();i++)
            if( absTrait<T>::function( (*this)(i)) >ret) ret = absTrait<T>::function((*this)(i));
        return ret;
        }


    T getMax() const {
         T ret=(*this)(0);
		for(int i=1;i<get_dim1();i++)
            if( (*this)(i) >ret) ret = (*this)(i);
        return ret;
    }

    T getMin() const {
        T ret=(*this)(0);
        for (int i=0;i<get_dim1();i++)
            if((*this)(i)<ret) ret=(*this)(i);
	return ret;
        }


	void normalize(){

		double norm = getNormN(2);
		if(norm!=0)
			for (int i=0;i<get_dim1();i++)
				(*this)(i)/=(T)norm;

	}

    void quicksort(int left,int right){

        if (right>left){
            int pivotIndex=left;
            int pivotNewIndex= partition(left,right,pivotIndex);
            quicksort(left,pivotNewIndex-1);
            quicksort(pivotNewIndex+1,right);
        }

     }

private:
        int partition(int left,int right, int pivotIndex){
             T pivotValue=(*this)(pivotIndex);
             swapm<int>::doit( (*this)(pivotIndex),(*this)(right));
             int storeIndex=left;
             for(int i=left;i<right;i++){
                 if((*this)(i)<=pivotValue){
                     swapm<int>::doit((*this)(i),(*this)(storeIndex));
                     storeIndex+=1;
                     }
                 }
             swapm<int>::doit((*this)(storeIndex),(*this)(right));
            return storeIndex;

         }



public:
	   void inline operator=(const char *data){

			char delim;
		   std::stringstream stream;
		   stream<<data;

		   for(int i=0;i<get_dim1();i++){

			   stream>>(*this)(i);
			}

				stream>>delim;
				assert(delim==';');





	   }

	template <class fType,class b>
	Marray<T,rank,base>& operator=(const Marray<fType,rank,b> &rterm){
	     	base::operator=(*dynamic_cast<const b*>(&rterm));
       return *this;
    }



	Marray<T,rank,base>& operator=(const Marray<T,rank,base> &rterm){

		base::operator=(*dynamic_cast<const base*>(&rterm));
		return *this;

    }


template <class A, class U>
	Marray<T,rank,base>& operator=(const F_Expr1 < A, U> &expr){

		for(int i=0;i<expr.get_dim1();i++){
			(*this)(i)=expr(i);
		}
		return (*this);
}

Marray<T,2> reshape(const int n1, const int n2)
	{
#ifdef CHECK_OOB
		assert(n1*n2==get_dim1());
#endif
		Marray<T,2> b(n1,n2);
		for(int i=0;i<n1;i++)
		{
			for(int j=0;j<n2;j++)
			{ b(i,j)=(*this)(i*n2+j);
			}
		}
		return b;
    }


    inline const T operator() (const int n) const{
		#ifdef CHECK_OOB
			CHECK
		#endif
        return base::operator()(n);
        }

	inline T & operator() (const int n){
		#ifdef CHECK_OOB
			CHECK
		#endif

        return base::operator()(n);
	    }

	void resize(unsigned long dimension){

		base::resize(&dimension);
		(*this)=0;
	}

    inline int get_size() const{
        return base::size[0];
        }

	inline int get_dim1() const{
	    return base::size[0];
	    }

	inline void show_shape(){
	    std::cout<<get_dim1()<<std::endl;
	    }

    	void sucesion(int init,int stride=1)
	{
		int count=init;
		for(int i=0;i<get_dim1();i++)
		{
			(*this)(i)=count;
			count = count + stride;
		}
	}

	template <class fType,class b>
	void repVec(const Marray<fType,rank,b> &rterm, int times){
		int sizer=rterm.get_dim1();
		unsigned int total=sizer*times;
		resize(total);
		for(unsigned int i=0;i<total;++i)
			(*this)(i)=rterm(i%sizer);
	}

	void saveToSpacedFile(const char* fileName){

	    std::ofstream out(fileName,std::ios::out|std::ios::trunc);

        #ifdef CHECK_FILES
        assert(out.is_open());
        #endif

        out<<get_dim1()<<std::endl;
        for(int i=0;i<get_dim1();i++)
            out<<(*this)(i)<<" ";

        out.close();



	}

    void loadFromSpacedFile(const char* fileName){

		std::ifstream in(fileName);

        #ifdef CHECK_FILES
        assert(in.is_open());
        #endif

		int size;
        in>>size;

        #ifdef CHECK_FILES
        assert(in.is_open());
        assert(size>0);

        #endif

		resize(size);
        for(int i=0;i<size;i++)
			in>>(*this)(i);

        in.close();

        }




typedef  Marray<T,rank,base> CSELFTYPE;


friend std::ostream & operator<< <T,base>(std::ostream & os, const Marray<T,rank,base> & v);


//F_EXPRS

//last template parameter needed due to language restrictions


F_Expr1 < Encapsulate_to_Expr1<Marray<T,rank,base>,T,rank,1,IndexF>, T>
	operator() (const IndexF &index)
	{
		typedef Encapsulate_to_Expr1<Marray<T,rank,base>,T,rank,1,IndexF> Expr1_Obj;
		return F_Expr1 < Expr1_Obj, T> (Expr1_Obj(*this,index));
	}

F_Expr1 < Encapsulate_to_Expr1<const Marray<T,rank,base>,T,rank,1,IndexF>, T>
	operator() (const IndexF &index) const
	{
		typedef Encapsulate_to_Expr1<const Marray<T,rank,base>,T,rank,1,IndexF> Expr1_Obj;
		return F_Expr1 < Expr1_Obj, T> (Expr1_Obj(*this,index));
	}


//end FEXPRS


    template < char i,int iType1>
	inline Expr1 < Encapsulate_to_Expr1<Marray<T,rank,base>,T,rank,1, Index<i,iType1> >, T, i >
	operator() (const Index < i , iType1 > &index)
	{

		typedef Encapsulate_to_Expr1<Marray<T,rank,base>,T,rank,1,Index<i,iType1> > Expr1_Obj;
		return Expr1 < Expr1_Obj, T, i > (Expr1_Obj(*this,index));
	}


	template < char i,int iType1>
	inline Expr1 < Encapsulate_to_Expr1<const Marray<T,rank,base>,T,rank,1,Index<i,iType1> >, T, i >
	operator() (const Index < i , iType1 > &index)const
	{
		typedef Encapsulate_to_Expr1<const Marray<T,rank,base>,T,rank,1, Index<i,iType1> > Expr1_Obj;
		return Expr1 < Expr1_Obj, T, i > (Expr1_Obj(*this,index));
	}



	/////////////////////////////////////////////////////////////
	// Assignation Operations from other Marray<T,rank,base> object
	/////////////////////////////////////////////////////////////
    template <class U,class base2>
	inline Marray<T,rank,base> & operator+= (const Marray<U,rank,base2> & a)
	{

		for(int i=0;i<get_dim1();i++)
			(*this)(i)+=a(i);

    return *this;
	}

    template <class U,class base2>
	inline Marray<T,rank,base> & operator-= (const Marray<U,rank,base2> & a)
	{
		for(int i=0;i<get_dim1();i++)
			(*this)(i)-=a(i);

		return *this;
	}

	 template <class U,class base2>
	inline Marray<T,rank,base> & operator/= (const Marray<U,rank,base2> & a)
	{
		for(int i=0;i<get_dim1();i++)
			(*this)(i)/=a(i);

		return *this;
	}

	 template <class U,class base2>
	inline Marray<T,rank,base> & operator*= (const Marray<U,rank,base2> & a)
	{
		for(int i=0;i<get_dim1();i++)
			(*this)(i)*=a(i);

		return *this;
	}

	//Implemented for compatibility only, shouldn't be used in performance critic sections
	//will fix this with the new cxx specification
	
	inline Marray<T,rank,base> operator-()
	{

		Marray<T,rank,base> temp(get_dim1());

		for(int i=0;i<get_dim1();i++)
			temp(i)=-(*this)(i);

		return temp;
	}

    //Scalar Assignations
    template <class U>
    inline Marray<T,rank,base> & operator= (const U &u){
        for( int i=0;i<get_dim1();i++)
            (*this)(i)=(T)u;
        return *this;


    }

    template <class U>
    inline Marray<T,rank,base> & operator+= (const U &u){
        for( int i=0;i<get_dim1();i++)
            (*this)(i)+=u;
        return *this;

    }

    template <class U>
    inline Marray<T,rank,base> & operator-= (const U &u){
        for( int i=0;i<get_dim1();i++)
            (*this)(i)-=u;
        return *this;

    }

    template <class U>
    inline Marray<T,rank,base> & operator*= (const U &u){
        for( int i=0;i<get_dim1();i++)
            (*this)(i)*=(T)u;
        return *this;

    }

    template <class U>
    inline Marray<T,rank,base> & operator/= (const U &u){
        for( int i=0;i<get_dim1();i++)
            (*this)(i)/=u;
        return *this;

    }


	//Compatibility Operations
	template <class Type>
	void fromCArray(Type *data,int n){
		this->resize(n);
		for(int i=0;i<n;i++){
			(*this)(i)=data[i];
		}

	}


	template <class Type>
	void toCArray(Type *data){
		for(int i=0;i<get_dim1();i++){
			data[i]=(*this)(i);
		}
	}


	 bool operator==(const Marray<T,rank,base> &a){
        if(get_dim1()!=a.get_dim1())
            return false;


        for(int i=0;i<get_dim1();i++){
            if( (*this)(i)!=a(i))
                return false;
                }

        return true;

        }




};



template <class Type, class base>
std::ostream & operator<< (std::ostream & os, const Marray<Type,rank,base> & v)
{
	std::cout<<"MArray1["<< v.get_dim1() << "] = " << std::endl<<"\t[ "<<std::endl;
	div_t result;
	for (int i=0;i<v.get_dim1();i++)
	{
		result = div (i, 10);
		std::cout << std::setw(8) << std::setprecision(4) <<  v(i) <<" ";
		if (result.rem == 9)
		{
			std::cout << std::endl;
		}
	}
	std::cout<<std::endl<<"\t]";
	std::cout<<std::endl;

	return os;
}


#undef CHECK
#undef rank

#endif

