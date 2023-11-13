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
#ifndef Marray_rank3_H
#define Marray_rank3_H

//#define USE_ASSERT_Marray
#ifdef USE_ASSERT_Marray
#include <assert.h>
#endif

#define CHECK										\
	assert( (n1<get_dim1())&& (n2<get_dim2()) && (n3<get_dim3()) );	\
	assert( (n1>=0) && (n2>=0) && (n3>=0) );

#include <algorithm>
#include <iostream>


template < class A, class T, char i > class Expr1;
template < class A, class T, char i, char j > class Expr2;
template < class A, class T, char i, char j , char k> class Expr3;
template < class A, class T, char i, char j , char k, char l> class Expr4;
template < class A, class T> class F_Expr3;


#define rank 3



//
// FUNTIONS DECLARATIONS in order to FRIEND Function "ostream <<"  to work
//
template < class T,class base> class Marray <T,rank,base>;

template<class T,class base>
std::ostream & operator << (std::ostream & os, const Marray<T, rank,base> & v);


//
// Class Definition
//
template < class T,class base>
class Marray <T,rank,base> : public base
{


    ////////////////
    //Constructors
    ///////////////
    public:



    Marray (long dimension1,long dimension2,long dimension3):base(dimension1,dimension2,dimension3){
	    (*this)=(T)0;
        }

    Marray (long dimension1,long dimension2,long dimension3, T valor):base(dimension1,dimension2,dimension3){
        *this=valor;
        }

    Marray(){

        };

	template <class fType,class b>
	Marray<T,rank,base>& operator=(const Marray<fType,rank,b> &rterm){

	  	base::operator=(*dynamic_cast<const b*>(&rterm));
       return *this;
    }

	Marray<T,rank,base>& operator=(const Marray<T,rank,base> &rterm){

	   base::operator=(*dynamic_cast<const base*>(&rterm));
		return *this;

    }

	 //template <class U,class base2, int rank2, int perm>
	//Marray<T,rank,base>& operator=(const F_Expr3 < Encapsulate_to_Expr3<Marray<U,rank2,base2>,U,rank2,perm>, U> &expr){
	template <class A,class U>
	Marray<T,rank,base>& operator=(const F_Expr3 <A,U> &expr){
		for(int i=0;i<expr.get_dim1();i++)
			for(int j=0;j<expr.get_dim2();j++)
				for(int k=0;k<expr.get_dim3();k++)
					(*this)(i,j,k)=expr(i,j,k);

		return (*this);
	}

	template <class U>
    inline Marray<T,rank,base> & operator= (const U &u){
        for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    (*this)(i,j,k)=u;
        return *this;

    }

    template <class U>
    inline Marray<T,rank,base> & operator+= (const U &u){
        for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    (*this)(i,j,k)+=u;
        return *this;

    }


    template <class U,class base2>
    inline Marray<T,rank,base> & operator+= (const Marray<U,rank,base2> &u){
        for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    (*this)(i,j,k)+=u(i,j,k);
        return *this;

    }


        template <class U,class base2>
    inline Marray<T,rank,base> & operator-= (const Marray<U,rank,base2> &u){
        for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    (*this)(i,j,k)-=u(i,j,k);
        return *this;

    }


    template <class U>
    inline Marray<T,rank,base> & operator-= (const U &u){
        for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    (*this)(i,j,k)-=u;
        return *this;

    }

    template <class U>
    inline Marray<T,rank,base> & operator*= (const U &u){
        for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    (*this)(i,j,k)*=u;
        return *this;

    }


    template <class U>
    inline Marray<T,rank,base> & operator/= (const U &u){
        for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    (*this)(i,j,k)/=u;
        return *this;

    }


	
	//Implemented for compatibility only, shouldn't be used in performance critic sections
	//will fix this with the new cxx specification
	
	inline Marray<T,rank,base> operator-()
	{

		Marray<T,rank,base> temp(get_dim1(),get_dim2(),get_dim3());

		 for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    temp(i,j,k)=-(*this)(i,j,k);
		
		return temp;
	}



    Marray(long* dimensions):base(dimensions){
	    (*this)=0;

        }




    //copy constructor
     Marray(const Marray<T, rank,base> &R){
			(*this)=R;
        }




    inline const T operator() (const int n1,const int n2,const int n3) const {
        #ifdef CHECK_OOB
        CHECK
        #endif

        return base::operator()(n1,n2,n3);

        }

	inline T & operator() (const int n1,const int n2,const int n3){
	    #ifdef CHECK_OOB
        CHECK
        #endif

	    return base::operator()(n1,n2,n3);
	    }


    void LeviCivita()
	{
		resize(3,3,3);
		(*this)(0,1,2)=1;
		(*this)(2,0,1)=1;
		(*this)(1,2,0)=1;
		(*this)(2,1,0)=-1;
		(*this)(0,2,1)=-1;
		(*this)(1,0,2)=-1;
	}


	//////////////////////////////////////////////////////////////
	// MEMBER FUNCTIONS
	//////////////////////////////////////////////////////////////



	void resize(long d1,long d2,long d3){
		unsigned long  dim[3];
		dim[0]=d1;
		dim[1]=d2;
		dim[2]=d3;
		base::resize(dim);
		(*this)=(T)0;

	}

	inline int get_size() const{

        return base::size[0]*base::size[1];

	}


	inline int get_dim1() const{
	    return base::size[0];

	    }

	inline int get_dim2() const{
	    return base::size[1];
	    }
	inline int get_dim3() const{
	    return base::size[2];
	    }

	inline void show_shape(){
	    std::cout<<base::size[0]<<"x"<<base::size[1]<<std::endl;
	    };

	friend std::ostream & operator << <T>
	(std::ostream & os, const Marray<T, rank,base> & v);


public:


	///////////////////////////////////////////////////////////////
	// create IndexedExpressions
	///////////////////////////////////////////////////////////////

	//F_EXPRS

	inline F_Expr3 < Encapsulate_to_Expr3<Marray<T,rank>,T,rank,123,IndexF,IndexF,IndexF>, T >
operator() (	const IndexF & index1,
				const IndexF & index2,
				const IndexF & index3 )
{
	typedef Encapsulate_to_Expr3<Marray<T,rank>,T,rank,123,IndexF,IndexF,IndexF> Expr_Obj;
	return F_Expr3<Expr_Obj,T>(Expr_Obj(*this,index1,index2,index3));
}

inline F_Expr3 < Encapsulate_to_Expr3<const Marray<T,rank>,T,rank,123,IndexF,IndexF,IndexF>, T >
operator() (	const IndexF & index1,
				const IndexF & index2,
				const IndexF & index3 ) const
{
	typedef Encapsulate_to_Expr3<const Marray<T,rank>,T,rank,123,IndexF,IndexF,IndexF> Expr_Obj;
	return F_Expr3<Expr_Obj,T>(Expr_Obj(*this,index1,index2,index3));
}


inline F_Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,12,IndexF,IndexF>, T >
operator() (const IndexF & index1,const IndexF & index2, const int N3)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,12,IndexF,IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N3,index1,index2));
}

inline F_Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,12,IndexF,IndexF>, T >
operator() (const IndexF & index1,const IndexF & index2, const int N3) const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,12,IndexF, IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N3,index1,index2));
}


inline F_Expr2 < Encapsulate_to_Expr2< Marray<T,rank>,T,rank,13,IndexF,IndexF>, T >
operator() (const IndexF  & index1, const int N3,const IndexF & index2)
{
	typedef Encapsulate_to_Expr2< Marray<T,rank>,T,rank,13,IndexF,IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N3,index1,index2));
}

inline F_Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,13,IndexF,IndexF>, T>
operator() (const IndexF & index1, const int N3,const IndexF & index2) const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,13,IndexF, IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N3,index1,index2));
}
inline F_Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,23,IndexF,IndexF>, T>
operator() (const int N3, const IndexF & index1, const IndexF  & index2)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,23,IndexF,IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N3),index1,index2);
}

inline F_Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,23,IndexF,IndexF>, T>
operator() (const int N3, const IndexF & index1, const IndexF  & index2)const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,23,IndexF,IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N3,index1,index2));
}

inline F_Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,1,IndexF,IndexF>, T>
operator() (const IndexF & index1, const int N2, const int N3)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,1,IndexF,IndexF> Expr_Obj;
	return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N2,N3,index1));
}


inline F_Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,1,IndexF>, T>
operator() (const IndexF & index1, const int N2, const int N3) const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,1,IndexF> Expr_Obj;
	return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N2,N3,index1));
}

inline F_Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,2,IndexF>, T>
operator() (const int N2, const IndexF  & index1, const int N3)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,2,IndexF> Expr_Obj;
	return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N2,N3,index1));
}


inline F_Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,2,IndexF>, T>
operator() (const int N2, const IndexF  & index1, const int N3) const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,2,IndexF> Expr_Obj;
	return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N2,N3,index1));
}

inline F_Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,3,IndexF>, T>
operator() (const int N2, const int N3, const IndexF & index1)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,3,IndexF> Expr_Obj;
	return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N2,N3,index1));
}

inline F_Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,3,IndexF>, T>
operator() (const int N2, const int N3, const IndexF & index1) const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,3,IndexF> Expr_Obj;
	return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N2,N3,index1));
}

	//End F_EXPRS

	//three indexes
//ggg
template <char i, char j, char k, int iType1, int iType2, int iType3 >
inline Expr3 <	Encapsulate_to_Expr3<Marray<T,rank>,T,rank,
				123,Index<i,iType1> ,Index<j,iType2>, Index<k,iType3> >,
				T, i, j, k >
operator() (	const Index < i, iType1 > & index1,
				const Index < j, iType2 > & index2,
				const Index < k, iType3> & index3 )
{
	typedef Encapsulate_to_Expr3<Marray<T,rank>,T,rank,
			123,Index<i,iType1> ,Index<j,iType2>, Index<k,iType3> > Expr_Obj;
	return Expr3<Expr_Obj,T,i,j,k>(Expr_Obj(*this,index1,index2,index3));
}

template <char i, char j, char k, int iType1, int iType2, int iType3 >
inline Expr3 <	Encapsulate_to_Expr3<const Marray<T,rank>,T,rank,
				123,Index<i,iType1> ,Index<j,iType2>, Index<k,iType3> >,
				T, i, j, k >
operator() (	const Index < i, iType1 > & index1,
				const Index < j, iType2 > & index2,
				const Index < k, iType3> & index3 )const
{
	typedef Encapsulate_to_Expr3<const Marray<T,rank>,T,rank,
			123,Index<i,iType1> ,Index<j,iType2>, Index<k,iType3> > Expr_Obj;
	return Expr3<Expr_Obj,T,i,j,k>(Expr_Obj(*this,index1,index2,index3));
}

//2 indexes


template <char i, char j , int iType1, int iType2>
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,12,Index<i,iType1>, Index<j,iType2> >, T, i ,j >
operator() (const Index < i,iType1 > & index1,const Index < j, iType2 > & index2, const int N3)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,12, Index<i, iType1> , Index<j,iType2> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,index1,index2));
}

template <char i, char j , int iType1, int iType2>
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,12,Index<i,iType1>, Index<j,iType2> >, T, i ,j >
operator() (const Index < i,iType1 > & index1,const Index < j, iType2 > & index2, const int N3)const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,12, Index<i, iType1> , Index<j,iType2> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,index1,index2));
}

//ended 2 indexes and const
//starting index const index

template <char i, char j , int iType1 , int iType3>
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,13,Index<i,iType1>, Index<j,iType3> >, T, i ,j >
operator() (const Index < i,iType1 > & index1, const int N3,const Index < j,iType3 > & index2)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,13,Index<i,iType1>, Index<j,iType3> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,index1,index2));
}

template <char i, char j , int iType1 , int iType3>
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,13,Index<i,iType1>, Index<j,iType3> >, T, i ,j >
operator() (const Index < i,iType1 > & index1, const int N3,const Index < j,iType3 > & index2) const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,13,Index<i,iType1>, Index<j,iType3> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,index1,index2));
}

//ended index const index

template <char i, char j , int iType2, int iType3>
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,23,Index<i, iType2> , Index<j,iType3> >, T, i ,j >
operator() (const int N3, const Index < i,iType2 > & index1, const Index < j,iType3 > & index2)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,23,Index<i, iType2> , Index<j,iType3> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,index1,index2));
}

template <char i, char j , int iType2, int iType3>
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,23,Index<i, iType2> , Index<j,iType3> >, T, i ,j >
operator() (const int N3, const Index < i,iType2 > & index1, const Index < j,iType3 > & index2)const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,23,Index<i, iType2> , Index<j,iType3> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,index1,index2));
}


//ended ended constant index index

//to rank 1

//starting index constant constant

template<char i, int iType1>
inline Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,1,Index<i,iType1> >, T, i>
operator() (const Index < i,iType1 > & index1, const int N2, const int N3)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,1,Index<i,iType1> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N2,N3,index1));
}

template<char i, int iType1>
inline Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,1,Index<i,iType1> >, T, i>
operator() (const Index < i,iType1 > & index1, const int N2, const int N3) const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,1,Index<i,iType1> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N2,N3,index1));
}


template <char i, int iType2>
inline Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,2,Index<i,iType2> >, T, i>
operator() (const int N2, const Index < i,iType2 > & index1, const int N3)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,2,Index<i,iType2> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N2,N3,index1));
}


template <char i, int iType2>
inline Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,2,Index<i,iType2> >, T, i>
operator() (const int N2, const Index < i,iType2 > & index1, const int N3)const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,2,Index<i,iType2> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N2,N3,index1));
}



template <char i, int iType3>
inline Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,3,Index<i,iType3> >, T, i>
operator() (const int N2, const int N3, const Index < i ,iType3> & index1)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,3,Index<i,iType3> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N2,N3,index1));
}

template <char i, int iType3>
inline Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,3,Index<i,iType3> >, T, i>
operator() (const int N2, const int N3, const Index < i ,iType3> & index1) const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,3,Index<i,iType3> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N2,N3,index1));
}

//ended this stuff

//now start repeating indexes


//ijj
//ggg
template <char i, char j , int iType1, int iType2, int iType3>
inline Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,122,
				Index <i,iType1> , Index<j,iType2> , Index<j,iType3> >, T, i>
operator() (	const Index < i, iType1 > & index1,
				const Index < j, iType2> & index2,
				const Index < j, iType3 > & index3)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,122,
								Index <i,iType1> , Index<j,iType2> , Index<j,iType3> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,index1,index2,index3));
}

template <char i, char j , int iType1, int iType2, int iType3>
inline Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,122,
				Index <i,iType1> , Index<j,iType2> , Index<j,iType3> >, T, i>
operator() (	const Index < i, iType1 > & index1,
				const Index < j, iType2> & index2,
				const Index < j, iType3 > & index3) const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,122,
								Index <i,iType1> , Index<j,iType2> , Index<j,iType3> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,index1,index2,index3));
}

//ended ijj

//start jij


template <char i, char j , int iType1, int iType2, int iType3>
inline Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,212,
				Index <j,iType1> , Index<i, iType2> , Index<j,iType3> >, T, i>
operator() (	const Index < j,iType1 > & index1,
				const Index < i,iType2 > & index2,
				const Index < j,iType3 > & index3)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,212,
								Index <j,iType1> , Index<i, iType2> , Index<j,iType3> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,index1,index2,index3));
}

template <char i, char j , int iType1, int iType2, int iType3>
inline Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,212,
				Index <j,iType1> , Index<i, iType2> , Index<j,iType3> >, T, i>
operator() (	const Index < j,iType1 > & index1,
				const Index < i,iType2 > & index2,
				const Index < j,iType3 > & index3)const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,212,
								Index <j,iType1> , Index<i, iType2> , Index<j,iType3> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,index1,index2,index3));
}

//ended jij

//jji
template <char i, char j, int iType1, int iType2, int iType3>
inline Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,221,
				Index<j,iType1>, Index<j,iType2>, Index<i,iType3> >, T, i>
operator() (	const Index < j,iType1 > & index1,
				const Index < j,iType2 > & index2,
				const Index < i,iType3 > & index3)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,221,
								Index<j,iType1>, Index<j,iType2>, Index<i,iType3> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,index1,index2,index3));
}

template <char i, char j, int iType1, int iType2, int iType3>
inline Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,221,
				Index<j,iType1>, Index<j,iType2>, Index<i,iType3> >, T, i>
operator() (	const Index < j,iType1 > & index1,
				const Index < j,iType2 > & index2,
				const Index < i,iType3 > & index3)const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,221,
								Index<j,iType1>, Index<j,iType2>, Index<i,iType3> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,index1,index2,index3));
}

//done

template <class Type>
void fromCArray(Type *data, int dim1, int dim2, int dim3){

	resize(dim1,dim2,dim3);

	for(int i=0;i<dim1;i++)
		for(int j=0;j<dim2;j++)
			for(int k=0;k<dim3;k++)
				(*this)(i,j,k)=data[i*dim1*dim2+j*dim2+k];

}

template <class Type>
void toCArray(Type *data){

	int dim1=get_dim1();
	int dim2=get_dim2();
	int dim3=get_dim3();
	for(int i=0;i<dim1;i++)
		for(int j=0;j<dim2;j++)
			for(int k=0;k<dim3;k++)
				data[i*dim1*dim2+j*dim2+k]=(*this)(i,j,k);

}



    	void sucesion(int init,int stride=1)
	{
		int count=init;

		for(int i=0;i<get_dim1();i++)
		{
		for(int j=0;j<get_dim2();j++)
		{
			for(int k=0;k<get_dim3();k++){
		    			(*this)(i,j,k)=count;
						count = count + stride;
			}

		}
		}
	}



    bool operator==(const Marray<T,rank,base> &a){
        if(get_dim1()!=a.get_dim1())
            return false;
        if(get_dim2()!=a.get_dim2())
            return false;
        if(get_dim3()!=a.get_dim3())
            return false;

        for(int i=0;i<get_dim1();i++){
            for(int j=0;j<get_dim2();j++){
                for(int k=0;k<get_dim3();k++){
                    if( (*this)(i,j,k)!=a(i,j,k))
                        return false;
                    }
                }

            }
        return true;

        }


};




template <class Type, class base>
std::ostream & operator<< (std::ostream & os,const  Marray<Type,rank,base> & v){


        std::cout<<"MArray3["<< v.get_dim1() <<"," << v.get_dim2() <<"," << v.get_dim3()<<"] = " << std::endl <<"\t[ "<<std::endl;

		for( int i=0;i<v.get_dim1();i++)
		{
			std::cout<<"index("<<i<<")"<< std::endl << "\t["<<std::endl;
			for(int j=0;j<v.get_dim2();j++){
				for(int k=0;k<v.get_dim3();k++)
					std::cout << std::setw(8) << std::setprecision(4)<< v(i,j,k) <<" ";
				std::cout<<std::endl;
			}
			std::cout<<"\t]"<<std::endl;
		}




    return os;


}

#undef CHECK
#undef rank

#endif
