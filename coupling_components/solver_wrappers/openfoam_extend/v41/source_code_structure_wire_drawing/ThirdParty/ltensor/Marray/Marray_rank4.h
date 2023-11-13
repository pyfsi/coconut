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
#ifndef Marray_rank4_H
#define Marray_rank4_H

#define CHECK										\
	assert( (n1<get_dim1())&& (n2<get_dim2()) && (n3<get_dim3()) &&( n4<get_dim4()) );	\
	assert( (n1>=0) && (n2>=0) && (n3>=0) && (n4>=0) );

template < class A, class T, char i > class Expr1;
template < class A, class T, char i, char j > class Expr2;
template < class A, class T, char i, char j , char k> class Expr3;
template < class A, class T, char i, char j , char k, char l> class Expr4;
template < class A, class T> class F_Expr4;


#define rank 4

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



    Marray (long dimension1,long dimension2,long dimension3,long dimension4):base(dimension1,dimension2,dimension3,dimension4){
	    (*this)=0;

        }

    Marray (long dimension1,long dimension2,long dimension3,long dimension4, T valor):base(dimension1,dimension2,dimension3,dimension4){
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

	template <class A, class U>
	Marray<T,rank,base>& operator=(const F_Expr4 < A, U> &expr){
		for(int i=0;i<expr.get_dim1();i++)
			for(int j=0;j<expr.get_dim2();j++)
				for(int k=0;k<expr.get_dim3();k++)
                    for(int l=0;l<expr.get_dim4();l++)
                        (*this)(i,j,k,l)=expr(i,j,k,l);

		return (*this);
	}

	template <class U>
    inline Marray<T,rank,base> & operator= (const U &u){
        for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    for(int l=0;l<get_dim4();l++)
                        (*this)(i,j,k,l)=u;
        return *this;

    }


    template <class U>
    inline Marray<T,rank,base> & operator+= (const U &u){
        for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    for(int l=0;l<get_dim4();l++)
                        (*this)(i,j,k,l)+=u;
        return *this;

    }

  template <class U,class base2>
    inline Marray<T,rank,base> & operator+= (const Marray<U,rank,base2> &u){

        for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    for(int l=0;l<get_dim4();l++)
                        (*this)(i,j,k,l)+=u(i,j,k,l);
        return *this;

    }

    template <class U,class base2>
    inline Marray<T,rank,base> & operator-= (const Marray<U,rank,base2> &u){
        for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    for(int l=0;l<get_dim4();l++)
                        (*this)(i,j,k,l)-=u(i,j,k,l);
        return *this;

    }

    template <class U>
    inline Marray<T,rank,base> & operator-= (const U &u){
        for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    for(int l=0;l<get_dim4();l++)
                        (*this)(i,j,k,l)-=u;
        return *this;

    }


    template <class U>
    inline Marray<T,rank,base> & operator*= (const U &u){
        for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    for(int l=0;l<get_dim4();l++)
                        (*this)(i,j,k,l)*=u;
        return *this;


    }

		//Implemented for compatibility only, shouldn't be used in performance critic sections
	//will fix this with the new cxx specification
	
	inline Marray<T,rank,base> operator-()
	{

		Marray<T,rank,base> temp(get_dim1(),get_dim2(),get_dim3(),get_dim4());

		  for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    for(int l=0;l<get_dim4();l++)
                    temp(i,j,k,l)=-(*this)(i,j,k,l);
		
		return temp;
	}



    template <class U>
    inline Marray<T,rank,base> & operator/= (const U &u){
        for( int i=0;i<get_dim1();i++)
            for(int j=0;j<get_dim2();j++)
                for(int k=0;k<get_dim3();k++)
                    for(int l=0;l<get_dim4();l++)
                        (*this)(i,j,k,l)/=u;
        return *this;

    }

    Marray(long* dimensions):base(dimensions){
	    (*this)=0;

        }



    //copy constructor
     Marray(const Marray<T, rank,base> &R){
		(*this)=R;
        }


    inline const T operator() (const int n1,const int n2,const int n3,const int n4) const {
        #ifdef CHECK_OOB
        CHECK
        #endif

        return base::operator()(n1,n2,n3,n4);

        }

	inline T & operator() (const int n1,const int n2,const int n3,const int n4){
	    #ifdef CHECK_OOB
        CHECK
        #endif

	    return base::operator()(n1,n2,n3,n4);
	    }


	//////////////////////////////////////////////////////////////
	// MEMBER FUNCTIONS
	//////////////////////////////////////////////////////////////



	void resize(long d1,long d2,long d3,long d4){
		unsigned long  dim[4];
		dim[0]=d1;
		dim[1]=d2;
		dim[2]=d3;
		dim[3]=d4;
		base::resize(dim);
		(*this)=0;

	}

	inline int get_size() const{

        return get_dim1()*get_dim2()*get_dim3()*get_dim4();

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

    inline int get_dim4() const{
	    return base::size[3];
	    }

	inline void show_shape(){
	    std::cout<<get_dim1()<<"x"<<get_dim2()<<"x"<<get_dim3()<<"x"<<get_dim4()<<std::endl;
	    };

	friend std::ostream & operator << <T>
	(std::ostream & os, const Marray<T, rank,base> & v);


public:


	///////////////////////////////////////////////////////////////
	// create IndexedExpressions
	///////////////////////////////////////////////////////////////


	//F_EXPRs

	inline F_Expr4 < Encapsulate_to_Expr4<Marray<T,rank>,T,rank,1234,IndexF,IndexF,IndexF,IndexF>, T >
operator() (	const IndexF & index1,
				const IndexF & index2,
				const IndexF & index3,
				const IndexF & index4)
{
	typedef Encapsulate_to_Expr4<Marray<T,rank>,T,rank,1234,IndexF,IndexF,IndexF,IndexF> Expr_Obj;
	return F_Expr4<Expr_Obj,T>(Expr_Obj(*this,index1,index2,index3,index4));
}

inline F_Expr4 < Encapsulate_to_Expr4<Marray<T,rank>,T,rank,1234,IndexF,IndexF,IndexF,IndexF>, T >
operator() (	const IndexF & index1,
				const IndexF & index2,
				const IndexF & index3,
				const IndexF & index4)const
{
	typedef Encapsulate_to_Expr4<const Marray<T,rank>,T,rank,1234,IndexF,IndexF,IndexF,IndexF> Expr_Obj;
	return F_Expr4<Expr_Obj,T>(Expr_Obj(*this,index1,index2,index3,index4));
}

//To expression 3


inline F_Expr3 < Encapsulate_to_Expr3<Marray<T,rank>,T,rank,123,IndexF,IndexF,IndexF>, T >
operator() (	const IndexF & index1,
				const IndexF & index2,
				const IndexF & index3,
				const int N4)
{
	typedef Encapsulate_to_Expr3<Marray<T,rank>,T,rank,123,IndexF,IndexF,IndexF> Expr_Obj;
	return F_Expr3<Expr_Obj,T>(Expr_Obj(*this,N4,index1,index2,index3));
}


inline F_Expr3 < Encapsulate_to_Expr3<const Marray<T,rank>,T,rank,123,IndexF,IndexF,IndexF>, T >
operator() (	const IndexF & index1,
				const IndexF & index2,
				const IndexF & index3,
				const int N4)const
{
	typedef Encapsulate_to_Expr3<const Marray<T,rank>,T,rank,123,IndexF,IndexF,IndexF> Expr_Obj;
	return F_Expr3<Expr_Obj,T>(Expr_Obj(*this,N4,index1,index2,index3));
}


inline F_Expr3 < Encapsulate_to_Expr3<Marray<T,rank>,T,rank,124,IndexF,IndexF,IndexF>, T >
operator() (	const IndexF & index1,
				const IndexF & index2,
				const int N3,
				const IndexF & index4)
{
	typedef Encapsulate_to_Expr3<Marray<T,rank>,T,rank,124,IndexF,IndexF,IndexF> Expr_Obj;
	return F_Expr3<Expr_Obj,T>(Expr_Obj(*this,N3,index1,index2,index4));
}

inline F_Expr3 < Encapsulate_to_Expr3<const Marray<T,rank>,T,rank,124,IndexF,IndexF,IndexF>, T >
operator() (	const IndexF & index1,
				const IndexF & index2,
				const int N3,
				const IndexF & index4)const
{
	typedef Encapsulate_to_Expr3<const Marray<T,rank>,T,rank,124,IndexF,IndexF,IndexF> Expr_Obj;
	return F_Expr3<Expr_Obj,T>(Expr_Obj(*this,N3,index1,index2,index4));
}

inline F_Expr3 < Encapsulate_to_Expr3<Marray<T,rank>,T,rank,134,IndexF,IndexF,IndexF>, T >
operator() (	const IndexF & index1,
                const int N2,
				const IndexF & index3,
				const IndexF & index4)
{
	typedef Encapsulate_to_Expr3<Marray<T,rank>,T,rank,134,IndexF,IndexF,IndexF> Expr_Obj;
	return F_Expr3<Expr_Obj,T>(Expr_Obj(*this,N2,index1,index3,index4));
}

inline F_Expr3 < Encapsulate_to_Expr3<const Marray<T,rank>,T,rank,134,IndexF,IndexF,IndexF>, T >
operator() (	const IndexF & index1,
                const int N2,
				const IndexF & index3,
				const IndexF & index4) const
{
	typedef Encapsulate_to_Expr3<const Marray<T,rank>,T,rank,134,IndexF,IndexF,IndexF> Expr_Obj;
	return F_Expr3<Expr_Obj,T>(Expr_Obj(*this,N2,index1,index3,index4));
}


inline F_Expr3 < Encapsulate_to_Expr3<Marray<T,rank>,T,rank,234,IndexF,IndexF,IndexF>, T >
operator() (	const int N1,
                const IndexF & index2,
				const IndexF & index3,
				const IndexF & index4)
{
	typedef Encapsulate_to_Expr3< Marray<T,rank>,T,rank,234,IndexF,IndexF,IndexF> Expr_Obj;
	return F_Expr3<Expr_Obj,T>(Expr_Obj(*this,N1,index2,index3,index4));
}


inline F_Expr3 < Encapsulate_to_Expr3<const Marray<T,rank>,T,rank,234,IndexF,IndexF,IndexF>, T >
operator() (	const int N1,
                const IndexF & index2,
				const IndexF & index3,
				const IndexF & index4)const
{
	typedef Encapsulate_to_Expr3<const Marray<T,rank>,T,rank,234,IndexF,IndexF,IndexF> Expr_Obj;
	return F_Expr3<Expr_Obj,T>(Expr_Obj(*this,N1,index2,index3,index4));
}


//To expression 2


inline F_Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,12,IndexF,IndexF>, T >
operator() (	const IndexF & index1,
                const IndexF & index2,
				const int N3,
				const int N4)
{
	typedef Encapsulate_to_Expr2< Marray<T,rank>,T,rank,12,IndexF,IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,index1,index2));
}

inline F_Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,12,IndexF,IndexF>, T >
operator() (	const IndexF & index1,
                const IndexF & index2,
				const int N3,
				const int N4) const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,12,IndexF,IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N3,N4,index1,index2));
}

inline F_Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,13,IndexF,IndexF>, T >
operator() (	const IndexF & index1,
                const int N2,
                const IndexF & index3,
				const int N4)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,13,IndexF,IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N2,N4,index1,index3));
}

inline F_Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,13,IndexF,IndexF>, T >
operator() (	const IndexF & index1,
                const int N2,
                const IndexF & index3,
				const int N4) const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,13,IndexF,IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N2,N4,index1,index3));
}

inline F_Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,23,IndexF,IndexF>, T >
operator() (	const int N1,
                const IndexF & index2,
                const IndexF & index3,
				const int N4)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,23,IndexF,IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N1,N4,index2,index3));
}

inline F_Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,23,IndexF,IndexF>, T >
operator() (	const int N1,
                const IndexF & index2,
                const IndexF & index3,
				const int N4) const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,23,IndexF,IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N1,N4,index2,index3));
}

inline F_Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,24,IndexF,IndexF>, T >
operator() (	const int N1,
                const IndexF & index2,
                const int N3,
                const IndexF & index4)
{
	typedef Encapsulate_to_Expr2< Marray<T,rank>,T,rank,24,IndexF,IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N1,N3,index2,index4));
}

inline F_Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,24,IndexF,IndexF>, T >
operator() (	const int N1,
                const IndexF & index2,
                const int N3,
                const IndexF & index4) const
{
	typedef Encapsulate_to_Expr2< const Marray<T,rank>,T,rank,24,IndexF,IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N1,N3,index2,index4));
}


inline F_Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,34,IndexF,IndexF>, T >
operator() (	const int N1,
                const int N2,
                const IndexF & index3,
                const IndexF & index4)
{
	typedef Encapsulate_to_Expr2< Marray<T,rank>,T,rank,34,IndexF,IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N1,N2,index3,index4));
}

inline F_Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,34,IndexF,IndexF>, T >
operator() (	const int N1,
                const int N2,
                const IndexF & index3,
                const IndexF & index4) const
{
	typedef Encapsulate_to_Expr2< const Marray<T,rank>,T,rank,34,IndexF,IndexF> Expr_Obj;
	return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,N1,N2,index3,index4));
}

// to expression 1

inline F_Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,4,IndexF>, T >
operator() (	const int N1,
                const int N2,
                const int N3,
                const IndexF & index4)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,4,IndexF> Expr_Obj;
	return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N1,N2,N3,index4));
}


inline F_Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,4,IndexF>, T >
operator() (	const int N1,
                const int N2,
                const int N3,
                const IndexF & index4) const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,4,IndexF> Expr_Obj;
	return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N1,N2,N3,index4));
}



inline F_Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,3,IndexF>, T >
operator() (	const int N1,
                const int N2,
                const IndexF &index3,
                const int N4)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,3,IndexF> Expr_Obj;
	return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N1,N2,N4,index3));
}

inline F_Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,3,IndexF>, T >
operator() (	const int N1,
                const int N2,
                const IndexF & index3,
                const int N4) const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,3,IndexF> Expr_Obj;
	return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N1,N2,N4,index3));
}

inline F_Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,2,IndexF>, T >
operator() (	const int N1,
                const IndexF & index2,
                const int N3,
                const int N4)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,2,IndexF> Expr_Obj;
	return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N1,N3,N4,index2));
}

inline F_Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,2,IndexF>, T >
operator() (	const int N1,
                const IndexF & index2,
                const int N3,
                const int N4)const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,2,IndexF> Expr_Obj;
	return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N1,N3,N4,index2));
}

inline F_Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,1,IndexF>, T >
operator() (	const IndexF & index1,
                const int N2,
                const int N3,
                const int N4)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,1,IndexF> Expr_Obj;
	return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N2,N3,N4,index1));
}

inline F_Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,1,IndexF>, T >
operator() (	const IndexF & index1,
                const int N2,
                const int N3,
                const int N4) const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,1,IndexF> Expr_Obj;
	return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N2,N3,N4,index1));
}



	//END F_EXPRS

//checked 1/sept
//gggg

    template < char i, char j, char k, char l, int iType1, int iType2, int iType3, int iType4 >
    inline Expr4 < Encapsulate_to_Expr4<Marray<T,rank>,T,rank,1234,
					Index<i,iType1>, Index<j,iType2>, Index<k,iType3>, Index<l,iType4> >, T, i, j, k,l >
    operator() (	const Index < i,iType1 > & index1, const Index < j,iType2 > & index2,
					const Index < k, iType3 > & index3, const Index < l, iType4 > & index4 )
    {
        typedef Encapsulate_to_Expr4<Marray<T,rank>,T,rank,1234,
											Index<i,iType1>, Index<j,iType2>,
											Index<k,iType3>, Index<l,iType4> > Expr_Obj;
        return Expr4<Expr_Obj,T,i,j,k,l>(Expr_Obj(*this,index1,index2,index3,index4));
        }

    template < char i, char j, char k, char l, int iType1, int iType2, int iType3, int iType4 >
    inline Expr4 < Encapsulate_to_Expr4<Marray<T,rank>,T,rank,1234,
					Index<i,iType1>, Index<j,iType2>, Index<k,iType3>, Index<l,iType4> >, T, i, j, k,l >
    operator() (	const Index < i,iType1 > & index1, const Index < j,iType2 > & index2,
					const Index < k, iType3 > & index3, const Index < l, iType4 > & index4 )const
    {
        typedef Encapsulate_to_Expr4<Marray<T,rank>,T,rank,1234,
											Index<i,iType1>, Index<j,iType2>,
											Index<k,iType3>, Index<l,iType4> > Expr_Obj;
        return Expr4<Expr_Obj,T,i,j,k,l>(Expr_Obj(*this,index1,index2,index3,index4));
        }

//done

// Create 3D-Tensor Expressions

template < char i, char j, char k , int iType1, int iType2, int iType3>
inline Expr3 < Encapsulate_to_Expr3<Marray<T,rank>,T,rank,123,
				Index<i,iType1>, Index<j,iType2>, Index<k,iType3> >, T, i ,j, k >
operator() (const Index < i,iType1 > & index1, const Index < j,iType2 > & index2,
					const Index < k,iType3 > & index3, const int N4)
{
	typedef Encapsulate_to_Expr3<Marray<T,rank>,T,rank,123,
								Index<i,iType1>, Index<j,iType2>, Index<k,iType3> > Expr_Obj;
	return Expr3<Expr_Obj,T,i,j,k>(Expr_Obj(*this,N4,index1,index2,index3));
}

template < char i, char j, char k , int iType1, int iType2, int iType3>
inline Expr3 < Encapsulate_to_Expr3<Marray<T,rank>,T,rank,123,
				Index<i,iType1>, Index<j,iType2>, Index<k,iType3> >, T, i ,j, k >
operator() (const Index < i,iType1 > & index1, const Index < j,iType2 > & index2,
					const Index < k,iType3 > & index3, const int N4)const
{
	typedef Encapsulate_to_Expr3<Marray<T,rank>,T,rank,123,
								Index<i,iType1>, Index<j,iType2>, Index<k,iType3> > Expr_Obj;
	return Expr3<Expr_Obj,T,i,j,k>(Expr_Obj(*this,N4,index1,index2,index3));
}


//iini

template < char i, char j, char k, int iType1, int iType2, int iType4 >
inline Expr3 < Encapsulate_to_Expr3<Marray<T,rank>,T,rank,124,
				Index<i,iType1>, Index<j,iType2>, Index<k,iType4> >, T, i ,j, k >
operator() (const Index < i,iType1 > & index1, const Index < j,iType2 > & index2,
					const int N4, const Index < k,iType4> & index3 )
{
	typedef Encapsulate_to_Expr3<Marray<T,rank>,T,rank,124,
								Index<i,iType1>, Index<j,iType2>, Index<k,iType4> > Expr_Obj;
	return Expr3<Expr_Obj,T,i,j,k>(Expr_Obj(*this,N4,index1,index2,index3));
}


template < char i, char j, char k, int iType1, int iType2, int iType4 >
inline Expr3 < Encapsulate_to_Expr3<Marray<T,rank>,T,rank,124,
				Index<i,iType1>, Index<j,iType2>, Index<k,iType4> >, T, i ,j, k >
operator() (const Index < i,iType1 > & index1, const Index < j,iType2 > & index2,
					const int N4, const Index < k,iType4> & index3 )const
{
	typedef Encapsulate_to_Expr3<Marray<T,rank>,T,rank,124,
								Index<i,iType1>, Index<j,iType2>, Index<k,iType4> > Expr_Obj;
	return Expr3<Expr_Obj,T,i,j,k>(Expr_Obj(*this,N4,index1,index2,index3));
}



//done


template < char i, char j, char k, int iType1, int iType3, int iType4 >
inline Expr3 < Encapsulate_to_Expr3<Marray<T,rank>,T,rank,134,
				Index<i,iType1>, Index<j,iType3>, Index<k,iType4> >, T, i ,j, k >
operator() (const Index < i, iType1 > & index1, const int N4,
					const Index < j, iType3 > & index2, const Index < k,iType4 > & index3 )
{
	typedef Encapsulate_to_Expr3<Marray<T,rank>,T,rank,134,
								Index<i,iType1>, Index<j,iType3>, Index<k,iType4> > Expr_Obj;
	return Expr3<Expr_Obj,T,i,j,k>(Expr_Obj(*this,N4,index1,index2,index3));
}

template < char i, char j, char k, int iType1, int iType3, int iType4 >
inline Expr3 < Encapsulate_to_Expr3<Marray<T,rank>,T,rank,134,
				Index<i,iType1>, Index<j,iType3>, Index<k,iType4> >, T, i ,j, k >
operator() (const Index < i, iType1 > & index1, const int N4,
					const Index < j, iType3 > & index2, const Index < k,iType4 > & index3 )const
{
	typedef Encapsulate_to_Expr3<Marray<T,rank>,T,rank,134,
								Index<i,iType1>, Index<j,iType3>, Index<k,iType4> > Expr_Obj;
	return Expr3<Expr_Obj,T,i,j,k>(Expr_Obj(*this,N4,index1,index2,index3));
}

//done

template < char i, char j, char k, int iType2, int iType3, int iType4 >
inline Expr3 < Encapsulate_to_Expr3<Marray<T,rank>,T,rank,234,
				Index<i,iType2>, Index<j,iType3>, Index<k,iType4> >, T, i ,j, k >
operator() (const int N4, const Index < i,iType2 > & index1,
					const Index < j,iType3 > & index2, const Index < k,iType4 > & index3 )
    {
	typedef Encapsulate_to_Expr3<Marray<T,rank>,T,rank,234,
								Index<i,iType2>, Index<j,iType3>, Index<k,iType4> > Expr_Obj;
	return Expr3<Expr_Obj,T,i,j,k>(Expr_Obj(*this,N4,index1,index2,index3));
    }


template < char i, char j, char k, int iType2, int iType3, int iType4 >
inline Expr3 < Encapsulate_to_Expr3<Marray<T,rank>,T,rank,234,
				Index<i,iType2>, Index<j,iType3>, Index<k,iType4> >, T, i ,j, k >
operator() (const int N4, const Index < i,iType2 > & index1,
					const Index < j,iType3 > & index2, const Index < k,iType4 > & index3 )const
    {
	typedef Encapsulate_to_Expr3<Marray<T,rank>,T,rank,234,
								Index<i,iType2>, Index<j,iType3>, Index<k,iType4> > Expr_Obj;
	return Expr3<Expr_Obj,T,i,j,k>(Expr_Obj(*this,N4,index1,index2,index3));
    }

//done

// Create 2D-Tensor Expressions




template < char i, char j, int iType1, int iType2 >
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,12,Index<i,iType1>, Index<j,iType2> >, T, i ,j >
operator() (const Index < i, iType1 > & index1, const Index < j, iType2 > & index2,
				    const int N3, const int N4)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,12,Index<i,iType1>, Index<j,iType2> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,N4,index1, index2));
}


template < char i, char j, int iType1, int iType2 >
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,12,Index<i,iType1>, Index<j,iType2> >, T, i ,j >
operator() (const Index < i, iType1 > & index1, const Index < j, iType2 > & index2,
				    const int N3, const int N4)const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,12,Index<i,iType1>, Index<j,iType2> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,N4,index1, index2));
}
//ended


template < char i, char j, int iType1, int iType3 >
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,13,Index<i, iType1> , Index<j,iType3> >, T, i ,j >
operator() (const Index < i, iType1 > & index1, const int N3,
					const Index < j, iType3 > & index2, const int N4)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,13,Index<i, iType1> , Index<j,iType3> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,N4,index1,index2));
}



template < char i, char j, int iType1, int iType3 >
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,13,Index<i, iType1> , Index<j,iType3> >, T, i ,j >
operator() (const Index < i, iType1 > & index1, const int N3,
					const Index < j, iType3 > & index2, const int N4) const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,13,Index<i, iType1> , Index<j,iType3> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,N4,index1,index2));
}


//ended

//inni

template < char i, char j, int iType1, int iType4>
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,14,Index<i,iType1>, Index<j,iType4> >, T, i ,j >
operator() (const Index < i, iType1 > & index1, const int N3,
					const int N4, const Index < j, iType4 > & index2)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,14, Index<i,iType1>, Index<j,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,N4,index1,index2));
}


template < char i, char j, int iType1, int iType4>
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,14,Index<i,iType1>, Index<j,iType4> >, T, i ,j >
operator() (const Index < i, iType1 > & index1, const int N3,
					const int N4, const Index < j, iType4 > & index2) const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,14, Index<i,iType1>, Index<j,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,N4,index1,index2));
}

//niin
//gg

template < char i, char j, int iType2, int iType3>
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,23,Index<i,iType2>, Index<j,iType3> >, T, i ,j >
operator() (const int N3, const Index < i,iType2 > & index1,
					const Index < j,iType3 > & index2, const int N4)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,23,Index<i,iType2>, Index<j,iType3> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,N4,index1, index2));
}

template < char i, char j, int iType2, int iType3>
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,23,Index<i,iType2>, Index<j,iType3> >, T, i ,j >
operator() (const int N3, const Index < i,iType2 > & index1,
					const Index < j,iType3 > & index2, const int N4)const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,23,Index<i,iType2>, Index<j,iType3> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,N4,index1, index2));
}

//ended


//nini
//gg
template < char i, char j, int iType2, int iType4>
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,24,Index<i,iType2>, Index<j,iType4> >, T, i ,j >
operator() (const int N3, const Index < i, iType2 > & index1,
					const int N4, const Index < j, iType4 > & index2)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,24,Index<i,iType2>, Index<j,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,N4,index1, index2));
}


template < char i, char j, int iType2, int iType4>
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,24,Index<i,iType2>, Index<j,iType4> >, T, i ,j >
operator() (const int N3, const Index < i, iType2 > & index1,
					const int N4, const Index < j, iType4 > & index2)const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,24,Index<i,iType2>, Index<j,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,N4,index1, index2));
}

//ended

//nnii
template < char i, char j, int iType3, int iType4 >
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,34,Index <i, iType3>, Index<j, iType4> >, T, i ,j >
operator() (const int N3, const int N4,
					const Index < i, iType3 > & index1, const Index < j,iType4 > & index2)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,34,Index <i, iType3>, Index<j, iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,N4,index1,index2));
}

template < char i, char j, int iType3, int iType4 >
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,34,Index <i, iType3>, Index<j, iType4> >, T, i ,j >
operator() (const int N3, const int N4,
					const Index < i, iType3 > & index1, const Index < j,iType4 > & index2)const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,34,Index <i, iType3>, Index<j, iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,N3,N4,index1,index2));
}

//hh
//self contracting
//iiii
//gggg

template < char i, char j, char k, int iType1, int iType2, int iType3, int iType4 >
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,1233,
				Index<i,iType1>, Index<j,iType2>, Index<k,iType3>, Index<k,iType4> >, T, i ,j >
operator() (const Index < i,iType1 > & index1, const Index < j,iType2 > & index2,
		const Index < k,iType3 > & index3, const Index < k,iType4 > & index4)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,1233,
								Index<i,iType1>, Index<j,iType2>,
								Index<k,iType3>, Index<k,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,index1,index2,index3,index4));
}

template < char i, char j, char k, int iType1, int iType2, int iType3, int iType4 >
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,1233,
				Index<i,iType1>, Index<j,iType2>, Index<k,iType3>, Index<k,iType4> >, T, i ,j >
operator() (const Index < i,iType1 > & index1, const Index < j,iType2 > & index2,
		const Index < k,iType3 > & index3, const Index < k,iType4 > & index4)const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,1233,
								Index<i,iType1>, Index<j,iType2>,
								Index<k,iType3>, Index<k,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,index1,index2,index3,index4));
}
//hhhh


//iiii
//gggg

template < char i, char j, char k , int iType1, int iType2, int iType3, int iType4>
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,1323,
				Index<i,iType1>, Index<k,iType2>, Index<j,iType3>, Index<k,iType4> >, T, i ,j >
operator() (const Index < i,iType1 > & index1, const Index < k,iType2 > & index3,
					const Index < j,iType3 > & index2, const Index < k,iType4 > & index4)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,1323,
								Index<i,iType1>, Index<k,iType2>,
								Index<j,iType3>, Index<k,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,index1,index3,index2,index4));
}


template < char i, char j, char k , int iType1, int iType2, int iType3, int iType4>
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,1323,
				Index<i,iType1>, Index<k,iType2>, Index<j,iType3>, Index<k,iType4> >, T, i ,j >
operator() (const Index < i,iType1 > & index1, const Index < k,iType2 > & index3,
					const Index < j,iType3 > & index2, const Index < k,iType4 > & index4) const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,1323,
								Index<i,iType1>, Index<k,iType2>,
								Index<j,iType3>, Index<k,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,index1,index3,index2,index4));
}

//iiii

//gggg
template < char i, char j, char k, int iType1, int iType2, int iType3, int iType4 >
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,1332,
				Index<i,iType1>, Index<k,iType2>, Index<k,iType3>, Index<j,iType4> >, T, i ,j >
operator() (const Index < i, iType1 > & index1, const Index < k, iType2 > & index3,
				const Index < k, iType3 > & index4, const Index < j, iType4 > & index2)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,1332,
								Index<i,iType1>, Index<k,iType2>,
								Index<k,iType3>, Index<j,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,index1,index3,index4,index2));
}

template < char i, char j, char k, int iType1, int iType2, int iType3, int iType4 >
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,1332,
				Index<i,iType1>, Index<k,iType2>, Index<k,iType3>, Index<j,iType4> >, T, i ,j >
operator() (const Index < i, iType1 > & index1, const Index < k, iType2 > & index3,
				const Index < k, iType3 > & index4, const Index < j, iType4 > & index2)const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,1332,
								Index<i,iType1>, Index<k,iType2>,
								Index<k,iType3>, Index<j,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,index1,index3,index4,index2));
}


//ended

//iiii
//gggg

template < char i, char j, char k , int iType1, int iType2, int iType3, int iType4>
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,3123,
				Index<k,iType1>, Index<i,iType2>, Index<j,iType3>, Index<k,iType4> >, T, i ,j >
operator() (const Index < k, iType1 > & index3,const Index < i, iType2 > & index1,
				 const Index < j, iType3 > & index2, const Index < k, iType4 > & index4)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,3123,
								Index<k,iType1>, Index<i,iType2>,
								Index<j,iType3>, Index<k,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,index3,index1,index2,index4));
}

template < char i, char j, char k , int iType1, int iType2, int iType3, int iType4>
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,3123,
				Index<k,iType1>, Index<i,iType2>, Index<j,iType3>, Index<k,iType4> >, T, i ,j >
operator() (const Index < k, iType1 > & index3,const Index < i, iType2 > & index1,
				 const Index < j, iType3 > & index2, const Index < k, iType4 > & index4)const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,3123,
								Index<k,iType1>, Index<i,iType2>,
								Index<j,iType3>, Index<k,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,index3,index1,index2,index4));
}


//done

//iiii

//gggg

template < char i, char j, char k, int iType1, int iType2, int iType3, int iType4 >
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,3132,
				Index<k,iType1>, Index<i,iType2>, Index<k, iType3>, Index<j,iType4> >, T, i ,j >
operator() (const Index < k, iType1 > & index3,const Index < i, iType2 > & index1,
				 const Index < k, iType3 > & index4, const Index < j, iType4 > & index2)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,3132,
								Index<k,iType1>, Index<i,iType2>, Index<k, iType3>, Index<j,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this),index3,index1,index4,index2);
}


template < char i, char j, char k, int iType1, int iType2, int iType3, int iType4 >
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,3132,
				Index<k,iType1>, Index<i,iType2>, Index<k, iType3>, Index<j,iType4> >, T, i ,j >
operator() (const Index < k, iType1 > & index3,const Index < i, iType2 > & index1,
				 const Index < k, iType3 > & index4, const Index < j, iType4 > & index2) const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,3132,
								Index<k,iType1>, Index<i,iType2>, Index<k, iType3>, Index<j,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this),index3,index1,index4,index2);
}


//iiii


////gggg

template < char i, char j, char k , int iType1, int iType2, int iType3, int iType4>
inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank>,T,rank,3312,
				Index<k,iType1>, Index<k,iType2>, Index<i,iType3>, Index<j,iType4> >, T, i ,j >
operator() (const Index < k , iType1> & index3, const Index < k,iType2 > & index4,
				 const Index < i, iType3 > & index1, const Index < j, iType4 > & index2)
{
	typedef Encapsulate_to_Expr2<Marray<T,rank>,T,rank,3312,Index<k,iType1>, Index<k,iType2>,
								Index<i,iType3>, Index<j,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,index3,index4,index1,index2));
}



template < char i, char j, char k , int iType1, int iType2, int iType3, int iType4>
inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,3312,
				Index<k,iType1>, Index<k,iType2>, Index<i,iType3>, Index<j,iType4> >, T, i ,j >
operator() (const Index < k , iType1> & index3, const Index < k, iType2 > & index4,
				 const Index < i, iType3 > & index1, const Index < j, iType4 > & index2)const
{
	typedef Encapsulate_to_Expr2<const Marray<T,rank>,T,rank,3312,Index<k,iType1>, Index<k,iType2>,
								Index<i,iType3>, Index<j,iType4> > Expr_Obj;
	return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,index3,index4,index1,index2));
}


// Create 1D-Tensor Expressions
//

template < char i, int iType1>
inline Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,1,Index<i,iType1> >, T, i>
operator() (const Index<i,iType1>  & index1, const int N2,
					const int N3, const int N4)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,1,Index<i,iType1> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N2,N3,N4,index1));
}

template < char i, int iType1>
inline Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,1,Index<i,iType1> >, T, i>
operator() (const Index<i,iType1>  & index1, const int N2,
					const int N3, const int N4)const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,1,Index<i,iType1> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N2,N3,N4,index1));
}




template < char i, int iType1>
inline Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,2,Index<i,iType1> >, T, i>
operator() ( const int N2,const Index<i,iType1> & index1,
					const int N3, const int N4)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,2,Index<i,iType1> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N2,N3,N4, index1));
}

template < char i, int iType1>
inline Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,2,Index<i,iType1> >, T, i>
operator() ( const int N2,const Index<i,iType1> & index1,
					const int N3, const int N4)const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,2,Index<i,iType1> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N2,N3,N4, index1));
}


template < char i, int iType>
inline Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,3,Index<i,iType> >, T, i>
operator() (  const int N2, const int N3,
					const Index < i, iType > & index1, const int N4)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,3,Index < i, iType > > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N2,N3,N4, index1));
}

template < char i, int iType>
inline Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,3,Index<i,iType> >, T, i>
operator() (  const int N2, const int N3,
					const Index < i, iType > & index1, const int N4)const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,3,Index < i, iType > > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N2,N3,N4, index1));
}

template < char i, int iType>
inline Expr1 < Encapsulate_to_Expr1<Marray<T,rank>,T,rank,4,Index<i,iType> >, T, i>
operator() ( const int N2, const int N3,
					const int N4, const Index<i,iType> & index1)
{
	typedef Encapsulate_to_Expr1<Marray<T,rank>,T,rank,4,Index<i,iType> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N2,N3,N4,index1));
}

template < char i, int iType>
inline Expr1 < Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,4,Index<i,iType> >, T, i>
operator() ( const int N2, const int N3,
					const int N4, const Index<i,iType> & index1)const
{
	typedef Encapsulate_to_Expr1<const Marray<T,rank>,T,rank,4,Index<i,iType> > Expr_Obj;
	return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N2,N3,N4,index1));
}

template <class Type>
void fromCArray(Type *data, int dim1, int dim2, int dim3,int dim4){

	resize(dim1,dim2,dim3,dim4);

	for(int i=0;i<dim1;i++)
		for(int j=0;j<dim2;j++)
			for(int k=0;k<dim3;k++)
				for(int l=0;l<dim4;l++)
					(*this)(i,j,k,l)=data[i*dim1*dim2*dim3+j*dim3*dim2+k*dim1+l];

}

template <class Type>
void toCArray(Type *data){

	int dim1=get_dim1();
	int dim2=get_dim2();
	int dim3=get_dim3();
	int dim4=get_dim4();
	for(int i=0;i<dim1;i++)
		for(int j=0;j<dim2;j++)
			for(int k=0;k<dim3;k++)
				for(int l=0;l<dim4;l++)
					data[i*dim1*dim2*dim3+j*dim3*dim2+k*dim1+l]=(*this)(i,j,k,l);

}


void sucesion(int init,int stride=1)
	{
		int count=init;

		for(int i=0;i<get_dim1();i++)
		{
			for(int j=0;j<get_dim2();j++)
			{
				for(int k=0;k<get_dim3();k++){
					for(int l=0;l<get_dim4();l++){
		    				(*this)(i,j,k,l)=count;
							count = count + stride;
					}
				}

			}
		}
	}




	void convert_to_Tensor2(Marray<T,2> & v)
	{
		int dim1=get_dim1();
		int dim2=get_dim2();
		int dim3=get_dim3();
		int dim4=get_dim4();
		for (int i=0;i<dim1;i++)
		{
			for (int j=0;j<dim2;j++)
			{
				for (int k=0;k<dim3;k++)
				{
					for (int l=0;l<dim4;l++)
					{
						v(i*dim2+j,k*dim4+l)=(*this)(i,j,k,l);
					}
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

        if(get_dim4()!=a.get_dim4())
            return false;

        for(int i=0;i<get_dim1();i++){
            for(int j=0;j<get_dim2();j++){
                for(int k=0;k<get_dim3();k++){
                    for(int l=0;l<get_dim4();l++)
                        if( (*this)(i,j,k,l)!=a(i,j,k,l))
                            return false;
                    }
                }

            }
        return true;

        }




};




template <class Type, class base>
std::ostream & operator<< (std::ostream & os,const  Marray<Type,rank,base> & v){

        std::cout<<"MArray4["<< v.get_dim1() <<"," << v.get_dim2() <<"," << v.get_dim3() <<"," << v.get_dim4()<<"] = " << std::endl <<"\t[ "<<std::endl;

		for( int i=0;i<v.get_dim1();i++)
		{
			for(int j=0;j<v.get_dim2();j++)
			{
			        std::cout<<"index("<<i<<","<<j<<")" << std::endl << "\t["<<std::endl;
				for(int k=0;k<v.get_dim3();k++)
				{
					for(int l=0;l<v.get_dim4();l++)
					{
						std::cout << std::setw(8) << std::setprecision(4)<< v(i,j,k,l) <<" ";
					}
					std::cout<<std::endl;
				}

				std::cout<<"\t]"<<std::endl;
			}

		}




    return os;

    }





#undef rank

#endif
