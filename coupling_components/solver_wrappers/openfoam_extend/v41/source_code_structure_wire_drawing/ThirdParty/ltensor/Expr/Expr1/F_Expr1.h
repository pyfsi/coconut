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
#ifndef F_EXPR1
#define F_EXPR1
#include <assert.h>

template < class A, class T>
class F_Expr1
{
	typedef F_Expr1<A,T> IndExprA;

      private:
    A  a;


      public:


      ~F_Expr1 ()
	{

	}

	F_Expr1 (const A & arhs): a(arhs)
	{

	}

	T operator() (const int N) const
	{
		return a(N);
	}

	T & operator()(const int N)
	{
		return a(N);
	}

	int get_dim1 () const
	{
		return a.get_dim1();
	}

	inline  F_Expr1<A,T> &
	operator=(const F_Expr1<A,T> &rhs){
	    #ifdef CHECK_OOB
	    assert(this->get_dim1()==rhs.get_dim1());
	    #endif
	    Li_equals_Ri(*this,rhs);
	    return (*this);

	}

	template<class B,class U>
	inline  F_Expr1<A,T> &
	operator=(const F_Expr1<B,U> &rhs){
	    #ifdef CHECK_OOB
	    assert(this->get_dim1()==rhs.get_dim1());
	    #endif
	    Li_equals_Ri(*this,rhs);
	    return (*this);
	}

	template <class base>
	inline F_Expr1<A,T> &
	operator=(const Marray<T,1,base> &term)
	{
        #ifdef CHECK_OOB
            assert(this->get_dim1()<=term.get_dim1());
	    #endif
		Li_equals_Ri(*this,term);
		return (*this);
	}

	template <class U,class base>
	inline F_Expr1<A,T> &
	operator=(const Marray<U,1,base> &term)
	{
		#ifdef CHECK_OOB
	    		assert(this->get_dim1()<=term.get_dim1());
	    	#endif
		Li_equals_Ri(*this,term);
		return (*this);
	}


	template<class B,class U>
	inline  F_Expr1<A,T> &
	operator+=(const F_Expr1<B,U> &rhs){

    #ifdef CHECK_OOB
        assert(this->get_dim1()==rhs.get_dim1());
	#endif
        Li_plusequals_Ri(*this,rhs);
        return (*this);
	}

	template <class base>
	inline F_Expr1<A,T> &
	operator+=(const Marray<T,1,base> &term)
	{
        	#ifdef CHECK_OOB
			assert(this->get_dim1()<=term.get_dim1());
	    	#endif
		Li_plusequals_Ri(*this,term);
		return (*this);
	}

    template <class U,class base>
	inline F_Expr1<A,T> &
	operator+=(const Marray<U,1,base> &term)
	{
		#ifdef CHECK_OOB
	    		assert(this->get_dim1()<=term.get_dim1());
	    	#endif
		Li_plusequals_Ri(*this,term);
		return (*this);
	}


	inline  F_Expr1<A,T> &
	operator+=(const F_Expr1<A,T> &rhs){
	    #ifdef CHECK_OOB
	    assert(this->get_dim1()==rhs.get_dim1());
	    #endif
	    Li_plusequals_Ri(*this,rhs);
	    return (*this);

	}


	template<class B,class U>
	inline  F_Expr1<A,T> &
	operator-=(const F_Expr1<B,U> &rhs){

	 #ifdef CHECK_OOB
        assert(this->get_dim1()==rhs.get_dim1());
	#endif
        Li_minusequals_Ri(*this,rhs);
        return (*this);
	}

	template <class base>
	inline F_Expr1<A,T> &
	operator-=(const Marray<T,1,base> &term)
	{
        	#ifdef CHECK_OOB
			assert(this->get_dim1()<=term.get_dim1());
	    	#endif
		Li_minusequals_Ri(*this,term);
		return (*this);
	}

	template <class U,class base>
	inline F_Expr1<A,T> &
	operator-=(const Marray<U,1,base> &term)
	{
		#ifdef CHECK_OOB
	    		assert(this->get_dim1()<=term.get_dim1());
	    	#endif
		Li_minusequals_Ri(*this,term);
		return (*this);
	}


	inline  F_Expr1<A,T> &
	operator-=(const F_Expr1<A,T> &rhs){
	    #ifdef CHECK_OOB
	    assert(this->get_dim1()==rhs.get_dim1());
	    #endif
	    Li_minusequals_Ri(*this,rhs);
	    return (*this);

	}
	template<class U>
	inline  F_Expr1<A,T> &
	operator=(const U &rhs)
	{
		for(int n1 = 0; n1<this->get_dim1();n1++)
		{
			(*this)(n1)=rhs;
		}
		return (*this);
	}

	template<class U>
	inline  F_Expr1<A,T> &
	operator+=(const U &rhs)
	{
		for(int n1 = 0; n1<this->get_dim1();n1++)
		{
			(*this)(n1)+=rhs;
		}
		return (*this);
	}


	template <class U>
    inline const F_Expr1<A,T> &
	operator*=(const U &val){
	      for(int i=0;i<this->get_dim1();i++)
                (*this)(i)*=val;

        return *this;

	}

template<class U>
	inline  F_Expr1<A,T> &
	operator-=(const U &rhs)
	{
		for(int n1 = 0; n1<this->get_dim1();n1++)
		{
			(*this)(n1)-=rhs;
		}
		return (*this);
	}


	friend std::ostream & operator<< (std::ostream & os, const F_Expr1<A,T> & v)
	{
		int n1;
		os << std::endl << "[" << std::endl;
		for (n1 = 0; n1 < v.get_dim1(); ++n1)
		{
			os << v(n1) << " ";
		}
		os << "]" << std::endl;

		return os;
	}





};

#endif
