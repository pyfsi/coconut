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
#ifndef F_EXPR2
#define F_EXPR2
template < class A, class T>
class F_Expr2
{


	typedef F_Expr2<A,T> IndExprA;

	private:

	A  a;



      public:

    ~F_Expr2 ()
	{
	}

	  F_Expr2 (const A & arhs): a(arhs)
	{

	}


	T operator() (const int N1, const int N2) const
	{

		return a(N1,N2);
	}

	T & operator()(const int N1, const int N2)
	{

		return a(N1,N2);
	}

	int get_dim1 () const
	{
		return a.get_dim1();
	}

	int get_dim2 () const
	{
		return a.get_dim2();
	}

	// T(i,j) = B(i,j)....
	//
	template<class B,class U>
	inline const F_Expr2<A,T> &
	operator=(const F_Expr2<B,U> &rhs){
	 #ifdef CHECK_OOB
        assert(this->get_dim1()==rhs.get_dim1());
        assert(this->get_dim2()==rhs.get_dim2());
	 #endif
        Lij_equals_Rij(*this,rhs);
        return (*this);

	}

	inline const F_Expr2<A,T> &
	operator=(const F_Expr2<A,T> &rhs){
	    #ifdef CHECK_OOB
        assert(this->get_dim1()==rhs.get_dim1());
        assert(this->get_dim2()==rhs.get_dim2());
	 #endif
        Lij_equals_Rij(*this,rhs);
        return (*this);
	}


    template<class B,class U>
    inline const F_Expr2<A,T> &
	operator+=(const F_Expr2<B,U> &rhs){
	    #ifdef CHECK_OOB
        assert(this->get_dim1()==rhs.get_dim1());
        assert(this->get_dim2()==rhs.get_dim2());
        #endif
        Lij_plusequals_Rij(*this,rhs);
        return (*this);
	}

	template <class base>
	inline F_Expr2<A,T> &
	operator+=(const Marray<T,2,base> &term){

        #ifdef CHECK_OOB
	    assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2()
              );
	    #endif
        Lij_plusequals_Rij(*this,term);
        return *this;
	}

	template <class U,class base>
	 F_Expr2<A,T> &
	operator+=(const Marray<U,2,base> &term){
        #ifdef CHECK_OOB
	    assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2()
              );
	    #endif
        Lij_plusequals_Rij(*this,term);
        return *this;
	}

		inline const F_Expr2<A,T> &
	operator+=(const F_Expr2<A,T> &rhs){
	    #ifdef CHECK_OOB
        assert(this->get_dim1()==rhs.get_dim1());
        assert(this->get_dim2()==rhs.get_dim2());
	 #endif
        Lij_plusequals_Rij(*this,rhs);
        return (*this);
	}


	template<class B,class U>
    inline const F_Expr2<A,T> &
	operator-=(const F_Expr2<B,U> &rhs){
	    #ifdef CHECK_OOB
        assert(this->get_dim1()==rhs.get_dim1());
        assert(this->get_dim2()==rhs.get_dim2());
        #endif
        Lij_minusequals_Rij(*this,rhs);
        return (*this);
	}

	template <class base>
	inline F_Expr2<A,T> &
	operator-=(const Marray<T,2,base> &term){

        #ifdef CHECK_OOB
	    assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2()
              );
	    #endif
        Lij_minusequals_Rij(*this,term);
        return *this;
	}

	template <class U,class base>
	 F_Expr2<A,T> &
	operator-=(const Marray<U,2,base> &term){
        #ifdef CHECK_OOB
	    assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2()
              );
	    #endif
        Lij_minusequals_Rij(*this,term);
        return *this;
	}

	inline const F_Expr2<A,T> &
	operator-=(const F_Expr2<A,T> &rhs){
	    #ifdef CHECK_OOB
        assert(this->get_dim1()==rhs.get_dim1());
        assert(this->get_dim2()==rhs.get_dim2());
	 #endif
        Lij_minusequals_Rij(*this,rhs);
        return (*this);
	}


	template <class base>
	inline F_Expr2<A,T> &
	operator=(const Marray<T,2,base> &term){

        #ifdef CHECK_OOB
	    assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2()
              );
	    #endif
        Lij_equals_Rij(*this,term);
        return *this;
	}

	template <class U,class base>
	 F_Expr2<A,T> &
	operator=(const Marray<U,2,base> &term){
        #ifdef CHECK_OOB
	    assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2()
              );
	    #endif
        Lij_equals_Rij(*this,term);
        return *this;
	}


template <class U>
	 F_Expr2<A,T> &
	operator=(const U &term){
	    for(int i=0;i<this->get_dim1();i++)
            for(int j=0;j<this->get_dim2();j++)
                (*this)(i,j)=term;

        return *this;
	}

	template <class U>
	 F_Expr2<A,T> &
	operator+=(const U &term){
	    for(int i=0;i<this->get_dim1();i++)
            for(int j=0;j<this->get_dim2();j++)
                (*this)(i,j)+=term;

        return *this;
	}


	template <class U>
	 F_Expr2<A,T> &
	operator-=(const U &term){
	    for(int i=0;i<this->get_dim1();i++)
            for(int j=0;j<this->get_dim2();j++)
                (*this)(i,j)-=term;

        return *this;
	}


	template <class U>
    inline const F_Expr2<A,T> &
	operator*=(const U &val){
	      for(int i=0;i<this->get_dim1();i++)
            for(int j=0;j<this->get_dim2();j++)
                    (*this)(i,j)*=val;

        return *this;

	}


	friend  std::ostream & operator<< (std::ostream & os, const F_Expr2<A,T> & v)
	{
		int n1;
		int n2;
		os << std::endl << "[" << std::endl;
		for (n1 = 0; n1 < v.get_dim1(); ++n1)
		{
			for (n2 = 0; n2 < v.get_dim2(); ++n2)
			{
				os << v(n1,n2) << " ";
			}
			os << std::endl;
		}
		os << "]" << std::endl;

		return os;
	}

};


#endif

