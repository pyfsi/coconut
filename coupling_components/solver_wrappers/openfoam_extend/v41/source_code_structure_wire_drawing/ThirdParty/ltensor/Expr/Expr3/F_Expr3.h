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
#ifndef F_EXPR3
#define F_EXPR3

template < class A, class T>
class F_Expr3
{
	typedef F_Expr3<A,T> IndExprA;
      private:
	A  a;


      public:
    ~F_Expr3 ()
	{


	}

	F_Expr3 (const A & arhs): a(arhs)
	{

	}


	T operator() (const int N1, const int N2, const int N3) const
	{
		return a(N1,N2,N3);
	}

	T & operator() (const int N1, const int N2, const int N3)
	{
		return a(N1,N2,N3);
	}

	int get_dim1() const
	{
		return a.get_dim1();
	}

	int get_dim2() const
	{
		return a.get_dim2();
	}

	int get_dim3() const
	{
		return a.get_dim3();
	}


	friend std::ostream & operator<< (std::ostream & os, const F_Expr3<A,T> & v)
	{
		int n1;
		int n2;
		int n3;
		os << std::endl << "[" << std::endl<< std::endl;
		for (n1 = 0; n1 < v.get_dim1(); ++n1)
		{
			os << "( "<< n1 << " , j , k )"<< std::endl;
			for (n2 = 0; n2 < v.get_dim2(); ++n2)
			{
				for (n3 = 0; n3 < v.get_dim3(); ++n3)
				{
					os << v(n1,n2,n3) << " ";
				}
				os << std::endl;
			}
			os << std::endl;
		}
		os << "]" << std::endl;
		return os;
	}


	template <class U>
    inline const F_Expr3<A,T> &
	operator*=(const U &val){
	      for(int i=0;i<this->get_dim1();i++)
            for(int j=0;j<this->get_dim2();j++)
                for(int k=0;k<this->get_dim3();k++)
                    (*this)(i,j,k)*=val;

        return *this;

	}

	// T(i,j,k) = B(i,j,k)....
	//
	template<class B,class U>
	inline const F_Expr3<A,T> &
	operator=(const F_Expr3<B,U> &rhs){
	    #ifdef CHECK_OOB
            assert(this->get_dim1()==rhs.get_dim1());
            assert(this->get_dim2()==rhs.get_dim2());
            assert(this->get_dim3()==rhs.get_dim3());
	    #endif
	    Lijk_equals_Rijk(*this,rhs);
        return (*this);
	}

	inline const F_Expr3<A,T> &
	operator=(const F_Expr3<A,T> &rhs){
        #ifdef CHECK_OOB
            assert(this->get_dim1()==rhs.get_dim1());
            assert(this->get_dim2()==rhs.get_dim2());
            assert(this->get_dim3()==rhs.get_dim3());
	    #endif
	    Lijk_equals_Rijk(*this,rhs);
        return (*this);

	}

	template <class base>
	 F_Expr3<A,T> &
	operator=(const Marray<T,3,base> &rhs){
        #ifdef CHECK_OOB
            assert(this->get_dim1()<=rhs.get_dim1());
            assert(this->get_dim2()<=rhs.get_dim2());
            assert(this->get_dim3()<=rhs.get_dim3());
	    #endif
	    Lijk_equals_Rijk(*this,rhs);
        return (*this);
	}


	template<class B,class U>
	inline const F_Expr3<A,T> &
	operator+=(const F_Expr3<B,U> &rhs){
	    #ifdef CHECK_OOB
            assert(this->get_dim1()==rhs.get_dim1());
            assert(this->get_dim2()==rhs.get_dim2());
            assert(this->get_dim3()==rhs.get_dim3());
	    #endif
	    Lijk_plusequals_Rijk(*this,rhs);
        return (*this);
	}

	inline const F_Expr3<A,T> &
	operator+=(const F_Expr3<A,T> &rhs){
        #ifdef CHECK_OOB
            assert(this->get_dim1()==rhs.get_dim1());
            assert(this->get_dim2()==rhs.get_dim2());
            assert(this->get_dim3()==rhs.get_dim3());
	    #endif
	    Lijk_plusequals_Rijk(*this,rhs);
        return (*this);

	}

	template <class base>
	 F_Expr3<A,T> &
	operator+=(const Marray<T,3,base> &term){
        #ifdef CHECK_OOB
            assert(this->get_dim1()<=term.get_dim1());
            assert(this->get_dim2()<=term.get_dim2());
            assert(this->get_dim3()<=term.get_dim3());
	    #endif
	    Lijk_plusequals_Rijk(*this,term);
        return (*this);
	}


	template <class U,class base>
	 F_Expr3<A,T> &
	operator+=(const Marray<U,3,base> &term){
        #ifdef CHECK_OOB
	    assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2() &&
               this->get_dim3()<=term.get_dim3()
              );
	    #endif
	  Lijk_plusequals_Rijk(*this,term);
        return *this;
	}



	template<class B,class U>
	inline const F_Expr3<A,T> &
	operator-=(const F_Expr3<B,U> &rhs){
	    #ifdef CHECK_OOB
            assert(this->get_dim1()==rhs.get_dim1());
            assert(this->get_dim2()==rhs.get_dim2());
            assert(this->get_dim3()==rhs.get_dim3());
	    #endif
	    Lijk_minusequals_Rijk(*this,rhs);
        return (*this);
	}

	inline const F_Expr3<A,T> &
	operator-=(const F_Expr3<A,T> &rhs){
        #ifdef CHECK_OOB
            assert(this->get_dim1()==rhs.get_dim1());
            assert(this->get_dim2()==rhs.get_dim2());
            assert(this->get_dim3()==rhs.get_dim3());
	    #endif
	    Lijk_minusequals_Rijk(*this,rhs);
        return (*this);

	}

	template <class base>
	 F_Expr3<A,T> &
	operator-=(const Marray<T,3,base> &rhs){
        #ifdef CHECK_OOB
            assert(this->get_dim1()<=rhs.get_dim1());
            assert(this->get_dim2()<=rhs.get_dim2());
            assert(this->get_dim3()<=rhs.get_dim3());
	    #endif
	    Lijk_minusequals_Rijk(*this,rhs);
        return (*this);
	}

	template <class U,class base>
	 F_Expr3<A,T> &
	operator-=(const Marray<U,3,base> &term){
        #ifdef CHECK_OOB
	    assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2() &&
               this->get_dim3()<=term.get_dim3()
              );
	    #endif
	  Lijk_minusequals_Rijk(*this,term);
        return *this;
	}



	template <class U,class base>
	 F_Expr3<A,T> &
	operator=(const Marray<U,3,base> &term){
        #ifdef CHECK_OOB
	    assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2() &&
               this->get_dim3()<=term.get_dim3()
              );
	    #endif
	  Lijk_equals_Rijk(*this,term);
        return *this;
	}

    template <class U>
	F_Expr3<A,T> &
	operator=(const U &val){
	     for(int i=0;i<this->get_dim1();i++)
            for(int j=0;j<this->get_dim2();j++)
                for(int k=0;k<this->get_dim3();k++)
                    (*this)(i,j,k)=val;
        return *this;
	}

	template <class U>
	F_Expr3<A,T> &
	operator+=(const U &val){
	     for(int i=0;i<this->get_dim1();i++)
            for(int j=0;j<this->get_dim2();j++)
                for(int k=0;k<this->get_dim3();k++)
                    (*this)(i,j,k)+=val;
        return *this;
	}

template <class U>
	F_Expr3<A,T> &
	operator-=(const U &val){
	     for(int i=0;i<this->get_dim1();i++)
            for(int j=0;j<this->get_dim2();j++)
                for(int k=0;k<this->get_dim3();k++)
                    (*this)(i,j,k)-=val;
        return *this;
	}



};



#endif
