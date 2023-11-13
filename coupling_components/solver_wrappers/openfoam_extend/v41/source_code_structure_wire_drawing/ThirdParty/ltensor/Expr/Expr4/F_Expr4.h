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
#ifndef F_EXPR4
#define F_EXPR4

template < class A, class T>
class F_Expr4
{
	typedef F_Expr4<A,T> IndExprA;
      private:
	A  a;

      public:

      ~F_Expr4 ()
	{
		#ifdef CHECK_Expr4
		std::cout << "Expr4<Obj>.Destruct = " << (this) <<std::endl;
		#endif
    }

	F_Expr4 (const A & arhs): a(arhs)
	{

	}


	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return a(N1,N2,N3,N4);
	}

	T & operator() (const int N1, const int N2, const int N3, const int N4)
	{
		return a(N1,N2,N3,N4);
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

	int get_dim4() const
	{
		return a.get_dim4();
	}


	friend std::ostream & operator<< (std::ostream & os, const F_Expr4<A,T> & v)
	{
		int n1;
		int n2;
		int n3;
		int n4;
		os << std::endl << "[" << std::endl<< std::endl;
		for (n1 = 0; n1 < v.get_dim1(); ++n1)
		{
			for (n2 = 0; n2 < v.get_dim2(); ++n2)
			{
				os << "( "<< n1 << " , " << n2 << " , j , k )"<< std::endl;
				for (n3 = 0; n3 < v.get_dim3(); ++n3)
				{
					for (n4 = 0; n4 < v.get_dim4(); ++n4)
					{
						os << v(n1,n2,n3,n4) << " ";
					}
					os << std::endl;
				}
				os << std::endl;
			}
		}
		os << "]" << std::endl;

		return os;

	}



	// L(i,j,k,l) = R(i,j,k,l)....
	//
	inline const F_Expr4<A,T> &
	operator=(const F_Expr4<A,T> &rhs){
	    #ifdef CHECK_OOB
	   assert(this->get_dim1()==rhs.get_dim1() &&
               this->get_dim2()==rhs.get_dim2() &&
               this->get_dim3()==rhs.get_dim3() &&
               this->get_dim4()==rhs.get_dim4()
                );
	    #endif
	    Lijkl_equals_Rijkl(*this,rhs);

        return *this;


	}

	inline const F_Expr4<A,T> &
	operator+=(const F_Expr4<A,T> &rhs){
	    #ifdef CHECK_OOB
	   assert(this->get_dim1()==rhs.get_dim1() &&
               this->get_dim2()==rhs.get_dim2() &&
               this->get_dim3()==rhs.get_dim3() &&
               this->get_dim4()==rhs.get_dim4()
                );
	    #endif
	    Lijkl_plusequals_Rijkl(*this,rhs);

        return *this;


	}

	inline const F_Expr4<A,T> &
	operator-=(const F_Expr4<A,T> &rhs){
	    #ifdef CHECK_OOB
	   assert(this->get_dim1()==rhs.get_dim1() &&
               this->get_dim2()==rhs.get_dim2() &&
               this->get_dim3()==rhs.get_dim3() &&
               this->get_dim4()==rhs.get_dim4()
                );
	    #endif
	    Lijkl_minusequals_Rijkl(*this,rhs);

        return *this;


	}

	template<class B, class U>
	inline const F_Expr4<A,T> &
	operator=(const F_Expr4<B,U> &rhs){
	    #ifdef CHECK_OOB
	   assert(this->get_dim1()==rhs.get_dim1() &&
               this->get_dim2()==rhs.get_dim2() &&
               this->get_dim3()==rhs.get_dim3() &&
               this->get_dim4()==rhs.get_dim4()
                );
	    #endif
	    Lijkl_equals_Rijkl(*this,rhs);

        return *this;
    }


	template<class B, class U>
	inline const F_Expr4<A,T> &
	operator+=(const F_Expr4<B,U> &rhs){
	    #ifdef CHECK_OOB
	   assert(this->get_dim1()==rhs.get_dim1() &&
               this->get_dim2()==rhs.get_dim2() &&
               this->get_dim3()==rhs.get_dim3() &&
               this->get_dim4()==rhs.get_dim4()
                );
	    #endif
	    Lijkl_plusequals_Rijkl(*this,rhs);

        return *this;
    }


	template<class B, class U>
	inline const F_Expr4<A,T> &
	operator-=(const F_Expr4<B,U> &rhs){
	    #ifdef CHECK_OOB
	   assert(this->get_dim1()==rhs.get_dim1() &&
               this->get_dim2()==rhs.get_dim2() &&
               this->get_dim3()==rhs.get_dim3() &&
               this->get_dim4()==rhs.get_dim4()
                );
	    #endif
	    Lijkl_minusequals_Rijkl(*this,rhs);

        return *this;
    }


    template <class U>
    inline const F_Expr4<A,T> &
	operator=(const U &val){
	      for(int i=0;i<this->get_dim1();i++)
            for(int j=0;j<this->get_dim2();j++)
                for(int k=0;k<this->get_dim3();k++)
                    for(int l=0;l<this->get_dim4();l++)
                        (*this)(i,j,k,l)=val;

        return *this;


	}


	template <class base>
	 F_Expr4<A,T> &
	operator=(const Marray<T,4,base> &term){
	     #ifdef CHECK_OOB
        assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2() &&
               this->get_dim3()<=term.get_dim3() &&
               this->get_dim4()<=term.get_dim4()
                );
	    #endif
	     Lijkl_equals_Rijkl(*this,term);

        return *this;
	}


	template <class U>
    inline const F_Expr4<A,T> &
	operator*=(const U &val){
	      for(int i=0;i<this->get_dim1();i++)
            for(int j=0;j<this->get_dim2();j++)
                for(int k=0;k<this->get_dim3();k++)
                    for(int l=0;l<this->get_dim4();l++)
                        (*this)(i,j,k,l)*=val;

        return *this;

	}


	template <class base>
	 F_Expr4<A,T> &
	operator+=(const Marray<T,4,base> &term){
	     #ifdef CHECK_OOB
        assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2() &&
               this->get_dim3()<=term.get_dim3() &&
               this->get_dim4()<=term.get_dim4()
                );
	    #endif
	     Lijkl_plusequals_Rijkl(*this,term);

        return *this;
	}

	template <class base>
	 F_Expr4<A,T> &
	operator-=(const Marray<T,4,base> &term){
	     #ifdef CHECK_OOB
        assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2() &&
               this->get_dim3()<=term.get_dim3() &&
               this->get_dim4()<=term.get_dim4()
                );
	    #endif
	     Lijkl_minusequals_Rijkl(*this,term);

        return *this;
	}



	 template <class U,class base>
	 F_Expr4<A,T> &
	 operator=(const Marray<U,4,base> &term){
	     #ifdef CHECK_OOB
	   assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2() &&
               this->get_dim3()<=term.get_dim3() &&
               this->get_dim4()<=term.get_dim4()
                );
	    #endif
	     Lijkl_equals_Rijkl(*this,term);

        return *this;
	}


	template <class U,class base>
	 F_Expr4<A,T> &
	 operator+=(const Marray<U,4,base> &term){
	     #ifdef CHECK_OOB
	   assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2() &&
               this->get_dim3()<=term.get_dim3() &&
               this->get_dim4()<=term.get_dim4()
                );
	    #endif
	     Lijkl_plusequals_Rijkl(*this,term);

        return *this;
	}

	template <class U,class base>
	 F_Expr4<A,T> &
	 operator-=(const Marray<U,4,base> &term){
	     #ifdef CHECK_OOB
	   assert(this->get_dim1()<=term.get_dim1() &&
               this->get_dim2()<=term.get_dim2() &&
               this->get_dim3()<=term.get_dim3() &&
               this->get_dim4()<=term.get_dim4()
                );
	    #endif
	     Lijkl_minusequals_Rijkl(*this,term);

        return *this;
	}

};



#endif
