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
#ifndef Expr4_H
#define Expr4_H


#ifdef USE_ASSERT_Expr4
#include <assert.h>
#endif

/*! \brief A class that represents an expression
*
*	This class acts as a container of any indexable object. It usually contains Encapsulates objects
*	and operator objects. Every object that participates in indexed expressions must be inside an Expr Object.
*	This way trees of Expr Objects are built when complex operations are required. To perform this, each object
*	calls the operator() of its contained object, which applies modifications and returns the data. This call is
*	done recursively through the tree.
*	This objects is the one that contains the chars that represents the indexes of the tensor expression. In this case
*	this class represents expressions with 4 indexes
*	This is the only type of object that take part in indexed expressions.
*	\tparam A A type indicating what kind of object it holds.
*	\tparam T The type of the data
*	\tparam i It's the char that represents the letter of the first index of the tensor expression Eg: A(i,j)=B(i,j)
*	\tparam j It's the char that represents the letter of the second index of the tensor expression Eg: A(i,j)=B(i,j)
*	\tparam k It's the char that represents the letter of the third index of the tensor expression Eg: A(i,j)=B(i,j)
*/

template < class A, class T, char i, char j , char k, char l>
class Expr4
{
	typedef Expr4<A,T,i,j,k,l> IndExprA;
      private:

	/*! \brief A copy to the Object it contains
	*	\remarks Although this is a copy, all the objects this can hold, have references to arrays, not arrays.
	*/
	A  a;
      public:
/*! \brief A simple destructor*/
      ~Expr4 ()
	{
		#ifdef CHECK_Expr4
		std::cout << "Expr4<Obj>.Destruct = " << (this) <<std::endl;
		#endif
	}
/*! \brief A simple constructor
	*
	*	\param arhs The object to hold
	*/
	Expr4 (const A & arhs): a(arhs)
	{
		#ifdef CHECK_Expr4
		std::cout << "Expr4D.Constructor(obj,index), " << (*this)<< std::endl;
		#endif
	}
/*! \brief Indexing operator*/

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return a(N1,N2,N3,N4);
	}
/*! \brief Indexing operator*/
	T & operator() (const int N1, const int N2, const int N3, const int N4)
	{
		return a(N1,N2,N3,N4);
	}
/*!\brief Returns the first dimension*/
	int get_dim1() const
	{
		return a.get_dim1();
	}
/*!\brief Returns the second dimension*/
	int get_dim2() const
	{
		return a.get_dim2();
	}
/*!\brief Returns the third dimension*/
	int get_dim3() const
	{
		return a.get_dim3();
	}
/*!\brief Returns the fourth dimension*/
	int get_dim4() const
	{
		return a.get_dim4();
	}
/*! \brief Debug Function*/
	void show_indexeddata()
	{
		std::cout 	<< "BEGIN.Expr4<Obj> = " << this << std::endl;
		std::cout 	<< "index1 = " << i << ", index2 = " << j
				<< ", index3 = " << k <<  ", index4 = " << l <<std::endl;
		std::cout 	<< "Obj = " << (*this);
		std::cout 	<< "END.Expr4<Obj>"<< std::endl;
		return;
	}
/*! Streams function overload */
	friend std::ostream & operator<< (std::ostream & os, const Expr4 & v)
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


	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(i,j,k,l)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,i,j,k,l> &rhs);

	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(i,j,k,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,i,j,k,l> &rhs);
/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(i,j,k,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,i,j,k,l> &rhs);
/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(i,j,k,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,i,j,k,l> &rhs);



	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(j,i,k,l)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,j,i,k,l> &rhs);
/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(j,i,k,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,j,i,k,l> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(j,i,k,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,j,i,k,l> &rhs);
/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(j,i,k,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,j,i,k,l> &rhs);


	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(i,j,l,k)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,i,j,l,k> &rhs);
/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(i,j,l,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,i,j,l,k> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(i,j,l,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,i,j,l,k> &rhs);
/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(i,j,l,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,i,j,l,k> &rhs);

	// L(i,j,k,l) = R(j,i,l,k)....04

	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(j,i,l,k)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,j,i,l,k> &rhs);
/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(j,i,l,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,j,i,l,k> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(j,i,l,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,j,i,l,k> &rhs);

	/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(j,i,l,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,j,i,l,k> &rhs);


	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(k,l,i,j)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,k,l,i,j> &rhs);

	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(k,l,i,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,k,l,i,j> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(k,l,i,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,k,l,i,j> &rhs);

	/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(k,l,i,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,k,l,i,j> &rhs);


	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(l,k,i,j)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,l,k,i,j> &rhs);
/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(l,k,i,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,l,k,i,j> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(l,k,i,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,l,k,i,j> &rhs);

	/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(l,k,i,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,l,k,i,j> &rhs);


	// L(i,j,k,l) = R(l,k,j,i)....07

	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(l,k,j,i)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,l,k,j,i> &rhs);
/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(l,k,j,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,l,k,j,i> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(l,k,j,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,l,k,j,i> &rhs);

	/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(l,k,j,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,l,k,j,i> &rhs);

	// L(i,j,k,l) = R(k,l,j,i)....08

	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(k,l,j,i)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,k,l,j,i> &rhs);
/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(k,l,j,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,k,l,j,i> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(k,l,j,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,k,l,j,i> &rhs);

	/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(k,l,j,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,k,l,j,i> &rhs);



	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(l,j,k,i)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,l,j,k,i> &rhs);
/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(l,j,k,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,l,j,k,i> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(l,j,k,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,l,j,k,i> &rhs);
/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(l,j,k,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,l,j,k,i> &rhs);



	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(i,k,j,l)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,i,k,j,l> &rhs);


	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(i,k,j,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,i,k,j,l> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(i,k,j,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,i,k,j,l> &rhs);

	/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(i,k,j,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,i,k,j,l> &rhs);


	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(i,k,l,j)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,i,k,l,j> &rhs);

	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(i,k,l,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,i,k,l,j> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(i,k,l,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,i,k,l,j> &rhs);

	/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(i,k,l,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,i,k,l,j> &rhs);


	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(i,l,j,k)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,i,l,j,k> &rhs);
/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(i,l,j,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,i,l,j,k> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(i,l,j,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,i,l,j,k> &rhs);

	/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(i,l,j,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,i,l,j,k> &rhs);


	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(i,l,k,j)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,i,l,k,j> &rhs);
/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(i,l,k,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,i,l,k,j> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(i,l,k,j)
	*/

	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,i,l,k,j> &rhs);

	/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(i,l,k,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,i,l,k,j> &rhs);


	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(j,k,i,l)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,j,k,i,l> &rhs);
	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(j,k,i,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,j,k,i,l> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(j,k,i,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,j,k,i,l> &rhs);

	/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(j,k,i,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,j,k,i,l> &rhs);




	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(j,k,l,i)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,j,k,l,i> &rhs);
/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(j,k,l,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,j,k,l,i> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(j,k,l,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,j,k,l,i> &rhs);

	/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(j,k,l,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,j,k,l,i> &rhs);




	/*! \brief Equal operator
	*
	*	 L(i,j,k,l) = R(j,l,i,k)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,j,l,i,k> &rhs);
/*! \brief Equal operator
	*
	*	 L(i,j,k,l) = R(j,l,i,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,j,l,i,k> &rhs);

	/*! \brief Plus Equal operator
	*
	*	 L(i,j,k,l) += R(j,l,i,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,j,l,i,k> &rhs);
/*! \brief Minus Equal operator
	*
	*	 L(i,j,k,l) -= R(j,l,i,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,j,l,i,k> &rhs);


	// L(i,j,k,l) = R(j,l,k,i)....18

	/*! \brief Equal operator
	*
	*	 L(i,j,k,l) = R(j,l,k,i)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,j,l,k,i> &rhs);
/*! \brief Equal operator
	*
	*	 L(i,j,k,l) = R(j,l,k,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,j,l,k,i> &rhs);
/*! \brief Plus Equal operator
	*
	*	 L(i,j,k,l) += R(j,l,k,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,j,l,k,i> &rhs);

	/*! \brief Minus Equal operator
	*
	*	 L(i,j,k,l) -= R(j,l,k,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,j,l,k,i> &rhs);


	/*! \brief Equal operator
	*
	*	 L(i,j,k,l) = R(k,i,j,l)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,k,i,j,l> &rhs);
	/*! \brief Equal operator
	*
	*	 L(i,j,k,l) = R(k,i,j,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,k,i,j,l> &rhs);
/*! \brief Plus Equal operator
	*
	*	 L(i,j,k,l) += R(k,i,j,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,k,i,j,l> &rhs);
/*! \brief Minus Equal operator
	*
	*	 L(i,j,k,l) -= R(k,i,j,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,k,i,j,l> &rhs);


	/*! \brief Equal operator
	*
	*	 L(i,j,k,l) = R(k,i,l,j)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,k,i,l,j> &rhs);
/*! \brief Equal operator
	*
	*	 L(i,j,k,l) = R(k,i,l,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,k,i,l,j> &rhs);
	/*! \brief Plus Equal operator
	*
	*	 L(i,j,k,l) += R(k,i,l,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,k,i,l,j> &rhs);
	/*! \brief Minus  Equal operator
	*
	*	 L(i,j,k,l) -= R(k,i,l,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,k,i,l,j> &rhs);




	/*! \brief Equal operator
	*
	*	 L(i,j,k,l) = R(k,j,i,l)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,k,j,i,l> &rhs);
	/*! \brief Equal operator
	*
	*	 L(i,j,k,l) = R(k,j,i,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,k,j,i,l> &rhs);

	/*! \brief Plus Equal operator
	*
	*	 L(i,j,k,l) += R(k,j,i,l)
	*/

	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,k,j,i,l> &rhs);

	/*! \brief Minus Equal operator
	*
	*	 L(i,j,k,l) -= R(k,j,i,l)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,k,j,i,l> &rhs);



	/*! \brief Equal operator
	*
	*	 L(i,j,k,l) = R(k,j,l,i)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,k,j,l,i> &rhs);

	/*! \brief Equal operator
	*
	*	 L(i,j,k,l) = R(k,j,l,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,k,j,l,i> &rhs);

	/*! \brief Plus Equal operator
	*
	*	 L(i,j,k,l) += R(k,j,l,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,k,j,l,i> &rhs);


	/*! \brief Minus Equal operator
	*
	*	 L(i,j,k,l) -= R(k,j,l,i)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,k,j,l,i> &rhs);


	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(l,i,j,k)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,l,i,j,k> &rhs);
/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(l,i,j,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,l,i,j,k> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(l,i,j,k)
	*/

	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,l,i,j,k> &rhs);

	/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(l,i,j,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,l,i,j,k> &rhs);



	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(l,i,k,j)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,l,i,k,j> &rhs);
	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(l,i,k,j)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,l,i,k,j> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) = R(l,i,k,j)
	*/

	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,l,i,k,j> &rhs);

	/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) = R(l,i,k,j)
	*/

	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,l,i,k,j> &rhs);




	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(l,j,i,k)
	*/
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<A,T,l,j,i,k> &rhs);
	/*! \brief Equal operator
	*
	*	L(i,j,k,l) = R(l,j,i,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator=(const Expr4<B,U,l,j,i,k> &rhs);

	/*! \brief Plus Equal operator
	*
	*	L(i,j,k,l) += R(l,j,i,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator+=(const Expr4<B,U,l,j,i,k> &rhs);

	/*! \brief Minus Equal operator
	*
	*	L(i,j,k,l) -= R(l,j,i,k)
	*/
	template<class B, class U>
	inline const Expr4<A,T,i,j,k,l> &
	operator-=(const Expr4<B,U,l,j,i,k> &rhs);

};




#endif
