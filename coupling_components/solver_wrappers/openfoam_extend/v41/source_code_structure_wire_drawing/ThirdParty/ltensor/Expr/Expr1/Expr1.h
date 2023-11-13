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
#ifndef Expr1_H
#define Expr1_H

//#define CHECK_Expr1

//#define USE_ASSERT_Expr1
#ifdef USE_ASSERT_Expr1
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
*	this class represents expressions with only 1 index
*	This is the only type of object that take part in indexed expressions.
*	\tparam A A type indicating what kind of object it holds.
*	\tparam T The type of the data
*	\tparam i It's the char that represents the letter of the tensor expression Eg: A(i)=B(i)
*/
template < class A, class T, char i >
class Expr1
{
	typedef Expr1<A,T,i> IndExprA;
      private:
	/*! \brief A copy to the Object it contains
	*	\remarks Although this is a copy, all the objects this can hold, have references to arrays, not arrays.
	*/
	A  a;

      public:
	/*! \brief A simple destructor*/
	~Expr1 ()
	{

		#ifdef CHECK_Expr1
		std::cout << "Expr1<Obj>.Destruct = " << (this) <<std::endl;
		#endif
	}
	/*! \brief A simple constructor
	*
	*	\param arhs The object to hold
	*/
	Expr1 (const A & arhs): a(arhs)
	{
		#ifdef CHECK_Expr1
		std::cout << "Expr1D.Constructor(obj,index), " << (*this)<< std::endl;
		#endif
	}

	/*! \brief Indexing operator*/
	T operator() (const int N) const
	{
		#ifdef USE_ASSERT_Expr1
			const int dim1 = a.get_dim1();
			assert (N >= 0 && N < dim1);
			assert (dim1 != 0);
		#endif
		return a(N);
	}

	/*! \brief Indexing operator*/
	T & operator()(const int N)
	{
		#ifdef USE_ASSERT_Expr1
			const int dim1 = a.get_dim1();
			assert (N >= 0 && N < dim1);
		      assert (dim1 != 0);
		#endif
		return a(N);
	}

	/*!\brief Returns the first dimension*/
	int get_dim1 () const
	{
		return a.get_dim1();
	}

	/*! \brief Equal operator*/
	inline  Expr1<A,T,i> &
	operator=(const Expr1<A,T,i> &rhs);

	/*! \brief Equal operator*/
	template<class B,class U>
	inline  Expr1<A,T,i> &
	operator=(const Expr1<B,U,i> &rhs);

	/*! \brief Plus Equal operator*/
	template<class B,class U>
	inline  Expr1<A,T,i> &
	operator+=(const Expr1<B,U,i> &rhs);

	/*!\brief  Minus Equal operator*/
	template<class B,class U>
	inline  Expr1<A,T,i> &
	operator-=(const Expr1<B,U,i> &rhs);



	/*! \brief Debug Function

		Shows debug info
	*/
	void show_indexeddata()
	{
		std::cout 	<< "BEGIN.Expr1<Obj> = " << this << std::endl;
		std::cout 	<< "index1 = " << i  <<std::endl;
		std::cout 	<< "Obj = " << (*this);
		std::cout 	<< "END.Expr1<Obj>"<< std::endl;
		return;
	}

	/*! Streams function overload */
	friend std::ostream & operator<< (std::ostream & os, const Expr1 & v)
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
