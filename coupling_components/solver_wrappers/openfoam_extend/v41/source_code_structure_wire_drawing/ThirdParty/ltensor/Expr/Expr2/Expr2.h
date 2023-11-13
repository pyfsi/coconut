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
#ifndef Expr2_H
#define Expr2_H

//#define CHECK_Expr2

//#define USE_ASSERT_Expr2
#ifdef USE_ASSERT_Expr2
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
*	this class represents expressions with 2 indexes
*	This is the only type of object that take part in indexed expressions. 
*	\tparam A A type indicating what kind of object it holds. 
*	\tparam T The type of the data
*	\tparam i It's the char that represents the letter of the first index of the tensor expression Eg: A(i,j)=B(i,j)
*	\tparam j It's the char that represents the letter of the second index of the tensor expression Eg: A(i,j)=B(i,j)
*/


template < class A, class T, char i , char j>
class Expr2
{
	typedef Expr2<A,T,i,j> IndExprA;

	private:
/*! \brief A copy to the Object it contains
	*	\remarks Although this is a copy, all the objects this can hold, have references to arrays, not arrays.
	*/
	A  a;

      public:
/*! \brief A simple destructor*/
      ~Expr2 ()
	{
		#ifdef CHECK_Expr2
		std::cout << "Expr2<Obj>.Destruct = " << (this) <<std::endl;
		#endif

	}
/*! \brief A simple constructor
	*
	*	\param arhs The object to hold
	*/
	Expr2 (const A & arhs): a(arhs)
	{
		#ifdef CHECK_Expr2
		std::cout << "Expr2D.Constructor(obj,index), " << (*this)<< std::endl;
		#endif
	}

/*! \brief Indexing operator*/
	T operator() (const int N1, const int N2) const
	{
		#ifdef USE_ASSERT_Expr2
			const int dim1 = a.get_dim1();
			const int dim2 = a.get_dim2();
			assert (N1 >= 0 && N1< dim1);
			assert (dim1 != 0);
			assert (N2 >= 0 && N2< dim2);
			assert (dim2 != 0);
		#endif
		return a(N1,N2);
	}
/*! \brief Indexing operator*/
	T & operator()(const int N1, const int N2)
	{
		#ifdef USE_ASSERT_Expr2
			const int dim1 = a.get_dim1();
			const int dim2 = a.get_dim2();
			assert (N1 >= 0 && N1< dim1);
			assert (dim1 != 0);
			assert (N2 >= 0 && N2< dim2);
			assert (dim2 != 0);
		#endif
		return a(N1,N2);
	}
/*!\brief Returns the first dimension*/
	int get_dim1 () const
	{
		return a.get_dim1();
	}

	/*!\brief Returns the second dimension*/
	int get_dim2 () const
	{
		return a.get_dim2();
	}

	/*! \brief Equal operator 
	*
	*	T(i,j) = B(i,j) 
	*/
	inline const Expr2<A,T,i,j> &
	operator=(const Expr2<A,T,i,j> &rhs);
/*! \brief Equal operator 
	*
	*	T(i,j) = B(i,j) 
	*/
	template<class B,class U>
	inline const Expr2<A,T,i,j> &
	operator=(const Expr2<B,U,i,j> &rhs);
/*! \brief Plus Equal operator 
	*
	*	T(i,j) += B(i,j) 
	*/
	template<class B,class U>
	inline const Expr2<A,T,i,j> &
	operator+=(const Expr2<B,U,i,j> &rhs);
/*! \brief Minus Equal operator 
	*
	*	T(i,j) -= B(i,j) 
	*/
	template<class B,class U>
	inline const Expr2<A,T,i,j> &
	operator-=(const Expr2<B,U,i,j> &rhs);



	/*! \brief Equal operator 
	*
	*	T(i,j) = B(j,i) 
	*	\remarks NOTE the permutted indexes
	*/
	inline const Expr2<A,T,i,j> &
	operator=(const Expr2<A,T,j,i> &rhs);

/*! \brief Equal operator 
	*
	*	T(i,j) = B(j,i) 
	*	\remarks NOTE the permutted indexes
	*/

	template<class B,class U>
	inline const Expr2<A,T,i,j> &
	operator=(const Expr2<B,U,j,i> &rhs);


	/*! \brief Plus Equal operator 
	*
	*	T(i,j) += B(j,i) 
	*	\remarks NOTE the permutted indexes
	*/
	template<class B,class U>
	inline const Expr2<A,T,i,j> &
	operator+=(const Expr2<B,U,j,i> &rhs);

	/*! \brief Minus Equal operator 
	*
	*	T(i,j) -= B(j,i) 
	*	\remarks NOTE the permutted indexes
	*/
	template<class B,class U>
	inline const Expr2<A,T,i,j> &
	operator-=(const Expr2<B,U,j,i> &rhs);


/*! \brief Debug Function
	
		Shows debug info
	*/
	void show_indexeddata()
	{
		std::cout 	<< "BEGIN.Expr2<Obj> = " << this << std::endl;
		std::cout 	<< "index1 = " << i << ", index2 = " << j << std::endl;
		std::cout 	<< "Obj = " << (*this);
		std::cout 	<< "END.Expr2<Obj>"<< std::endl;
		return;
	}
/*! Streams function overload */
	friend std::ostream & operator<< (std::ostream & os, const Expr2 & v)
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
