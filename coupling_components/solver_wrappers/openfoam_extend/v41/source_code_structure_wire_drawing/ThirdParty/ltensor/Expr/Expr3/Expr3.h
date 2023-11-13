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
#ifndef Expr3_H
#define Expr3_H


#ifdef USE_ASSERT_Expr3
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
*	this class represents expressions with 3 indexes
*	This is the only type of object that take part in indexed expressions. 
*	\tparam A A type indicating what kind of object it holds. 
*	\tparam T The type of the data
*	\tparam i It's the char that represents the letter of the first index of the tensor expression Eg: A(i,j)=B(i,j)
*	\tparam j It's the char that represents the letter of the second index of the tensor expression Eg: A(i,j)=B(i,j)
*	\tparam k It's the char that represents the letter of the third index of the tensor expression Eg: A(i,j)=B(i,j)
*/

template < class A, class T, char i, char j , char k>
class Expr3
{
	typedef Expr3<A,T,i,j,k> IndExprA;
      private:

		  /*! \brief A copy to the Object it contains
	*	\remarks Although this is a copy, all the objects this can hold, have references to arrays, not arrays.
	*/
	A  a;
      public:
/*! \brief A simple destructor*/
      ~Expr3 ()
	{
		#ifdef CHECK_Expr3
		std::cout << "Expr3<Obj>.Destruct = " << (this) <<std::endl;
		#endif
	}
/*! \brief A simple constructor
	*
	*	\param arhs The object to hold
	*/
	Expr3 (const A & arhs): a(arhs)
	{
		#ifdef CHECK_Expr3
		std::cout << "Expr3D.Constructor(obj,index), " << (*this)<< std::endl;
		#endif
	}
/*! \brief Indexing operator*/

	T operator() (const int N1, const int N2, const int N3) const
	{
		return a(N1,N2,N3);
	}
/*! \brief Indexing operator*/
	T & operator() (const int N1, const int N2, const int N3)
	{
		return a(N1,N2,N3);
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
/*! \brief Debug Function
	
		Shows debug info
	*/
	void show_indexeddata()
	{
		std::cout 	<< "BEGIN.Expr3<Obj> = " << this << std::endl;
		std::cout 	<< "index1 = " << i << ", index2 = " << j
				<< ", index3 = " << k << std::endl;
		std::cout 	<< "Obj = " << (*this);
		std::cout 	<< "END.Expr3<Obj>"<< std::endl;
		return;
	}
/*! Streams function overload */
	friend std::ostream & operator<< (std::ostream & os, const Expr3 & v)
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

	/*! \brief Equal operator 
	*
	*	T(i,j,k) = B(i,j,k)
	*/
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator=(const Expr3<B,T,i,j,k> &rhs);

	/*! \brief Plus Equal operator 
	*
	*	T(i,j,k) += B(i,j,k)
	*/
	
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator+=(const Expr3<B,T,i,j,k> &rhs);

	/*! \brief Minus Equal operator 
	*
	*	T(i,j,k) -= B(i,j,k)
	*/
	
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator-=(const Expr3<B,T,i,j,k> &rhs);

	inline const Expr3<A,T,i,j,k> &
	operator=(const Expr3<A,T,i,j,k> &rhs);

	/*! \brief Equal operator 
	*
	*	T(i,j,k) = B(i,k,j)
	*	\remarks NOTE the permutted indexes
	*/
	
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator=(const Expr3<B,T,i,k,j> &rhs);

	/*! \brief Equal operator 
	*
	*	T(i,j,k) = B(i,k,j)
	*	\remarks NOTE the permutted indexes
	*/
	
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator+=(const Expr3<B,T,i,k,j> &rhs);

	/*! \brief Minus Equal operator 
	*
	*	T(i,j,k) -= B(i,k,j)
	*	\remarks NOTE the permutted indexes
	*/
	
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator-=(const Expr3<B,T,i,k,j> &rhs);


	/*! \brief Equal operator 
	*
	*	T(i,j,k) = B(i,k,j)
	*	\remarks NOTE the permutted indexes
	*/
	
	inline const Expr3<A,T,i,j,k> &
	operator=(const Expr3<A,T,i,k,j> &rhs);

	
	/*! \brief Equal operator 
	*
	*	T(i,j,k) = B(j,i,k)
	*	\remarks NOTE the permutted indexes
	*/
	
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator=(const Expr3<B,T,j,i,k> &rhs);

	/*! \brief Plus Equal operator 
	*
	*	T(i,j,k) += B(j,i,k)
	*	\remarks NOTE the permutted indexes
	*/
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator+=(const Expr3<B,T,j,i,k> &rhs);

	/*! \brief Minus Equal operator 
	*
	*	T(i,j,k) -= B(j,i,k)
	*	\remarks NOTE the permutted indexes
	*/
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator-=(const Expr3<B,T,j,i,k> &rhs);

	/*! \brief Equal operator 
	*
	*	T(i,j,k) = B(j,i,k)
	*	\remarks NOTE the permutted indexes
	*/
	inline const Expr3<A,T,i,j,k> &
	operator=(const Expr3<A,T,j,i,k> &rhs);

	
	/*! \brief Equal operator 
	*
	*	T(i,j,k) = B(k,i,j)
	*	\remarks NOTE the permutted indexes
	*/
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator=(const Expr3<B,T,k,i,j> &rhs);

	/*! \brief Plus Equal operator 
	*
	*	T(i,j,k) += B(k,i,j)
	*	\remarks NOTE the permutted indexes
	*/
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator+=(const Expr3<B,T,k,i,j> &rhs);

	/*! \brief Minus Equal operator 
	*
	*	T(i,j,k) -= B(k,i,j)
	*	\remarks NOTE the permutted indexes
	*/
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator-=(const Expr3<B,T,k,i,j> &rhs);

	/*! \brief Equal operator 
	*
	*	T(i,j,k) = B(k,i,j)
	*	\remarks NOTE the permutted indexes
	*/
	inline const Expr3<A,T,i,j,k> &
	operator=(const Expr3<A,T,k,i,j> &rhs);

	/*! \brief Equal operator 
	*
	*	T(i,j,k) = B(j,k,i)
	*	\remarks NOTE the permutted indexes
	*/
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator=(const Expr3<B,T,j,k,i> &rhs);

	/*! \brief Plus Equal operator 
	*
	*	T(i,j,k) += B(j,k,i)
	*	\remarks NOTE the permutted indexes
	*/
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator+=(const Expr3<B,T,j,k,i> &rhs);

	/*! \brief Minus Equal operator 
	*
	*	T(i,j,k) -= B(j,k,i)
	*	\remarks NOTE the permutted indexes
	*/
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator-=(const Expr3<B,T,j,k,i> &rhs);

	/*! \brief Equal operator 
	*
	*	T(i,j,k) = B(j,k,i)
	*	\remarks NOTE the permutted indexes
	*/
	inline const Expr3<A,T,i,j,k> &
	operator=(const Expr3<A,T,j,k,i> &rhs);

	
	/*! \brief Equal operator 
	*
	*	T(i,j,k) = B(k,j,i)
	*	\remarks NOTE the permutted indexes
	*/
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator=(const Expr3<B,T,k,j,i> &rhs);

	/*! \brief Plus Equal operator 
	*
	*	T(i,j,k) += B(k,j,i)
	*	\remarks NOTE the permutted indexes
	*/
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator+=(const Expr3<B,T,k,j,i> &rhs);

	/*! \brief Minus Equal operator 
	*
	*	T(i,j,k) -= B(k,j,i)
	*	\remarks NOTE the permutted indexes
	*/
	template<class B>
	inline const Expr3<A,T,i,j,k> &
	operator-=(const Expr3<B,T,k,j,i> &rhs);

	/*! \brief Equal operator 
	*
	*	T(i,j,k) = B(k,j,i)
	*	\remarks NOTE the permutted indexes
	*/
	inline const Expr3<A,T,i,j,k> &
	operator=(const Expr3<A,T,k,j,i> &rhs);


};




#endif

