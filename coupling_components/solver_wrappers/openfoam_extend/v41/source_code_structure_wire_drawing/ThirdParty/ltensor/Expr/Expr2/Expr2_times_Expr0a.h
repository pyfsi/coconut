/* Multipliess a Tensor1 by a generic (or vice versa), yielding a Tensor1.
   Usually used for doubles, but could be used for complex, etc.  All
   that it requires is that you can add an element of the Tensor1 to
   it.  */

#ifndef Expr2_times_Expr0_H
#define Expr2_times_Expr0_H

//#define CHECK_Expr2_times_Expr0

/* A(i) * d0 -> Tensor1 */

template < class A, class T, class U, char i > class Expr2_times_Expr0
{
	const Expr2 < A, T, i > iterA;
	const U d;
      public:
	const int dim1;

	Expr2_times_Expr0 (const Expr2 < A, T, i > &a, const U & d0):
	iterA(a), d(d0),dim1(a.get_dim1())
	{
	}

	typename promote < T, U >::V operator () (const int N) const
	{
		return iterA (N) * d;
	}

	int get_dim1 () const
	{
		return dim1;
	}

};

template < class A, class T, class U, char i >
inline const Expr2 < const Expr2_times_Expr0 < A, T, U, i >, 
				typename promote < T, U >::V, i >
operator* (const Expr2 < A, T, i > &a, const U & d0)
{
	typedef const Expr2_times_Expr0 < A, T, U, i > ExprObj;
	return Expr2 < ExprObj, typename promote < T, U >::V,i > (ExprObj (a, d0));
}

/* d0 * A(i) -> Tensor1 */

template < class A, class T, class U, char i >
inline const Expr2 < const Expr2_times_Expr0 < A, T, U, i >,
				typename promote < T, U >::V, i >
operator* (const U & d0, const Expr2 < A, T, i > &a)
{
	typedef const Expr2_times_Expr0 < A, T, U, i > ExprObj;
	return Expr2 < ExprObj, typename promote < T, U >::V,i > (ExprObj (a, d0));
}


#endif
