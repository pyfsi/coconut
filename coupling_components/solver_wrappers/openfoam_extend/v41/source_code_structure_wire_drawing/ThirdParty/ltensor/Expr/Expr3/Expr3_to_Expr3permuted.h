/////////////////////////////////////////////////////
// Given Expr2=T(i,j) define a new permuted Expr2 P object
// so to define: T(i,j)=P(j,i)
/////////////////////////////////////////////////////

#ifndef Expr3_to_Expr3permuted_H
#define Expr3_to_Expr3permuted_H

////////////////////////////////////////////
// Define Object B(j,i)=A(i,j)  
///////////////////////////////////////////

template < class A, class T, char i , char j> 
class Expr3ji_to_Expr3ij
{
	const Expr2 < A, T, j, i > TA;
      public:

	Expr2ji_to_Expr2ij (const Expr2 < A, T, j, i > &a):
	TA(a)
	{
	}

	int get_dim1 () const
	{
		return TA.get_dim2();
	}
	
	int get_dim2 () const
	{
		return TA.get_dim1();
	}
	
	T operator () (const int N1,const int N2) const
	{
		return TA(N2,N1);
	}
};

#endif
