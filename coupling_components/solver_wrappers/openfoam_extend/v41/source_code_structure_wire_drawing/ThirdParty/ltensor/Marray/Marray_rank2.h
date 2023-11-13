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
// -*- c++ -*-
#ifndef Marray_rank2_H
#define Marray_rank2_H

//#define USE_ASSERT_Marray
#ifdef USE_ASSERT_Marray
#include <assert.h>
#endif

#define CHECK										\
	assert( (n1<get_dim1())&& (n2<get_dim2()) );	\
	assert( (n1>=0) && (n2>=0));


//#define CHECK_Marray_rank2


#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "../algorithms/algorithms.h"

#include "../Tensor_Operations/Tensor_Operations.h"
#include "../Expr/Expr1/Expr1s.h"
#include "../Expr/Expr2/Expr2s.h"


template < class A, class T, char i > class Expr1;
template < class A, class T, char i, char j > class Expr2;
template < class A, class T> class F_Expr2;



#define rank 2

//
// FUNTIONS DECLARATIONS in order to FRIEND Function "ostream <<"  to work
//
template < class T,class base> class Marray <T,rank,base>;

template<class T,class base>
std::ostream & operator << (std::ostream & os, const Marray<T, rank,base> & v);

//
// Class Definition
//
template < class T,class base>
class Marray <T,rank,base> :public base
{


	////////////////
	//Constructors
	///////////////
public:



	Marray (long dimension1,long dimension2):base(dimension1,dimension2){
		(*this)=(T)0.0;

	}

	Marray (long dimension1,long dimension2, T valor):base(dimension1,dimension2){
		*this=valor;
	}

	Marray(){
	};



	template <class fType,class b>
	Marray<T,rank,base>& operator=(const Marray<fType,rank,b> &rterm){
		//	   //make sure sizes are compatible
		//	   //base::equals(*dynamic_cast<base*>(const_cast<Marray<T,rank,base>*> (&rterm)) );
		//// 	    #ifdef CHECK_OOB
		//// 	    assert( (get_dim1() !=0 ) && (get_dim2()!=0));
		//// 	    #endif
		//	resize(rterm.get_dim1(),rterm.get_dim2());
		//		for(int i=0;i<get_dim1();i++)
		//			for(int j=0;j<get_dim2();j++)
		//			(*this)(i,j)=rterm(i,j);
		base::operator=(*dynamic_cast<const b*>(&rterm));
		return *this;
	}

	Marray<T,rank,base>& operator=(const Marray<T,rank,base> &rterm){
		//make sure sizes are compatible
		//base::equals(*dynamic_cast<base*>(const_cast<Marray<T,rank,base>*> (&rterm)) );
		//resize(rterm.get_dim1(),rterm.get_dim2());
		//for(int i=0;i<get_dim1();i++)
		//	for(int j=0;j<get_dim2();j++)
		//	(*this)(i,j)=rterm(i,j);
		base::operator=(*dynamic_cast<const base*>(&rterm));
		return *this;

	}

	template <class A, class U>
	Marray<T,rank,base>& operator=(const F_Expr2 < A, U> &expr){

		for(int i=0;i<expr.get_dim1();i++)
			for(int j=0;j<expr.get_dim1();j++)
				(*this)(i,j)=expr(i,j);


		return (*this);
	}

	template <class U>
	inline Marray<T,rank,base> & operator= (const U &u){
		for( int i=0;i<get_dim1();i++)
			for(int j=0;j<get_dim2();j++)
				(*this)(i,j)=u;
		return *this;

	}

	template <class U>
	inline Marray<T,rank,base> & operator+= (const U &u){
		for( int i=0;i<get_dim1();i++)
			for(int j=0;j<get_dim2();j++)
				(*this)(i,j)+=u;
		return *this;

	}

	template <class U>
	inline Marray<T,rank,base> & operator-= (const U &u){
		for( int i=0;i<get_dim1();i++)
			for(int j=0;j<get_dim2();j++)
				(*this)(i,j)-=u;
		return *this;

	}


	template <class U>
	inline Marray<T,rank,base> & operator*= (const U &u){
		for( int i=0;i<get_dim1();i++)
			for(int j=0;j<get_dim2();j++)
				(*this)(i,j)*=u;
		return *this;

	}


	template <class U>
	inline Marray<T,rank,base> & operator/= (const U &u){
		for( int i=0;i<get_dim1();i++)
			for(int j=0;j<get_dim2();j++)
				(*this)(i,j)/=u;
		return *this;

	}

	//#else
	//	using base::operator=;
	//#endif


	Marray(long* dimensions):base(dimensions){
		(*this)=0;
	}


	//copy constructor
	Marray(const Marray<T, rank,base> &R){
		(*this)=R;
	}

	void sucesion(int init,int stride=1)
	{
		int count=init;
		for(int i=0;i<get_dim1();i++)
		{
			for(int j=0;j<get_dim2();j++)
			{			(*this)(i,j)=count;
			count = count + stride;
			}
		}
	}
	inline const T operator() (const int n1,const int n2) const {
#ifdef CHECK_OOB
		CHECK
#endif
			return base::operator()(n1,n2);

	}

	inline T & operator() (const int n1,const int n2){
#ifdef CHECK_OOB
		CHECK
#endif
			return base::operator()(n1,n2);
	}


	//////////////////////////////////////////////////////////////
	// MEMBER FUNCTIONS
	//////////////////////////////////////////////////////////////



	void resize(long d1,long d2){
		unsigned long  dim[2];
		dim[0]=d1;
		dim[1]=d2;
		base::resize(dim);
		(*this)=(T)0;

	}

	inline int get_size() const{

		return base::size[0]*base::size[1];

	}


	inline int get_dim1() const{
		return base::size[0];

	}

	inline int get_dim2() const{
		return base::size[1];
	}

	inline void show_shape(){
		std::cout<<base::size[0]<<"x"<<base::size[1]<<std::endl;
	};

	void inline operator=(const char *data){

		char delim;
		std::stringstream stream;
		stream<<data;

		for(int i=0;i<get_dim1();i++){
			for(int j=0;j<get_dim2();j++){
				stream>>(*this)(i,j);
			}
			stream>>delim;
			assert(delim==';');
		}





	}

	friend std::ostream & operator << <T>
		(std::ostream & os, const Marray<T, rank,base> & v);


	//////////////////////////////////////////////////////////////
	// Common Matrix Operations
	/////////////////////////////////////////////////////////////

	void Identity(){
		for (unsigned int i=0;i<base::size[0];i++)
			for(unsigned int j=0;j<base::size[1];j++)
				(*this)(i,j)= (T)(i==j?1:0);

	}
private:



	inline void NBackSubstitution(int *nrow, Marray<T,rank,base> &Mtr){

		T acum;

		for (int k=Mtr.size[1]-Mtr.size[0]-1;k>=0;k--){

			for(int i=Mtr.size[0]-1;i>=0;i--){
				acum=Mtr(nrow[i],Mtr.size[0]+k);
				for(int j=Mtr.size[0]-1;j>i;j--){
					acum-=Mtr(nrow[i],j) * (*this)(j,k);
					//cout<<Mtr(nrow[i],j)<<endl;
					//cout<<(*this)(nrow[j],Mtr.size[0]+k)<<endl;
				}
				(*this)(i,k)=acum/Mtr(nrow[i],i);
				//cout<<"i:"<<i<<" " << "k:"<<k<<"="<<(*this)(i,k)<<endl;
			}

		}

	}

	//---------------------------------------
	//Helper function for SVD
	//---------------------------------------
	float pythag(double a, double b)
		//Computes (a2 + b2)1=2 without destructive underflow or overflow.
	{
		double absa,absb;
		absa=abs(a);
		absb=abs(b);
		if (absa > absb) return (T)(absa*sqrt(1.0+pow(absb/absa,2)));
		else return (T)(absb == 0.0 ? 0.0 : absb*sqrt(1.0+pow(absa/absb,2)));
	}


	inline Marray<T,rank,base> NBackSubstitutionRet(int *nrow, Marray<T,rank,base> &Mtr){

		T acum;
		Marray<T,rank,base> ret(Mtr.size[0],Mtr.size[0]);
		for (int k=Mtr.size[1]-Mtr.size[0]-1;k>=0;k--){
			for(int i=Mtr.size[0]-1;i>=0;i--){
				acum=Mtr(nrow[i],Mtr.size[0]+k);
				for(int j=Mtr.size[0]-1;j>i;j--){
					acum-=Mtr(nrow[i],j) * ret(j,k);
					//cout<<Mtr(nrow[i],j)<<endl;
					//cout<<(*this)(nrow[j],Mtr.size[0]+k)<<endl;
				}
				ret(i,k)=acum/Mtr(nrow[i],i);
				//cout<<"i:"<<i<<" " << "k:"<<k<<"="<<(*this)(i,k)<<endl;
			}

		}
		return ret;

	}





public:

	void compute_deviatoric_tensor(Marray<T, rank,base> & v, double & mean_trace)
	{

		Index<'i'> iG;
		Index<'j'> jG;
		mean_trace=((*this)(iG,iG))/get_dim1();
		v(iG,jG)=(*this)(iG,jG);
		for (int j=0;j<get_dim1();j++)
		{	v(j,j)-=mean_trace;
		}
	}
	void convert_to_Tensor1(Marray<T,1> & v)
	{
		int dim1=get_dim1();
		int dim2=get_dim2();
		for (int i=0;i<dim1;i++)
		{
			for (int j=0;j<dim2;j++)
			{
				v(i*dim2+j)=(*this)(i,j);
			}
		}
	}


	Marray<T,rank,base>& RotationMatrixDeg(const double theta1)
	{
		double pi = 2.0*atan2(1,0);
		double theta =pi*theta1/180.0;
		(*this)(0,0) =	cos(theta);
		(*this)(1,0) =	sin(theta);
		(*this)(1,1) =	cos(theta);
		(*this)(0,1) =	-sin(theta);
		(*this)(2,2)= 1.0;
		return *this;
	}
	Marray<T,rank,base>& RotationMatrix(const double theta)
	{
		(*this)(0,0) =	cos(theta);
		(*this)(1,0) =	sin(theta);
		(*this)(1,1) =	cos(theta);
		(*this)(0,1) =	-sin(theta);
		(*this)(2,2)= 1.0;
		return *this;
	}



	// Matrix Invertion Found on the Net
	// Code was transformed to LTensors standards
	//
	// Compute inverse of matrix
	Marray<T,rank,base> Inv(double * determ1 = NULL)
		// Input
		//    A    -    Matrix A (N by N)
		// Outputs
		//   Ainv  -    Inverse of matrix A (N by N)
		//  determ -    Determinant of matrix A	(return value)
	{
		int N = (*this).get_dim1();
		assert( N == (*this).get_dim2() );
		Marray<T,rank,base> A(N,N),Ainv(N,N);
		A = (*this);
		double  determ;
		int i, j, k;
		Marray<T,1 > scale(N); // Scale factor
		Marray<T,rank,base> b(N,N);	 // work array
		int *index;  index = new int [N+1];

		//* Matrix b is initialized to the identity matrix
		b.Identity();

		//* Set scale factor, scale(i) = max( |a(i,j)| ), for each row
		for( i=0; i<N; i++ )
		{
			index[i] = i;			  // Initialize row index list
			T scalemax = 0.;
			for( j=0; j<N; j++ )
			{
				scalemax = (scalemax > fabs(A(i,j))) ? scalemax : fabs(A(i,j));
			}
			scale(i) = scalemax;
			//assert(scalemax>0);
			if(!(scalemax>0))return Ainv;
		}

		//* Loop over rows k = 1, ..., (N-1)
		int signDet = 1;
		for( k=0; k<N-1; k++ ) {
			//* Select pivot row from max( |a(j,k)/s(j)| )
			double ratiomax = 0.0;
			int jPivot = k;
			for( i=k; i<N; i++ ) {
				double ratio = fabs(A(index[i],k))/scale(index[i]);
				if( ratio > ratiomax ) {
					jPivot=i;
					ratiomax = ratio;
				}
			}
			//* Perform pivoting using row index list
			int indexJ = index[k];
			if( jPivot != k )
			{	          // Pivot
				indexJ = index[jPivot];
				index[jPivot] = index[k];   // Swap index jPivot and k
				index[k] = indexJ;
				signDet *= -1;			  // Flip sign of determinant
			}
			assert(fabs(A(indexJ,k))>0);  //matrix is not invertible
			//* Perform forward elimination
			for( i=k+1; i<N; i++ ) {
				T coeff = A(index[i],k)/A(indexJ,k);
				for( j=k+1; j<N; j++ )
					A(index[i],j) -= coeff*A(indexJ,j);
				A(index[i],k) = coeff;
				for( j=0; j<N; j++ )
					b(index[i],j) -= A(index[i],k)*b(indexJ,j);
			}
		}
		//* Compute determinant as product of diagonal elements
		determ = signDet;	   // Sign of determinant
		for( i=0; i<N; i++ )
			determ *= A(index[i],i);

		if (determ1 != NULL)
		{
			*determ1 = determ;
		}
		if (fabs(determ)>1.e-30)
		{
			//* Perform backsubstitution
			for( k=0; k<N; k++ ) {
				Ainv(N-1,k) = b(index[N-1],k)/A(index[N-1],N-1);
				for( i=N-2; i>=0; i--) {
					T sum = b(index[i],k);
					for( j=i+1; j<N; j++ )
						sum -= A(index[i],j)*Ainv(j,k);
					Ainv(i,k) = sum/A(index[i],i);
				}
			}
		} else
		{
			std::cout << "WARNING: matrix is not invertible, determinant = " << determ <<  std::endl;
		}
		delete [] index;	// Release allocated memory
		return(Ainv);
	}

	double compute_Determinant(void)
		// Input
		//    A    -    Matrix A (N by N)
		// Outputs
		//   Ainv  -    Inverse of matrix A (N by N)
		//  determ -    Determinant of matrix A	(return value)
	{
		int N = (*this).get_dim1();
		assert( N == (*this).get_dim2() );
		Marray<T,rank,base> A(N,N);
		A = (*this);

		int i, j, k;
		Marray<T,1 > scale(N); // Scale factor
		Marray<T,rank,base> b(N,N);	 // work array
		int *index;  index = new int [N+1];

		//* Matrix b is initialized to the identity matrix
		b.Identity();

		//* Set scale factor, scale(i) = max( |a(i,j)| ), for each row
		for( i=0; i<N; i++ )
		{
			index[i] = i;			  // Initialize row index list
			double scalemax = 0.;
			for( j=0; j<N; j++ )
			{
				scalemax = (scalemax > fabs(A(i,j))) ? scalemax : fabs(A(i,j));
			}
			scale(i) = scalemax;
			//assert(scalemax>0);
			if(!(scalemax>0))return 0.0;
		}

		//* Loop over rows k = 1, ..., (N-1)
		int signDet = 1;
		for( k=0; k<N-1; k++ ) {
			//* Select pivot row from max( |a(j,k)/s(j)| )
			double ratiomax = 0.0;
			int jPivot = k;
			for( i=k; i<N; i++ ) {
				double ratio = fabs(A(index[i],k))/scale(index[i]);
				if( ratio > ratiomax ) {
					jPivot=i;
					ratiomax = ratio;
				}
			}
			//* Perform pivoting using row index list
			int indexJ = index[k];
			if( jPivot != k )
			{	          // Pivot
				indexJ = index[jPivot];
				index[jPivot] = index[k];   // Swap index jPivot and k
				index[k] = indexJ;
				signDet *= -1;			  // Flip sign of determinant
			}
			//assert(fabs(A(indexJ,k))>0);  //matrix is not invertible
			if(!fabs(A(indexJ,k))>0) return 0.0;
			//* Perform forward elimination
			for( i=k+1; i<N; i++ ) {
				double coeff = A(index[i],k)/A(indexJ,k);
				for( j=k+1; j<N; j++ )
					A(index[i],j) -= coeff*A(indexJ,j);
				A(index[i],k) = coeff;
				for( j=0; j<N; j++ )
					b(index[i],j) -= A(index[i],k)*b(indexJ,j);
			}
		}
		//* Compute determinant as product of diagonal elements
		double determ = signDet;	   // Sign of determinant
		for( i=0; i<N; i++ )
			determ *= A(index[i],i);

		delete [] index;	// Release allocated memory
		return(determ);
	}

	Marray<T,1> SolvebyGaussElim(const Marray<T,1> & b, bool *err = NULL)
		// Solve system of N linear equations with N unknowns
		// using Gaussian elimination with scaled partial pivoting
		// err returns true if process fails;
	{
		int N = (*this).get_dim1();
		assert( N == (*this).get_dim2() );
		Marray<T,rank,base> A(N,N);
		Marray<T,1> b1(N),result(N);
		Index<'i'> iG;
		Index<'j'> jG;
		A(iG,jG)=(*this)(iG,jG);
		b1(iG)=b(iG);
		Marray<int,1> indx(N);
		Marray<T,1> scale(N);
		double sum;
		int maxIndx;
		int tmpIndx;
		double ratio;
		double maxRatio;
		// index array initialization
		for (int i = 0; i < N; i++) indx(i) = i;

		// determine scale factors
		for (int row = 0; row < N; row++)
		{
			scale(row) = absTrait<T>::function((A(row,0)));
			for (int col = 1; col < N; col++)
			{
				if (absTrait<T>::function(A(row,col)) > scale(row)) scale(row) = absTrait<T>::function(A(row,col));
			}
		}

		// forward elimination
		for (int k = 0; k < N; k++)
		{
			// determine index of pivot row
			maxRatio = absTrait<T>::function(A(indx(k), k))/scale(indx(k));
			maxIndx = k;
			for (int i = k+1; i < N; i++)
			{
				if (absTrait<T>::function(A(indx(i),k))/scale(indx(i)) > maxRatio)
				{
					maxRatio = absTrait<T>::function(A(indx(i),k))/scale(indx(i));
					maxIndx = i;
				}
			}
			if (maxRatio == 0) // no pivot available
			{
				std::cout << "WARNING: matrix is not invertible linear system can not be solved"<<std::endl;
				if (err != NULL)
				{*err = true;}
				return result;
			}
			tmpIndx =indx(k); indx(k)=indx(maxIndx); indx(maxIndx) = tmpIndx;

			// use pivot row to eliminate kth variable in "lower" rows
			for (int i = k+1; i < N; i++)
			{
				ratio = -A(indx(i),k)/A(indx(k),k);
				for (int col = k; col <= N-1; col++)
				{
					A(indx(i), col) += ratio*A(indx(k),col);
				}
				b1(indx(i)) += ratio*b1(indx(k));
			}
		}

		// back substitution
		for (int k = N-1; k >= 0; k--)
		{
			sum = 0;

			for (int col = k+1; col < N; col++)
			{
				sum += A(indx(k),col) * b1(indx(col));
			}
			b1(indx(k)) = (b1(indx(k)) - sum)/A(indx(k),k);
		}

		/*
		cout << endl;
		for (int r=0; r<N; r++)
		{
		cout << indx[r];
		for (int c=0; c<=N; c++) cout<<"  " << A[r][c];
		cout << endl;
		}
		*/

		for (int k = 0; k < N; k++)
		{ result(k) = b1(indx(k));
		}
		return result;
	}

	Marray<T,rank,base> SolvebyGaussElim(const Marray<T,rank,base> & b, bool *err = NULL)
		// Solve system of N linear equations with N unknowns
		// using Gaussian elimination with scaled partial pivoting
		// First N rows and N+1 columns of A contain the system
		// with right-hand sides of equations in column N+1
		// err returns true if process fails; false if it is successful
		// original contents of A are destroyed
		// solution appears in column N
		// obtained from http://www.daniweb.com/forums/thread41803.html
	{
		int N = (*this).get_dim1();
		assert( N == (*this).get_dim2() );
		int Nrhs = b.get_dim2();
		Marray<T,rank,base> A(N,N),b1(N,Nrhs),result(N,Nrhs);
		Index<'i'> iG;
		Index<'j'> jG;
		A(iG,jG)=(*this)(iG,jG);
		b1(iG,jG)=b(iG,jG);
		Marray<int,1> indx(N);
		Marray<T,1> scale(N), sum(N);
		int maxIndx;
		int tmpIndx;
		double ratio;
		double maxRatio;
		// index array initialization
		for (int i = 0; i < N; i++) indx(i) = i;

		// determine scale factors
		for (int row = 0; row < N; row++)
		{
			scale(row) = absTrait<T>::function((A(row,0)));
			for (int col = 1; col < N; col++)
			{
				if (absTrait<T>::function(A(row,col)) > scale(row)) scale(row) = absTrait<T>::function(A(row,col));
			}
		}

		// forward elimination
		for (int k = 0; k < N; k++)
		{
			// determine index of pivot row
			maxRatio = absTrait<T>::function(A(indx(k), k))/scale(indx(k));
			maxIndx = k;
			for (int i = k+1; i < N; i++)
			{
				if (absTrait<T>::function(A(indx(i),k))/scale(indx(i)) > maxRatio)
				{
					maxRatio = absTrait<T>::function(A(indx(i),k))/scale(indx(i));
					maxIndx = i;
				}
			}
			if (maxRatio == 0) // no pivot available
			{
				std::cout << "WARNING: matrix is not invertible linear system can not be solved"<<std::endl;
				if (err != NULL)
				{*err = true;}
				return result;
			}
			tmpIndx =indx(k); indx(k)=indx(maxIndx); indx(maxIndx) = tmpIndx;

			// use pivot row to eliminate kth variable in "lower" rows
			for (int i = k+1; i < N; i++)
			{
				ratio = -A(indx(i),k)/A(indx(k),k);
				for (int col = k; col <= N-1; col++)
				{
					A(indx(i), col) += ratio*A(indx(k),col);
				}
				for (int col = 0; col <= Nrhs-1; col++)
				{
					b1(indx(i), col) += ratio*b1(indx(k),col);
				}
			}
		}

		// back substitution
		for (int j = 0; j<Nrhs; j++)
		{
			for (int k = N-1; k >= 0; k--)
			{
				sum(j) = 0;

				for (int col = k+1; col < N; col++)
				{
					sum(j) += A(indx(k),col) * b1(indx(col),j);
				}
				b1(indx(k),j) = (b1(indx(k),j) - sum(j))/A(indx(k),k);
			}
		}

		/*
		cout << endl;
		for (int r=0; r<N; r++)
		{
		cout << indx[r];
		for (int c=0; c<=N; c++) cout<<"  " << A[r][c];
		cout << endl;
		}
		*/
		for (int j = 0; j<Nrhs; j++)
		{
			for (int k = 0; k < N; k++)
			{ result(k,j) = b1(indx(k),j);
			}
		}
		return result;
	}



	bool LU(Marray<int,1> &indx, float &d){
		//Performs LU decomposition
		//indx its the record of row permutations vector
		//d is output assuming +-1 depending if the number of row interchanges
		//was even or odd respectively

		int i,imax,j,k;
		double big,dum,sum,temp;
		int n=get_dim1();
		indx.resize(n);
		Marray<double,1> vv(n);  //vv stores the implicit scaling of each row.
		d=1.0;                          //no row interchanges yet

		for(i=0;i<n;i++){      //Loop over rows to get the implicit scaling information
			big=0.0;
			for(j=0;j<n;j++)
				if( (temp=absTrait<T>::function((*this)(i,j)))>big) big=temp;
			if(big==0.0) return false;
			//No nonzero largest element
			vv(i)=1.0/big;              //save the scaling
		}
		for(j=0;j<n;j++){      //this is loop over columns of Crout's method
			for(i=0;i<j;i++){
				sum=(*this)(i,j);
				for(k=0;k<i;k++)
					sum-=(*this)(i,k)*(*this)(k,j);
				(*this)(i,j)=sum;
			}
			big=0.0;            //initialize for the search for largest pivot element
			for(i=j;i<n;i++){
				sum=(*this)(i,j);
				for(k=0;k<j;k++)
					sum-=(*this)(i,k)*(*this)(k,j);
				(*this)(i,j)=sum;
				if( (dum=vv(i)*absTrait<T>::function(sum) ) >= big) {
					//Is the figure of merit for the pivot bettern than the best so far?
					big=dum;
					imax=i;
				}
			}
			if(j!=imax){        //do we need to interchange rows?

				for(k=0;k<n;k++){   //yes, do so
					dum=(*this)(imax,k);
					(*this)(imax,k)=(*this)(j,k);
					(*this)(j,k)=dum;
				}
				d=-d;               //change parity
				vv(imax)=vv(j);     //interchange the scale factor
			}
			indx(j)=imax;
			if((*this)(j,j)==0.0) (*this)(j,j)=1.0e-20;
			//if the pivot element is ZERO the matrix is singular.
			if(j!=n-1){
				dum=1.0/((*this)(j,j));       //divide by pivot
				for(i=j+1;i<n;i++)
					(*this)(i,j) *=dum;
			}
		}   //go back for the next column in the reduction

		return true;
	}

	bool getLU(Marray<T,rank> &A, Marray<int,1> &indx, float *d=NULL){
		A=(*this);
		float dl;
		if(!A.LU(indx,dl)) return false;
		if(d!=NULL)
			*d=dl;
		return true;
	}

	void LU_backSubs(const Marray<int,1> &indx, Marray<T,1> &b ){
		//This algoritm assumes this is an LU factorized matrix
		int i,ii=0,ip,j;
		T sum;

		int n=get_dim1();
#ifdef CHECK_OOB
		assert(b.get_dim1()==n);
#endif

		for(i=1;i<=n;i++){
			ip=indx(i-1);
			sum=b(ip);
			b(ip)=b(i-1);
			if(ii)
				for(j=ii;j<=i-1;j++)
					sum-=(*this)(i-1,j-1)*b(j-1);
			else if(sum) ii=i;
			b(i-1)=sum;
		}
		for(i=n;i>=1;i--){
			sum=b(i-1);
			for(j=i+1;j<=n;j++) sum-=(*this)(i-1,j-1)*b(j-1);
			b(i-1)=sum/(*this)(i-1,i-1);
		}

	}

	bool  Cholesky(Marray<T,1> &p){
		//P represents the diagonal of the decomposed matrix
		//while in the lower diagonal of this is L (Lt*L=A)
		int i,j,k;
		double sum;
		int n=get_dim1();
		p.resize(n);

		for(i=0;i<n;i++){
			for(j=i;j<n;j++){
				for(sum=(*this)(i,j) , k=i-1 ; k>=0;k--)
					sum-=(*this)(i,k)*(*this)(j,k);
				if(i==j){
					if(sum<=0.0) //not positive definite
						if (sum<=0.0) return false;
					p(i)=sqrtTrait<T>::function(sum);
				}else{
					(*this)(j,i)=sum/p(i);
				}
			}
		}
		return true;
	}

	bool getCholesky(Marray<T,rank> &A,Marray<T,1> &p){
		A=(*this);
		if(!A.Cholesky(p)) return false;
		return true;
	}

	bool InvertByLU(Marray<T,rank> &y){
		//returns the inverse in y
		//keeps LU in this

		int n=get_dim1();
		y.resize(n,n);
		Marray<int,1> indx(n);
		Marray<T,1> col(n);
		Index <'i'> iG;
		float d;
		int j;

		//Decompose the matrix just once.
		if( !LU(indx,d) ) return false;

		for(j=0;j<n;j++) { //Find inverse by columns.
			col=0.0;
			col(j)=1.0;
			LU_backSubs(indx,col);
			y(iG,j)=col(iG);
		}

		return true;

	}


	bool getInvertByLU(Marray<T,rank> &A){
		Marray<T,rank> temp;
		temp=(*this);
		if(!temp.InvertByLU(A)) return false;

		return true;
	}




	void saveToSpacedFile(const char* fileName){

		std::ofstream out(fileName,std::ios::out|std::ios::trunc);

#ifdef CHECK_FILES
		assert(out.is_open());
#endif
		out<<get_dim1()<<" "<<get_dim2()<<std::endl;

		for(int i=0;i<get_dim1();i++){
			for(int j=0;j<get_dim2();j++){
				out<<(*this)(i,j)<<" ";
			}
			out<<std::endl;
		}
		out.close();


	}

	//Column x Rows ordered data

	void loadFromSpacedFile(const char* fileName){

		std::ifstream in(fileName);

#ifdef CHECK_FILES
		assert(in.is_open());
#endif
		int dim1,dim2;

		in>>dim1;
		in>>dim2;
#ifdef CHECK_FILES
		assert(dim1>0);

		assert(dim2>0);

#endif
		resize(dim1,dim2);
		for(int i=0;i<dim1;i++)
			for(int j=0;j<dim2;j++)
				in>>(*this)(i,j);

		in.close();

	}



	///////////////////////////////////////////////////////////////
	// create IndexedExpressions
	///////////////////////////////////////////////////////////////


	//F_Expres

	inline F_Expr2 < Encapsulate_to_Expr2< Marray<T,rank,base>,T,rank,12,IndexF,IndexF>, T>
		operator() (const IndexF  & index1,const IndexF & index2)
	{

		typedef Encapsulate_to_Expr2<Marray<T,rank,base>,T,rank,12,IndexF,IndexF> Expr_Obj;
		return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,index1,index2));
	}

	inline F_Expr2 < Encapsulate_to_Expr2<const Marray<T,rank,base>,T,rank,12,IndexF,IndexF>, T>
		operator() (const IndexF  & index1,const IndexF & index2) const
	{
		typedef Encapsulate_to_Expr2<const Marray<T,rank,base>,T,rank,12,IndexF,IndexF> Expr_Obj;
		return F_Expr2<Expr_Obj,T>(Expr_Obj(*this,index1,index2));
	}



	inline F_Expr1 < Encapsulate_to_Expr1<Marray<T,rank,base>,T,rank,1,IndexF>, T>
		operator() (const IndexF & index1, const int N)
	{
		typedef Encapsulate_to_Expr1<Marray<T,rank,base>,T,rank,1,IndexF> Expr_Obj;
		return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N,index1));
	}
	inline F_Expr1 < Encapsulate_to_Expr1<const Marray<T,rank,base>,T,rank,1,IndexF>, T>
		operator() (const IndexF & index1, const int N) const
	{
		typedef Encapsulate_to_Expr1<const Marray<T,rank,base>,T,rank,1,IndexF> Expr_Obj;
		return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N,index1));
	}


	inline F_Expr1 < Encapsulate_to_Expr1<Marray<T,rank,base>,T,rank,2,IndexF>, T>
		operator() (const int N, const IndexF & index1 )
	{
		typedef Encapsulate_to_Expr1<Marray<T,rank,base>,T,rank,2,IndexF> Expr_Obj;
		return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N,index1));
	}
	inline F_Expr1 < Encapsulate_to_Expr1<const Marray<T,rank,base>,T,rank,2,IndexF>, T>
		operator() (const int N, const IndexF & index1 ) const
	{
		typedef Encapsulate_to_Expr1<const Marray<T,rank,base>,T,rank,2,IndexF> Expr_Obj;
		return F_Expr1<Expr_Obj,T>(Expr_Obj(*this,N,index1));
	}


	//end FExpres

	template < char i, char j , int iType1, int iType2>
	inline Expr2 < Encapsulate_to_Expr2<Marray<T,rank,base>,T,rank,12,Index<i,iType1>,Index<j,iType2> >, T, i ,j >
		operator() (const Index < i , iType1> & index1,const Index < j ,iType2> & index2)
	{

		typedef Encapsulate_to_Expr2<Marray<T,rank,base>,T,rank,12,Index<i,iType1> , Index <j,iType2 > >Expr_Obj;
		return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,index1,index2));
	}

	template < char i, char j , int iType1, int iType2>
	inline Expr2 < Encapsulate_to_Expr2<const Marray<T,rank,base>,T,rank,12,Index<i,iType1>,Index<j,iType2> >, T, i ,j >
		operator() (const Index < i , iType1> & index1,const Index <j ,iType2> & index2)const
	{

		typedef Encapsulate_to_Expr2<const Marray<T,rank,base>,T,rank,12,Index<i,iType1> , Index <j, iType2 > >Expr_Obj;
		return Expr2<Expr_Obj,T,i,j>(Expr_Obj(*this,index1,index2));
	}

	//checked


	template < char i, int iType1>
	inline Expr1 < Encapsulate_to_Expr1<Marray<T,rank,base>,T,rank,1,Index<i,iType1> >, T, i>
		operator() (const Index < i, iType1 > & index1, const int N)
	{
		typedef Encapsulate_to_Expr1<Marray<T,rank,base>,T,rank,1,Index<i,iType1> > Expr_Obj;
		return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N,index1));
	}

	template < char i, int iType1>
	inline Expr1 < Encapsulate_to_Expr1<const Marray<T,rank,base>,T,rank,1,Index<i,iType1> >, T, i>
		operator() (const Index < i, iType1 > & index1, const int N)const
	{
		typedef Encapsulate_to_Expr1<const Marray<T,rank,base>,T,rank,1,Index<i,iType1> > Expr_Obj;
		return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N,index1));
	}


	//checked

	template < char i, int iType2>
	inline Expr1 < Encapsulate_to_Expr1<Marray<T,rank,base>,T,rank,2, Index<i,iType2> >, T, i>
		operator() (const int N, const Index < i, iType2 > & index2 )
	{
		typedef Encapsulate_to_Expr1<Marray<T,rank,base>,T,rank,2,Index<i,iType2> > Expr_Obj;
		return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N,index2));
	}

	template < char i, int iType2>
	inline Expr1 < Encapsulate_to_Expr1<const Marray<T,rank,base>,T,rank,2, Index<i,iType2> >, T, i>
		operator() (const int N, const Index < i, iType2 > & index2 )const
	{
		typedef Encapsulate_to_Expr1<const Marray<T,rank,base>,T,rank,2,Index<i,iType2> > Expr_Obj;
		return Expr1<Expr_Obj,T,i>(Expr_Obj(*this,N,index2));
	}

	//checked

	template < char i , int iType1, int iType2>
	inline T  operator() (const Index < i,iType1 > & index1,const Index < i,iType2 > & index2)
	{
		T res = 0;
		for (int n1 = 0; n1 < get_dim1(); ++n1)
		{
			res += (*this)(n1,n1);
		}
		return res;
	}

	template < char i , int iType1, int iType2>
	inline T  operator() (const Index < i,iType1 > & index1,const Index < i,iType2 > & index2) const
	{
		T res = 0;
		for (int n1 = 0; n1 < get_dim1(); ++n1)
		{
			res += (*this)(n1,n1);
		}
		return res;
	}

	//added 02/06





	/////////////////////////////////////////////////////////////
	// Assignation Operations from other Marray<T,rank,base> object
	/////////////////////////////////////////////////////////////

	template <class U,class base2>
	inline Marray<T,rank,base> & operator+= (const Marray<U,rank,base2> & a)
	{
		for(int i=0;i<get_dim1();i++)
		{
			for(int j=0;j<get_dim2();j++)
			{
				(*this)(i,j)+=a(i,j);
			}
		}
		return *this;
	}



	//Implemented for compatibility only, shouldn't be used in performance critic sections
	//will fix this with the new cxx specification

	inline Marray<T,rank,base> operator-()
	{

		Marray<T,rank,base> temp(get_dim1(),get_dim2());

			for(int i=0;i<get_dim1();i++)
		{
			for(int j=0;j<get_dim2();j++)
			{
				temp(i,j)=-(*this)(i,j);
			}
		}

		return temp;
	}



	//Compatibily operators
	//cstyle ordering

	template <class Type>
	void fromCArray(Type * datos,int dim1,int dim2){
		resize(dim1,dim2);

		for(int i=0;i<dim1;i++)
			for(int j=0;j<dim2;j++)
				(*this)(i,j)=datos[i*dim1+j];
	}

	template <class Type>
	void toCArray(Type * datos){


		int dim1=get_dim1();
		int dim2=get_dim2();
		for(int i=0;i<dim1;i++)
			for(int j=0;j<dim2;j++)
				datos[i*dim1+j]=(*this)(i,j);
	}

	template <class U,class base2>
	inline Marray<T,rank,base> & operator-= (const Marray<U,rank,base2> & a)
	{
		for(int i=0;i<get_dim1();i++)
		{
			for(int j=0;j<get_dim2();j++)
			{
				(*this)(i,j)-=a(i,j);
			}
		}
		return *this;
	}

	bool operator==(const Marray<T,rank,base> &a){
		if(get_dim1()!=a.get_dim1())
			return false;
		if(get_dim2()!=a.get_dim2())
			return false;

		for(int i=0;i<get_dim1();i++){
			for(int j=0;j<get_dim2();j++){

				if( (*this)(i,j)!=a(i,j))
					return false;
			}

		}
		return true;

	}


//------------------------------------------------------------------------------------------
//Computes the Moore-Penrose pseudo inverse of the matrix using SVDSorted
//------------------------------------------------------------------------------------------
	Marray<T,2,base> PInv(){

		Marray<float,2> U;
		Marray<float,2> V;
		Marray<float,2> W;
		Marray<float,2> inv;
		Marray<float,2> temp;

		Index<'i'> i;
		Index<'j'> j;
		Index<'k'> k;
		int sd;
		double eps=1e-6;
		//Compute the SVD decomposition with no sorting
		SVD(U,W,V);
		sd=W.get_dim1();
		//Compute the diagonal pseudo inverse of W
		//by taking the reciprocal of each diagonal element different from 0
		//and leaving the 0s be.
		for(int l=0;l<sd;++l)
			W(l,l)= (T) (W(l,l)>eps? 1.0/W(l,l) : 0.0);

		inv.resize(get_dim2(),get_dim1());
		temp.resize(V.get_dim1(),V.get_dim1());
		//The Moore penrose pseudo inverse can be achivied by doing
		//Pinv=V*W^*V''; where W^ is the pseudo inverse of W

		temp(i,k)=V(i,j)*W(k,j);
		inv(i,k)=temp(i,j)*U(k,j);

		return inv;


	}
//------------------------------------------------------------------------------------------
//Given the matrix A stored in u[0..m-1][0..n-1], this routine computes its singular value
//decomposition, A D U W  VT and stores the results in the matrices u and v, and the vector
//w in descending order.
//------------------------------------------------------------------------------------------
	void SVDSorted( Marray<T,2,base> &U, Marray<T,2,base> &W, Marray<T,2,base> &V){
		SVD(U,W,V);
		SVDReorder(U,W,V);
	}
//------------------------------------------------------------------------------------------
//Given the matrix A stored in u[0..m-1][0..n-1], this routine computes its singular value
//decomposition, A D U W  VT and stores the results in the matrices u and v, and the vector
//w.
//------------------------------------------------------------------------------------------
void SVD( Marray<T,2,base> &U, Marray<T,2,base> &W, Marray<T,2,base> &V)
{
	U=*this;
	int m=get_dim1();
	int n=get_dim2();
	W.resize(n,n);
	V.resize(n,n);
	double eps=1e-10;

	bool flag;
	int i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z;
	Marray<T,1,base> rv1(n);
	g = scale = anorm = 0.0; //Householder reduction to bidiagonal form.
	for (i=0;i<n;i++) {
		l=i+2;
		rv1(i)=(T)scale*(T)g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += abs(U(k,i));
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					U(k,i) /= (T)scale;
					s += U(k,i)*U(k,i);
				}
				f=U(i,i);
				g = -copysign(sqrt(s),f);
				h=f*g-s;
				U(i,i)=(T)(f-g);
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += U(k,i)*U(k,j);
					f=s/h;
					for (k=i;k<m;k++) U(k,j) += (T)f*U(k,i);
				}
				for (k=i;k<m;k++) U(k,i) *=(T) scale;
			}
		}
		W(i,i)=(T)scale *(T)g;
		g=s=scale=0.0;
		if (i+1 <= m && i+1 != n) {
			for (k=l-1;k<n;k++) scale += abs(U(i,k));
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					U(i,k) /= (T)scale;
					s += U(i,k)*U(i,k);
				}
				f=U(i,l-1);
				g = -copysign(sqrt(s),f);
				h=f*g-s;
				U(i,l-1)=(T)(f-g);
				for (k=l-1;k<n;k++) rv1(k)=(T)U(i,k)/(T)h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += U(j,k)*U(i,k);
					for (k=l-1;k<n;k++) U(j,k) +=(T) s*rv1(k);
				}
				for (k=l-1;k<n;k++) U(i,k) *= (T)scale;
			}
		}
		anorm=std::max((double)anorm,(double)(abs(W(i,i))+abs(rv1(i))));
	}
	for (i=n-1;i>=0;i--) { //Accumulation of right-hand transformations.
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++) //Double division to avoid possible underflow.
					V(j,i)=(T)(U(i,j)/U(i,l))/(T)g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += U(i,k)*V(k,j);
					for (k=l;k<n;k++) V(k,j) +=(T) s*V(k,i);
				}
			}
			for (j=l;j<n;j++) V(i,j)=V(j,i)=0.0;
		}
		V(i,i)=1.0;
		g=rv1(i);
		l=i;
	}
	for (i=std::min(m,n)-1;i>=0;i--) {// Accumulation of left-hand transformations.
		l=i+1;
		g=W(i,i);
		for (j=l;j<n;j++) U(i,j)=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += U(k,i)*U(k,j);
				f=(s/U(i,i))*g;
				for (k=i;k<m;k++) U(k,j) +=(T) f*U(k,i);
			}
			for (j=i;j<m;j++) U(j,i) *= (T)g;
		} else for (j=i;j<m;j++) U(j,i)=0.0;
		U(i,i)=U(i,i)+1;
	}
	for (k=n-1;k>=0;k--) {// Diagonalization of the bidiagonal form: Loop over
		for (its=0;its<30;its++) { //singular values, and over allowed iterations.
			flag=true;
			for (l=k;l>=0;l--) { //Test for splitting.
				nm=l-1;
				if (l == 0 || abs(rv1(l)) <= eps*anorm) {
					flag=false;
					break;
				}
				if (abs(W(nm,nm)) <= eps*anorm) break;
			}
			if (flag) {
				c=0.0; //Cancellation of rv1[l], if l > 0.
				s=1.0;
				for (i=l;i<k+1;i++) {
					f=s*rv1(i);
					rv1(i)=(T)c*rv1(i);
					if (abs(f) <= eps*anorm) break;
					g=W(i,i);
					h=pythag(f,g);
					W(i,i)=(T)h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=U(j,nm);
						z=U(j,i);
						U(j,nm)=(T)(y*c+z*s);
						U(j,i)=(T)(z*c-y*s);
					}
				}
			}
			z=W(k,k);
			if (l == k) { //Convergence.
				if (z < 0.0) { //Singular value is made nonnegative.
					W(k,k) = (T)-z;
					for (j=0;j<n;j++) V(j,k) = -V(j,k);
				}
				break;
			}
			if (its == 29) throw("no convergence in 30 svdcmp iterations");
			x=W(l,l); //Shift from bottom 2-by-2 minor.
			nm=k-1;
			y=W(nm,nm);
			g=rv1(nm);
			h=rv1(k);
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+copysign(g,f)))-h))/x;
			c=s=1.0;	//Next QR transformation:
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1(i);
				y=W(i,i);
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1(j)=(T)z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=V(jj,j);
					z=V(jj,i);
					V(jj,j)=(T)(x*c+z*s);
					V(jj,i)=(T)(z*c-x*s);
				}
				z=pythag(f,h);
				W(j,j)=(T)z;		//Rotation can be arbitrary if z D 0.
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=U(jj,j);
					z=U(jj,i);
					U(jj,j)=(T)(y*c+z*s);
					U(jj,i)=(T)(z*c-y*s);
				}
			}
			rv1(l)=(T)0.0;
			rv1(k)=(T)f;
			W(k,k)=(T)x;
		}
	}

}
//-------------------------------------------------------------------------------------------------
//Given the output of SVD, this routine sorts the singular values, and corresponding columns
//of u and v, by decreasing magnitude. Also, signs of corresponding columns are
//ipped so as to
//maximize the number of positive elements.
//-------------------------------------------------------------------------------------------------
void SVDReorder(Marray<T,2,base> &U, Marray<T,2,base> &W, Marray<T,2,base> &V) {

	int i,j,k,s,inc=1;
	double sw;
	int m=U.get_dim1();
	int n=U.get_dim2();
	Marray<double,1> su(m), sv(n);
	do { inc *= 3; inc++; } while (inc <= n); //Sort. The method is Shell's sort.
											  //(The work is negligible as compared
											  //to that already done in
		                                      //decompose.)
		do {
			inc /= 3;
			for (i=inc;i<n;i++) {
				sw = W(i,i);
				for (k=0;k<m;k++) su(k) = U(k,i);
				for (k=0;k<n;k++) sv(k) = V(k,i);
				j = i;
				while (W(j-inc,j-inc) < sw) {
					W(j,j) = W(j-inc,j-inc);
					for (k=0;k<m;k++) U(k,j) = U(k,j-inc);
					for (k=0;k<n;k++) V(k,j) = V(k,j-inc);
					j -= inc;
					if (j < inc) break;
				}
				W(j,j) = (T)sw;
				for (k=0;k<m;k++) U(k,j) =(T) su(k);
				for (k=0;k<n;k++) V(k,j) =(T) sv(k);
			}
		} while (inc > 1);
	for (k=0;k<n;k++) { //Flip signs.
		s=0;
	for (i=0;i<m;i++) if (U(i,k) < 0.) s++;
	for (j=0;j<n;j++) if (V(j,k) < 0.) s++;
	if (s > (m+n)/2) {
		for (i=0;i<m;i++) U(i,k) = -U(i,k);
		for (j=0;j<n;j++) V(j,k) = -V(j,k);
	}
	}
}



};




template <class Type, class base>
std::ostream & operator<< (std::ostream & os,const  Marray<Type,rank,base> & v)
{
	std::cout<<"MArray2["<< v.get_dim1() <<"," << v.get_dim2() <<"] = " << std::endl<<"\t[ "<<std::endl;

	if ((v.get_dim1()*v.get_dim2())>0)
	{
		for( int i=0;i<v.get_dim1();i++)
		{
			for(int j=0;j<v.get_dim2();j++)
				std::cout<<std::setw(10)<<std::setprecision(5) << v(i,j) <<" ";
			std::cout<<std::endl;
		}
	}
	std::cout<<"\t]"<<std::endl;

	return os;
}
#undef CHECK
#undef rank

#endif

