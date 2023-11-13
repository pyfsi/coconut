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
/////////////////////////////////////////////////////
// Given A(j,i,k) define a new permuted tensor object
// so to define: B(i,j,k)=A(j,i,k)
/////////////////////////////////////////////////////

#ifndef Permute_TensorsRank3_H
#define Permute_TensorsRank3_H

////////////////////////////////////////////
// Define B(i,j,k)=A(i,j,k)  ==> No permutation
///////////////////////////////////////////

template < class A, class T> 
class Aijk_to_Aijk
{
	const A TA;
      public:

	Aijk_to_Aijk (const A &a): TA(a)
	{
	}

	int get_dim1 () const
	{
		return TA.get_dim1();
	}
	
	int get_dim2 () const
	{
		return TA.get_dim2();
	}

	int get_dim3 () const
	{
		return TA.get_dim3();
	}	

	T operator () (const int N1, const int N2, const int N3) const
	{
		return TA(N1,N2,N3);
	}
};

////////////////////////////////////////////
// Define B(i,j,k)=A(i,k,j)  
///////////////////////////////////////////

template < class A, class T> 
class Aikj_to_Aijk
{
	const A TA;
      public:

	Aikj_to_Aijk (const A &a): TA(a)
	{
	}

	int get_dim1 () const
	{
		return TA.get_dim1();
	}
	
	int get_dim2 () const
	{
		return TA.get_dim3();
	}

	int get_dim3 () const
	{
		return TA.get_dim2();
	}	

	T operator () (const int N1, const int N2, const int N3) const
	{
		return TA(N1,N3,N2);
	}
};

////////////////////////////////////////////
// Define B(i,j,k)=A(j,i,k)  
///////////////////////////////////////////

template < class A, class T> 
class Ajik_to_Aijk
{
	const A TA;
      public:

	Ajik_to_Aijk (const A &a): TA(a)
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

	int get_dim3 () const
	{
		return TA.get_dim3();
	}	

	T operator () (const int N1, const int N2, const int N3) const
	{
		return TA(N2,N1,N3);
	}
};

////////////////////////////////////////////
// Define B(i,j,k)=A(k,i,j)  
///////////////////////////////////////////

template < class A, class T> 
class Akij_to_Aijk
{
	const A TA;
      public:

	Akij_to_Aijk (const A &a): TA(a)
	{
	}

	int get_dim1 () const
	{
		return TA.get_dim2();
	}
	
	int get_dim2 () const
	{
		return TA.get_dim3();
	}

	int get_dim3 () const
	{
		return TA.get_dim1();
	}	

	T operator () (const int N1, const int N2, const int N3) const
	{
		return TA(N3,N1,N2);
	}
};

////////////////////////////////////////////
// Define B(i,j,k)=A(j,k,i)  
///////////////////////////////////////////

template < class A, class T> 
class Ajki_to_Aijk
{
	const A TA;
      public:

	Ajki_to_Aijk (const A &a): TA(a)
	{
	}

	int get_dim1 () const
	{
		return TA.get_dim3();
	}
	
	int get_dim2 () const
	{
		return TA.get_dim1();
	}

	int get_dim3 () const
	{
		return TA.get_dim2();
	}	

	T operator () (const int N1, const int N2, const int N3) const
	{
		return TA(N2,N3,N1);
	}
};

////////////////////////////////////////////
// Define B(i,j,k)=A(k,j,i)  
///////////////////////////////////////////

template < class A, class T> 
class Akji_to_Aijk
{
	const A TA;
      public:

	Akji_to_Aijk (const A &a): TA(a)
	{
	}

	int get_dim1 () const
	{
		return TA.get_dim3();
	}
	
	int get_dim2 () const
	{
		return TA.get_dim2();
	}

	int get_dim3 () const
	{
		return TA.get_dim1();
	}	

	T operator () (const int N1, const int N2, const int N3) const
	{
		return TA(N3,N2,N1);
	}
};





#endif
