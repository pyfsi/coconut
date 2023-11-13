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
// Given A(j,i,k,l) define a new permuted tensor object
// so to define: B(i,j,k,l)=A(j,i,k,l)
/////////////////////////////////////////////////////

#ifndef Permute_TensorsRank4_H
#define Permute_TensorsRank4_H

////////////////////////////////////////////
// Define B(i,j,k,l)=A(i,j,k,l)  ==> No permutation
///////////////////////////////////////////

template < class A, class T>
class Aijkl_to_Aijkl
{
	const A TA;
      public:

	Aijkl_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim1();
	}

	int get_dim2() const
	{
		return TA.get_dim2();
	}

	int get_dim3() const
	{
		return TA.get_dim3();
	}

	int get_dim4() const
	{
		return TA.get_dim4();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N1,N2,N3,N4);
	}
};

///////////////////////////////////////////////////////////////////////
/// Define "MINOR" permutations
///////////////////////////////////////////////////////////////////////

////////////////////////////////////////////
// Define B(i,j,k,l)=A(j,i,k,l)
///////////////////////////////////////////

template < class A, class T>
class Ajikl_to_Aijkl
{
	const A TA;
      public:

	Ajikl_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim2();
	}

	int get_dim2() const
	{
		return TA.get_dim1();
	}

	int get_dim3() const
	{
		return TA.get_dim3();
	}

	int get_dim4() const
	{
		return TA.get_dim4();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N2,N1,N3,N4);
	}
};


////////////////////////////////////////////
// Define B(i,j,k,l)=A(i,l,k,j)
///////////////////////////////////////////

template < class A, class T>
class Ailkj_to_Aijkl
{
	const A TA;
      public:

	Ailkj_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim1();
	}

	int get_dim2() const
	{
		return TA.get_dim4();
	}

	int get_dim3() const
	{
		return TA.get_dim3();
	}

	int get_dim4() const
	{
		return TA.get_dim2();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N1,N4,N3,N2);
	}
};


////////////////////////////////////////////
// Define B(i,j,k,l)=A(i,j,l,k)
///////////////////////////////////////////

template < class A, class T>
class Aijlk_to_Aijkl
{
	const A TA;
      public:

	Aijlk_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim1();
	}

	int get_dim2() const
	{
		return TA.get_dim2();
	}

	int get_dim3() const
	{
		return TA.get_dim4();
	}

	int get_dim4() const
	{
		return TA.get_dim3();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N1,N2,N4,N3);
	}
};

////////////////////////////////////////////
// Define B(i,j,k,l)=A(j,i,l,k)
///////////////////////////////////////////

template < class A, class T>
class Ajilk_to_Aijkl
{
	const A TA;
      public:

	Ajilk_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim2();
	}

	int get_dim2() const
	{
		return TA.get_dim1();
	}

	int get_dim3() const
	{
		return TA.get_dim4();
	}

	int get_dim4() const
	{
		return TA.get_dim3();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N2,N1,N4,N3);
	}
};


////////////////////////////////////////////
// Define B(i,j,k,l)=A(j,i,l,k)
///////////////////////////////////////////

template < class A, class T>
class Aikjl_to_Aijkl
{
	const A TA;
      public:

	Aikjl_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim1();
	}

	int get_dim2() const
	{
		return TA.get_dim3();
	}

	int get_dim3() const
	{
		return TA.get_dim2();
	}

	int get_dim4() const
	{
		return TA.get_dim4();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N1,N3,N2,N4);
	}
};


////////////////////////////////////////////
// Define B(i,j,k,l)=A(i,k,l,j)
///////////////////////////////////////////

template < class A, class T>
class Aiklj_to_Aijkl
{
	const A TA;
      public:

	Aiklj_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim1();
	}

	int get_dim2() const
	{
		return TA.get_dim3();
	}

	int get_dim3() const
	{
		return TA.get_dim4();
	}

	int get_dim4() const
	{
		return TA.get_dim2();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N1,N3,N4,N2);
	}
};


////////////////////////////////////////////
// Define B(i,j,k,l)=A(j,k,i,l)
///////////////////////////////////////////

template < class A, class T>
class Ajkil_to_Aijkl
{
	const A TA;
      public:

	Ajkil_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim2();
	}

	int get_dim2() const
	{
		return TA.get_dim3();
	}

	int get_dim3() const
	{
		return TA.get_dim1();
	}

	int get_dim4() const
	{
		return TA.get_dim4();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N2,N3,N1,N4);
	}
};


////////////////////////////////////////////
// Define B(i,j,k,l)=A(j,k,l,i)
///////////////////////////////////////////

template < class A, class T>
class Ajkli_to_Aijkl
{
	const A TA;
      public:

	Ajkli_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim2();
	}

	int get_dim2() const
	{
		return TA.get_dim3();
	}

	int get_dim3() const
	{
		return TA.get_dim4();
	}

	int get_dim4() const
	{
		return TA.get_dim1();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N2,N3,N4,N1);
	}
};




////////////////////////////////////////////
// Define B(i,j,k,l)=A(k,i,j,l)
///////////////////////////////////////////

template < class A, class T>
class Akijl_to_Aijkl
{
	const A TA;
      public:

	Akijl_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim3();
	}

	int get_dim2() const
	{
		return TA.get_dim1();
	}

	int get_dim3() const
	{
		return TA.get_dim2();
	}

	int get_dim4() const
	{
		return TA.get_dim4();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N3,N1,N2,N4);
	}
};

////////////////////////////////////////////
// Define B(i,j,k,l)=A(k,i,l,j)
///////////////////////////////////////////

template < class A, class T>
class Akilj_to_Aijkl
{
	const A TA;
      public:

	Akilj_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim3();
	}

	int get_dim2() const
	{
		return TA.get_dim1();
	}

	int get_dim3() const
	{
		return TA.get_dim4();
	}

	int get_dim4() const
	{
		return TA.get_dim2();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N3,N1,N4,N2);
	}
};


////////////////////////////////////////////
// Define B(i,j,k,l)=A(k,j,i,l)
///////////////////////////////////////////

template < class A, class T>
class Akjil_to_Aijkl
{
	const A TA;
      public:

	Akjil_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim3();
	}

	int get_dim2() const
	{
		return TA.get_dim2();
	}

	int get_dim3() const
	{
		return TA.get_dim1();
	}

	int get_dim4() const
	{
		return TA.get_dim4();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N3,N2,N1,N4);
	}
};


////////////////////////////////////////////
// Define B(i,j,k,l)=A(k,j,l,i)
///////////////////////////////////////////

template < class A, class T>
class Akjli_to_Aijkl
{
	const A TA;
      public:

	Akjli_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim3();
	}

	int get_dim2() const
	{
		return TA.get_dim2();
	}

	int get_dim3() const
	{
		return TA.get_dim4();
	}

	int get_dim4() const
	{
		return TA.get_dim1();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N3,N2,N4,N1);
	}
};

///////////////////////////////////////////////////////////////////////
/// Define "MAYOR" permutations
///////////////////////////////////////////////////////////////////////

////////////////////////////////////////////
// Define B(i,j,k,l)=A(k,l,i,j)
///////////////////////////////////////////

template < class A, class T>
class Aklij_to_Aijkl
{
	const A TA;
      public:

	Aklij_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim3();
	}

	int get_dim2() const
	{
		return TA.get_dim4();
	}

	int get_dim3() const
	{
		return TA.get_dim1();
	}

	int get_dim4() const
	{
		return TA.get_dim2();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N3,N4,N1,N2);
	}
};

////////////////////////////////////////////
// Define B(i,j,k,l)=A(l,k,i,j)
///////////////////////////////////////////

template < class A, class T>
class Alkij_to_Aijkl
{
	const A TA;
      public:

	Alkij_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim3();
	}

	int get_dim2() const
	{
		return TA.get_dim4();
	}

	int get_dim3() const
	{
		return TA.get_dim2();
	}

	int get_dim4() const
	{
		return TA.get_dim1();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N4,N3,N1,N2);
	}
};

////////////////////////////////////////////
// Define B(i,j,k,l)=A(l,k,j,i)
///////////////////////////////////////////

template < class A, class T>
class Alkji_to_Aijkl
{
	const A TA;
      public:

	Alkji_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim4();
	}

	int get_dim2() const

	{
		return TA.get_dim3();
	}

	int get_dim3() const
	{
		return TA.get_dim2();
	}

	int get_dim4() const
	{
		return TA.get_dim1();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N4,N3,N2,N1);
	}
};

////////////////////////////////////////////
// Define B(i,j,k,l)=A(k,l,j,i)
///////////////////////////////////////////

template < class A, class T>
class Aklji_to_Aijkl
{
	const A TA;
      public:

	Aklji_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim4();
	}

	int get_dim2() const
	{
		return TA.get_dim3();
	}

	int get_dim3() const
	{
		return TA.get_dim1();
	}

	int get_dim4() const
	{
		return TA.get_dim2();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N3,N4,N2,N1);
	}
};

////////////////////////////////////////////
// Define B(i,j,k,l)=A(l,j,k,i)
///////////////////////////////////////////

template < class A, class T>
class Aljki_to_Aijkl
{
	const A TA;
      public:

	Aljki_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim4();
	}

	int get_dim2() const
	{
		return TA.get_dim2();
	}

	int get_dim3() const
	{
		return TA.get_dim3();
	}

	int get_dim4() const
	{
		return TA.get_dim1();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N4,N2,N3,N1);
	}
};

// Permute index l
//
////////////////////////////////////////////
// Define B(i,j,k,l)=A(i,l,j,k)
///////////////////////////////////////////

template < class A, class T>
class Ailjk_to_Aijkl
{
	const A TA;
      public:

	Ailjk_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim1();
	}

	int get_dim2() const
	{
		return TA.get_dim3();
	}

	int get_dim3() const
	{
		return TA.get_dim4();
	}

	int get_dim4() const
	{
		return TA.get_dim2();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N1,N4,N2,N3);
	}
};

////////////////////////////////////////////
// Define B(i,j,k,l)=A(i,l,j,k)
///////////////////////////////////////////

template < class A, class T>
class Alijk_to_Aijkl
{
	const A TA;
      public:

	Alijk_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim2();
	}

	int get_dim2() const
	{
		return TA.get_dim3();
	}

	int get_dim3() const
	{
		return TA.get_dim4();
	}

	int get_dim4() const
	{
		return TA.get_dim1();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N4,N1,N2,N3);
	}
};




////////////////////////////////////////////
// Define B(i,j,k,l)=A(j,l,i,k)
///////////////////////////////////////////

template < class A, class T>
class Ajlik_to_Aijkl
{
	const A TA;
      public:

	Ajlik_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim3();
	}

	int get_dim2() const
	{
		return TA.get_dim1();
	}

	int get_dim3() const
	{
		return TA.get_dim4();
	}

	int get_dim4() const
	{
		return TA.get_dim2();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N2,N4,N1,N3);
	}
};




////////////////////////////////////////////
// Define B(i,j,k,l)=A(j,l,k,i)
///////////////////////////////////////////

template < class A, class T>
class Ajlki_to_Aijkl
{
	const A TA;
      public:

	Ajlki_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim4();
	}

	int get_dim2() const
	{
		return TA.get_dim1();
	}

	int get_dim3() const
	{
		return TA.get_dim3();
	}

	int get_dim4() const
	{
		return TA.get_dim2();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N2,N4,N3,N1);
	}
};



////////////////////////////////////////////
// Define B(i,j,k,l)=A(l,i,k,j)
///////////////////////////////////////////

template < class A, class T>
class Alikj_to_Aijkl
{
	const A TA;
      public:

	Alikj_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim2();
	}

	int get_dim2() const
	{
		return TA.get_dim4();
	}

	int get_dim3() const
	{
		return TA.get_dim3();
	}

	int get_dim4() const
	{
		return TA.get_dim1();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N4,N1,N3,N2);
	}
};


////////////////////////////////////////////
// Define B(i,j,k,l)=A(l,j,i,k)
///////////////////////////////////////////

template < class A, class T>
class Aljik_to_Aijkl
{
	const A TA;
      public:

	Aljik_to_Aijkl (const A &a): TA(a)
	{
	}

	int get_dim1() const
	{
		return TA.get_dim4();
	}

	int get_dim2() const
	{
		return TA.get_dim2();
	}

	int get_dim3() const
	{
		return TA.get_dim1();
	}

	int get_dim4() const
	{
		return TA.get_dim3();
	}

	T operator() (const int N1, const int N2, const int N3, const int N4) const
	{
		return TA(N4,N2,N1,N3);
	}
};


#endif
