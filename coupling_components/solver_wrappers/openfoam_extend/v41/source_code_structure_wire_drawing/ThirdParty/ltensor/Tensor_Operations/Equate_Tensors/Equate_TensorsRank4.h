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
//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k,l) = R(i,j,k,l)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef Equate_TensorsRank4_H
#define Equate_TensorsRank4_H


//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k,l) = R(i,j,k,l)
//////////////////////////////////////////////////////////////////////
/* \brief 
*	A method to equal two indexable objects of rank 4. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used
*/

template<class L, class R>
inline void Lijkl_equals_Rijkl(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	const int dim4 = TLeft.get_dim4();
	int n1;
	int n2;
	int n3;
	int n4;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				for(n4 = 0; n4<dim4;++n4)
				{
					TLeft(n1,n2,n3,n4) = TRight(n1,n2,n3,n4);
				}
			}
		}
	}
}

template<class L, class R>
inline void Lijkl_plusequals_Rijkl(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	const int dim4 = TLeft.get_dim4();
	int n1;
	int n2;
	int n3;
	int n4;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				for(n4 = 0; n4<dim4;++n4)
				{
					TLeft(n1,n2,n3,n4) += TRight(n1,n2,n3,n4);
				}
			}
		}
	}
}


template<class L, class R>
inline void Lijkl_minusequals_Rijkl(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	const int dim4 = TLeft.get_dim4();
	int n1;
	int n2;
	int n3;
	int n4;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				for(n4 = 0; n4<dim4;++n4)
				{
					TLeft(n1,n2,n3,n4) -= TRight(n1,n2,n3,n4);
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k,l) = R(i,j,l,k)
//////////////////////////////////////////////////////////////////////

/* \brief 
*	A method to equal two indexable objects of rank 4. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/

template<class L, class R>
inline void Lijkl_equals_Rijlk(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	const int dim4 = TLeft.get_dim4();
	int n1;
	int n2;
	int n3;
	int n4;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				for(n4 = 0; n4<dim4;++n4)
				{
					TLeft(n1,n2,n3,n4) = TRight(n1,n2,n4,n3);
				}
			}
		}
	}
}
/* \brief 
*	A method to perform += operation between two indexable objects of rank 4. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used
*/
template<class L, class R>
inline void Lijkl_plusequals_Rijlk(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	const int dim4 = TLeft.get_dim4();
	int n1;
	int n2;
	int n3;
	int n4;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				for(n4 = 0; n4<dim4;++n4)
				{
					TLeft(n1,n2,n3,n4) += TRight(n1,n2,n4,n3);
				}
			}
		}
	}
}
/* \brief 
*	A method to perform -= operation between two indexable objects of rank 4. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used
*/
template<class L, class R>
inline void Lijkl_minusequals_Rijlk(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	const int dim4 = TLeft.get_dim4();
	int n1;
	int n2;
	int n3;
	int n4;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				for(n4 = 0; n4<dim4;++n4)
				{
					TLeft(n1,n2,n3,n4) -= TRight(n1,n2,n4,n3);
				}
			}
		}
	}
}



//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k,l) = R(j,i,l,k)
//////////////////////////////////////////////////////////////////////

/* \brief 
*	A method to equal two indexable objects of rank 4. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijkl_equals_Rjilk(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	const int dim4 = TLeft.get_dim4();
	int n1;
	int n2;
	int n3;
	int n4;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				for(n4 = 0; n4<dim4;++n4)
				{
					TLeft(n1,n2,n3,n4) = TRight(n2,n1,n4,n3);
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k,l) = R(j,i,k,l)
//////////////////////////////////////////////////////////////////////

/* \brief 
*	A method to equal two indexable objects of rank 4. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijkl_equals_Rjikl(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	const int dim4 = TLeft.get_dim4();
	int n1;
	int n2;
	int n3;
	int n4;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				for(n4 = 0; n4<dim4;++n4)
				{
					TLeft(n1,n2,n3,n4) = TRight(n2,n1,n3,n4);
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k,l) = R(k,l,i,j)
//////////////////////////////////////////////////////////////////////

/* \brief 
*	A method to equal two indexable objects of rank 4. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijkl_equals_Rklij(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	const int dim4 = TLeft.get_dim4();
	int n1;
	int n2;
	int n3;
	int n4;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				for(n4 = 0; n4<dim4;++n4)
				{
					TLeft(n1,n2,n3,n4) = TRight(n3,n4,n1,n2);
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k,l) = R(k,l,j,i)
//////////////////////////////////////////////////////////////////////

/* \brief 
*	A method to equal two indexable objects of rank 4. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijkl_equals_Rklji(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	const int dim4 = TLeft.get_dim4();
	int n1;
	int n2;
	int n3;
	int n4;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				for(n4 = 0; n4<dim4;++n4)
				{
					TLeft(n1,n2,n3,n4) = TRight(n3,n4,n2,n1);
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k,l) = R(l,k,j,i)
//////////////////////////////////////////////////////////////////////

/* \brief 
*	A method to equal two indexable objects of rank 4. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijkl_equals_Rlkji(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	const int dim4 = TLeft.get_dim4();
	int n1;
	int n2;
	int n3;
	int n4;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				for(n4 = 0; n4<dim4;++n4)
				{
					TLeft(n1,n2,n3,n4) = TRight(n4,n3,n2,n1);
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k,l) = R(l,k,i,j)
//////////////////////////////////////////////////////////////////////

/* \brief 
*	A method to equal two indexable objects of rank 4. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijkl_equals_Rlkij(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	const int dim4 = TLeft.get_dim4();
	int n1;
	int n2;
	int n3;
	int n4;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				for(n4 = 0; n4<dim4;++n4)
				{
					TLeft(n1,n2,n3,n4) = TRight(n4,n3,n1,n2);
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k,l) = R(l,j,k,i)
//////////////////////////////////////////////////////////////////////

/* \brief 
*	A method to equal two indexable objects of rank 4. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijkl_equals_Rljki(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	const int dim4 = TLeft.get_dim4();
	int n1;
	int n2;
	int n3;
	int n4;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				for(n4 = 0; n4<dim4;++n4)
				{
					TLeft(n1,n2,n3,n4) = TRight(n4,n2,n3,n1);
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k,l) = R(i,k,j,l)
//////////////////////////////////////////////////////////////////////

/* \brief 
*	A method to equal two indexable objects of rank 4. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijkl_equals_Rikjl(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	const int dim4 = TLeft.get_dim4();
	int n1;
	int n2;
	int n3;
	int n4;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				for(n4 = 0; n4<dim4;++n4)
				{
					TLeft(n1,n2,n3,n4) = TRight(n1,n3,n2,n4);
				}
			}
		}
	}
}



#endif
