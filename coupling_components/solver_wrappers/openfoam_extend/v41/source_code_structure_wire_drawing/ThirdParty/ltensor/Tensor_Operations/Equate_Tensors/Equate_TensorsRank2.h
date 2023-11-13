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
/////////////////////////////////////////////////////////
// Define FUNCTIONS to perform operations of type:
// Expr1 = Expr1
//////////////////////////////////////////////////////////

#ifndef Equate_TensorsRank2_H
#define Equate_TensorsRank2_H


//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j) = R(i,j)
//////////////////////////////////////////////////////////////////////



/* \brief 
*	A method to equal two indexable objects of rank 2. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used
*/

template<class L, class R>
inline void Lij_equals_Rij(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	int n1;
	int n2;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			TLeft(n1,n2) = TRight(n1,n2);
		}
	}
}


/* \brief 
*	A method to perform += operation between two indexable objects of rank 2. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used
*/

template<class L, class R>
inline void Lij_plusequals_Rij(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	int n1;
	int n2;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			TLeft(n1,n2) += TRight(n1,n2);
		}
	}
}

/* \brief 
*	A method to perform -= operation between two indexable objects of rank 2. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used
*/

template<class L, class R>
inline void Lij_minusequals_Rij(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	int n1;
	int n2;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			TLeft(n1,n2) -= TRight(n1,n2);
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j) = R(j,i)
//////////////////////////////////////////////////////////////////////


/* \brief 
*	A method to equal two indexable objects of rank 2. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/

template<class L, class R>
inline void Lij_equals_Rji(L & TLeft, const R & TRight)
{
	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	int n1;
	int n2;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			TLeft(n1,n2) = TRight(n2,n1);
		}
	}
}


/* \brief 
*	A method to perform += operation between two indexable objects of rank 2. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/

template<class L, class R>
inline void Lij_plusequals_Rji(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	int n1;
	int n2;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			TLeft(n1,n2) += TRight(n2,n1);
		}
	}
}


/* \brief 
*	A method to perform -= operation between two indexable objects of rank 2. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/

template<class L, class R>
inline void Lij_minusequals_Rji(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	int n1;
	int n2;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			TLeft(n1,n2) -= TRight(n2,n1);
		}
	}
}

#endif
