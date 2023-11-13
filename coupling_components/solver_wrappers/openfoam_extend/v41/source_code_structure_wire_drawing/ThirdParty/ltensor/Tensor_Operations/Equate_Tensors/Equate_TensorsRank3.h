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
// L(i,j,k) = R(i,j,k)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef Equate_TensorsRank3_H
#define Equate_TensorsRank3_H


//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k) = R(i,j,k)
//////////////////////////////////////////////////////////////////////

/* \brief 
*	A method to equal two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used
*/

template<class L, class R>
inline void Lijk_equals_Rijk(L & TLeft, const R & TRight)
{
	int n1;
	int n2;
	int n3;

 	for(n1 = 0; n1<TLeft.get_dim1();++n1)
	{
		for(n2 = 0; n2<TLeft.get_dim2();++n2)
		{
			for(n3 = 0; n3<TLeft.get_dim3();++n3)
			{
				TLeft(n1,n2,n3) = TRight(n1,n2,n3);
			}
		}
	}
}



/* \brief 
*	A method to perform += operation between two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used
*/
template<class L, class R>
inline void Lijk_plusequals_Rijk(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) += TRight(n1,n2,n3);
			}
		}
	}
}
/* \brief 
*	A method to perform -= operation between two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used
*/

template<class L, class R>
inline void Lijk_minusequals_Rijk(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) -= TRight(n1,n2,n3);
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k) = R(i,k,j)
//////////////////////////////////////////////////////////////////////
/* \brief 
*	A method to equal two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/

template<class L, class R>
inline void Lijk_equals_Rikj(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) = TRight(n1,n3,n2);
			}
		}
	}
}
/* \brief 
*	A method to perform += operation between two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijk_plusequals_Rikj(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) += TRight(n1,n3,n2);
			}
		}
	}
}
/* \brief 
*	A method to perform -= operation between two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijk_minusequals_Rikj(L & TLeft, const R & TRight)
{
	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) -= TRight(n1,n3,n2);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k) = R(j,i,k)
//////////////////////////////////////////////////////////////////////

/* \brief 
*	A method to equal two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijk_equals_Rjik(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) = TRight(n2,n1,n3);
			}
		}
	}
}
/* \brief 
*	A method to perform += operation between two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijk_plusequals_Rjik(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) += TRight(n2,n1,n3);
			}
		}
	}
}
/* \brief 
*	A method to perform -= operation between two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijk_minusequals_Rjik(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) -= TRight(n2,n1,n3);
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k) = R(k,i,j)
//////////////////////////////////////////////////////////////////////

/* \brief 
*	A method to equal two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijk_equals_Rkij(L & TLeft, const R & TRight)
{
	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) = TRight(n3,n1,n2);
			}
		}
	}
}
/* \brief 
*	A method to perform += operation between two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijk_plusequals_Rkij(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) += TRight(n3,n1,n2);
			}
		}
	}
}
/* \brief 
*	A method to perform -= operation between two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijk_minusequals_Rkij(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) -= TRight(n3,n1,n2);
			}
		}
	}
}



//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k) = R(j,k,i)
//////////////////////////////////////////////////////////////////////
/* \brief 
*	A method to equal two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/

template<class L, class R>
inline void Lijk_equals_Rjki(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) = TRight(n2,n3,n1);
			}
		}
	}
}
/* \brief 
*	A method to perform += operation between two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijk_plusequals_Rjki(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) += TRight(n2,n3,n1);
			}
		}
	}
}
/* \brief 
*	A method to perform -= operation between two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijk_mimusequals_Rjki(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) -= TRight(n2,n3,n1);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Definitions of Assignations of tensors such as
// L(i,j,k) = R(k,j,i)
//////////////////////////////////////////////////////////////////////
/* \brief 
*	A method to equal two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/

template<class L, class R>
inline void Lijk_equals_Rkji(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) = TRight(n3,n2,n1);
			}
		}
	}
}
/* \brief 
*	A method to perform += operation between two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijk_plusequals_Rkji(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) += TRight(n3,n2,n1);
			}
		}
	}
}
/* \brief 
*	A method to perform -= operation between two indexable objects of rank 3. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used. NOTE THE PERMUTED INDEXES
*/
template<class L, class R>
inline void Lijk_minusequals_Rkji(L & TLeft, const R & TRight)
{

	const int dim1 = TLeft.get_dim1();
	const int dim2 = TLeft.get_dim2();
	const int dim3 = TLeft.get_dim3();
	int n1;
	int n2;
	int n3;
 	for(n1 = 0; n1<dim1;++n1)
	{
		for(n2 = 0; n2<dim2;++n2)
		{
			for(n3 = 0; n3<dim3;++n3)
			{
				TLeft(n1,n2,n3) -= TRight(n3,n2,n1);
			}
		}
	}
}


#endif
