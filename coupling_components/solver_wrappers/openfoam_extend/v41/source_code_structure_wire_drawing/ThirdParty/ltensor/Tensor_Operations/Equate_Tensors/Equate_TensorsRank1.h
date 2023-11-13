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

#ifndef Equate_TensorsRank1_H
#define Equate_TensorsRank1_H



/* \brief 
*	A method to equal two indexable objects of rank 1. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used
*/
template<class L, class R>
inline void Li_equals_Ri(L & TLeft, const R & TRight)
{
	const int dim1 = TLeft.get_dim1();
 	for(int n1 = 0; n1<dim1;++n1)
	{
		
		TLeft(n1)=TRight(n1);

	}
}


/* \brief 
*	A method to perform += operation between two indexable objects of rank 1. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used
*/

template<class L, class R>
inline void Li_plusequals_Ri(L & TLeft, const R & TRight)
{
	const int dim1 = TLeft.get_dim1();
 	for(int n1 = 0; n1<dim1;++n1)
	{
		TLeft(n1)+=TRight(n1);
	}
}

/* \brief 
*	A method to perform -= operation between two indexable objects of rank 1. 
*	\tparam L The first Object Type
*	\tparam R The second Object Type
*	\param TLeft the object that receives the assignation
*	\param TRight the object the data is extracted
*	\remarks The size of the TLeft object is used
*/
template<class L, class R>
inline void Li_minusequals_Ri(L & TLeft, const R & TRight)
{
	const int dim1 = TLeft.get_dim1();
 	for(int n1 = 0; n1<dim1;++n1)
	{
		TLeft(n1)-=TRight(n1);
	}
}



#endif
