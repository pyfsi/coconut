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
#ifdef MSVC
#define __restrict__ __restrict
#endif



#include "./base/Array_base.h"
#include "./base/TinyArray_base.h"
#include "./Promote_Types/promote.h"
#include "./Marray/Marray.h"
#include "./Tensor_Operations/Tensor_Operations.h"
#include "./Expr/Exprs.h"


typedef Marray<double,1 > DTensor1;
typedef Marray<double,2 > DTensor2;
typedef Marray<double,3 > DTensor3;
typedef Marray<double,4 > DTensor4;
typedef Marray<int,1 > ITensor1;
typedef Marray<int,2 > ITensor2;
typedef Marray<int,3 > ITensor3;
typedef Marray<int,4 > ITensor4;
typedef Marray<float,1 > FTensor1;
typedef Marray<float,2 > FTensor2;
typedef Marray<float,3 > FTensor3;
typedef Marray<float,4 > FTensor4;


