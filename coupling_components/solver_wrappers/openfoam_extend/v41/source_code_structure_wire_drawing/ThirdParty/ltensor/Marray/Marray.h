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
#ifndef Marray_H
#define Marray_H

#ifndef InheritedBase
	#define InheritedBase TinyArray_base<T,rank>
#endif

#include "../meta/Metaprograms.h"
#include "../base/Array_base.h"

template <char i, int Type>
class Index;

class IndexF;


template < class T, int rank ,class base=InheritedBase > class Marray ;



#include "../Expr/Index.h"
#include "../Expr/IndexF.h"
#include "../Expr/IndexH.h"
#include "../Expr/IndexG.h"
#include "./Marray_rank1.h"
#include "./Marray_rank2.h"
#include "./Marray_rank3.h"
#include "./Marray_rank4.h"



#endif
