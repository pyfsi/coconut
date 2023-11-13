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
#include "./Expr3.h"
#include "./F_Expr3.h"
#include "./Expr3_Equals.h"
#include "./Expr3_plus_Expr3.h"
#include "./Expr3_minus_Expr3.h"
#include "./Expr3_times_Expr0.h"
#include "./Expr3_divided_Expr0.h"
#include "./minus_Expr3.h"
#include "./plus_Expr3.h"
#include "./Expr3_contract_Expr3.h"
#include "./Expr3_contract_Expr1.h"
#include "./Expr3_contract_Expr2.h"
#include "./Expr2_times_Expr1.h"
