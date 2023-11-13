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
#ifndef RANGE_H_INCLUDED
#define RANGE_H_INCLUDED
#include "../static/sVector.h"
#include <cmath>

//template only valid for integer types


template <class Type>
class Range{
    public:

    enum Predet{ALL,LHALF,RHALF,EVENS,ODDS,FIXED};
    Type init;
    Type stride;
    Type end;
    int pred;
    int size;
    bool rev;

    Range(Type init,Type stride,Type end){
        this->init=init;
        this->stride=stride;
        this->end=end;
        size=end-init>=0? (end-init+1)/abs(stride): (init-end+1)/abs(stride);
        if(this->end<=this->init)
            rev=true;
        else
            rev=false;
        }
    Range(){
        size=0;

    }

    Range(int pred){
        this->pred=pred;
    }

    Range (int pred,int fixed){
        this->pred=pred;
        this->size=fixed;

        }

    Type* getVector(){

        Type* vec=new Type[size];
        if (rev){
            for (Type i=init;i>=end;i+=stride)
                vec[init-i]=i;
            return vec;


        }else{
            for (Type i=init;i<=end;i+=stride)
                vec[i]=i;
            return vec;
        }


        }


    };


#endif // RANGE_H_INCLUDED
