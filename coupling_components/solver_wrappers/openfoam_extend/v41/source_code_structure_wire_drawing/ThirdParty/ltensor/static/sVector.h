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
#ifndef SVECTOR_H_INCLUDED
#define SVECTOR_H_INCLUDED
#include "../meta/Metaprograms.h"


#include "../base/Range.h"
//template metaprogram




template <class Type, int length>

class sVector{

    private:

    Type data[length];


    public:

    sVector (){


    }



    sVector(Type initVal){
        assignLoopConst<Type,Type,length>::loopAssign(data,initVal);
        }


    sVector(Range <Type> vals){

        if (vals.init<vals.end){
             for (int i=vals.init;i<=vals.end;i+=vals.stride){
                data[i]=i;
                }
        }else{
             for (int i=vals.init;i>=vals.end;i+=vals.stride){
                data[i]=vals.init-i;
                }

            }

        }
    ~sVector(){

        }

     inline Type& operator()(const int &pos){
         return data[pos];
         }

    inline int  getPos(const Type& elem){
        for (int i=0;i<length;i++)
            if (elem==data[i])
                return i;
        return -1;

        }

         //assignation metaprogram

    inline sVector& operator=(Type* data){

        assignLoop<Type,Type,length-1>::VloopAssign(this->data,data);
        return *this;

        }

    inline sVector& operator=(sVector<Type,length>  &data){

        assignLoop<Type,sVector<Type,length>,length-1>::loopAssign(this->data,data);
        return *this;

        }

        template <class U>
    inline sVector& operator=(U &data){

        assignLoopConst<Type,U,length-1>::loopAssign(this->data,data);
        return *this;

        }


    };




#endif // SVECTOR_H_INCLUDED
