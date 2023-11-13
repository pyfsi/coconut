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
#ifndef METAPROGRAMS_H_INCLUDED
#define METAPROGRAMS_H_INCLUDED
#include "../static/sVector.h"
#include "../base/Array_base.h"
#include <cstdlib>

template <class T,int rank>
class sVector;


template <class Type>
class swapm{
    public:
    static inline void doit(Type &a,Type &b){
        Type c=a;
        a=b;
        b=c;

        }


    };

template <class Type, int I >
class compareVec{
public:
	static inline bool comp(Type arr1[],Type arr2[]){
		if(compareVec<Type,I-1>::comp(arr1,arr2))
			return arr1[I]==arr2[I];
		return false;

	}
};
template <class Type>
class compareVec<Type,0>{
public:
	static inline bool comp(Type arr1[],Type arr2[]){
		return arr1[0]==arr2[0];

	}
};


//very specific watch out
template <class Type, int I>
class vecProd{
public:
	static inline void prod(Type arr1[],Type arr2[]){
		vecProd<Type,I-1>::prod(arr1,arr2);
		arr1[I]=arr1[I]*arr2[I]-1;
	}

};

template <class Type>
class vecProd<Type,0>{
public:
	static inline void prod(Type arr1[],Type arr2[]){
		arr1[0]=arr1[0]*arr2[0];
	}

};



template <class Type1,class Type2,int I>
class assignLoop{
    public:

    static inline void loopAssign(Type1 lterm[],Type2 &rterm){
        assignLoop<Type1,Type2,I-1>::loopAssign(lterm,rterm);
        lterm[I]=rterm(I);
        }

     static inline void VloopAssign(Type1 lterm[],Type2 rterm[]){
        assignLoop<Type1,Type2,I-1>::VloopAssign(lterm,rterm);
        lterm[I]=rterm[I];
        }


    static inline void SloopConcat(Type1 &lterm,Type2 rterm[]){
            assignLoop<Type1,Type2,I-1>::SloopConcat(lterm,rterm);

            char temp[10];
            sprintf(temp,"%d",rterm[I]);
            lterm +=temp;
            lterm += " " ;


        }

    static inline void SVloopAssign(Type1 &lterm,Type2 &rterm){
            assignLoop<Type1,Type2,I-1>::SVloopAssign(lterm,rterm);
            lterm(I)=rterm(I);


        }


    };

template <class Type1,class Type2 >
class assignLoop<Type1,Type2,0>{
    public:
    static inline void loopAssign(Type1 lterm[],Type2 &rterm){
            lterm[0]=rterm(0);
        }

    static inline void VloopAssign(Type1 lterm[],Type2 rterm[]){
            lterm[0]=rterm[0];
        }

    static inline void SloopConcat(Type1 &lterm,Type2 rterm[]){
            char temp[10];

            sprintf(temp,"%d",rterm[0]);

            lterm +=temp;
            lterm +=" " ;
            }

    static inline void SVloopAssign(Type1 &lterm,Type2 &rterm){
        lterm(0)=rterm(0);
        }

    };


template <class Type1,class Type2,int I>
class assignLoopConst{
    public:

    static inline void loopAssign(Type1 lterm[],Type2 &rterm){
        assignLoopConst<Type1,Type2,I-1>::loopAssign(lterm,rterm);
        lterm[I]=rterm;
        }




    };


template <class Type1,class Type2 >
class assignLoopConst<Type1,Type2,0>{
    public:
    static inline void loopAssign(Type1 lterm[],Type2 &rterm){
            lterm[0]=rterm;
        }

    };


template <int I,int rank>
class SizeCounter{
    public:

    static inline void count(sVector<long,rank> &vec, long &acum){

        SizeCounter<I-1,rank>::count(vec,acum);
        acum*=vec(I);

        }

    };

template <int rank>
class SizeCounter<0,rank>{
    public:

    static inline void count(sVector<long,rank> &vec, long &acum){

        acum=vec(0);

        }
    };

    template <int I,int rank>
class StrideCounter{
    public:
    static inline void count(sVector<long,rank> &vec, long acum[],sVector<short,rank> &ord){

        StrideCounter<I-1,rank>::count(vec,acum,ord);
        acum[I]=acum[I-1]*vec(ord.getPos(I));

        }

    };

template <int rank>
class StrideCounter<0,rank>{
    public:

    static inline void count(sVector<long,rank> &vec, long acum[],sVector<short,rank> &ord){

        acum[0]=1;

        }
    };







#endif // METAPROGRAMS_H_INCLUDED
