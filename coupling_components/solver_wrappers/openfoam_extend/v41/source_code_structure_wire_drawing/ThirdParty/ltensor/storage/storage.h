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
#ifndef STORAGE_H_INCLUDED
#define STORAGE_H_INCLUDED
#include "../static/sVector.h"
#include "../meta/Metaprograms.h"

#ifndef NULL
    #define NULL 0
#endif
#define RELEASE(x) if(x!=NULL) delete[] x;
//storage class
//used to manage memory blocks
//and storage configurations

//template Metaprogram

template <class T>
inline void copyVector(T* v1,T* v2,int n){
    for(int i=0;i<n;i++)
        v1[i]=v2[i];
    }

template <class T1,class T2>
inline void copyVector(T1* v1,T2* v2,int n){
    for(int i=0;i<n;i++)
        v1[i]=static_cast<T1>(v2[i]);
    }

template <class T>

int getPos(T* v,T val,int n){
    for(int i=0;i<n;i++){
            if(v[i]==val) return i;
        }
    return -1;

    }

template <class T1,class T2>

int getPos(T1* v,T2 val,int n){
    for(int i=0;i<n;i++){
            if(v[i]==val) return i;
        }
    return -1;

    }



template <class Type>
class GenericStorage{

public:

    //to handle different memory mappings
    short rank;
    short* ordering;
    unsigned long* dimensions;
    long* stride;

    Type* __restrict__ data ;

    int size;

    int refCount;

    void inline initVectors(){

        ordering=new short[rank];
        dimensions=new unsigned long[rank];
        stride=new long[rank];
        }

    void calcDimensions(){
        size=1;
        for (int i=0;i<rank;i++)
            size*=dimensions[i];
        }

    public:

    GenericStorage(unsigned long* dimensions,int n ){
		data=NULL;

        rank=n;
        initVectors();
        refCount=1;
        copyVector(this->dimensions,dimensions,rank);
        allocate();
        }

    GenericStorage(unsigned long* dimensions,int n,short* ordering=NULL){
		data=NULL;
        rank=n;
        initVectors();
        copyVector(this->dimensions,dimensions,rank);

        if (ordering==NULL)
            //default to cStyleOrdering
            for (int i=0;i<rank;i++)
                this->ordering[i]=rank-1-i;
        else
            copyVector(this->ordering,ordering,rank);

        this->refCount=1;
        allocate();
        }





    GenericStorage(){
        initVectors();
        data=NULL;

        this->refCount=1;
        }

        ~GenericStorage(){
            RELEASE(data);
			data=NULL;
            delete[] ordering;
            delete[] dimensions;
            delete[] stride;

        }

    void setDimension(unsigned long* dimensions){
        RELEASE(data);
		data=NULL;
        copyVector(this->dimensions,dimensions,rank);
        allocate();
        }

    GenericStorage& operator=(GenericStorage<Type> &sto){

         for (int i=0;i<rank;i++){
            ordering[i]=sto.ordering[i];
            dimensions[i]=sto.dimensions[i];
			stride[i]=sto.stride[i];
            }
       if(data!=NULL){
				if(size!=sto.size){
				delete[] data;
				data=new Type[sto.size];
				}
           }
		size=sto.size;
        for(int i=0;i<size;i++)
            data[i]=sto.data[i];
        return *this;

        }
    template <class Type2>
    GenericStorage& operator=(GenericStorage<Type2> &sto){


         for (int i=0;i<rank;i++){
            ordering[i]=sto.ordering[i];
            dimensions[i]=sto.dimensions[i];
			stride[i]=sto.stride[i];
            }
		 if(data!=NULL){
				if(size!=sto.size){
				delete[] data;
				data=new Type[sto.size];
				}
           }
		 size=sto.size;


        for(int i=0;i<size;i++)
            data[i]=sto.data[i];
        return *this;

        }

    long getBeginAddress(int offset,int dim){

        long gsize=0;
        gsize= (dim==0?0: stride(getPos(ordering,dim-1,rank) )) ;
        return gsize+offset;
        }


    void calcStride(){

        stride[ordering[0]]=1;
        for(int i=1;i<rank;i++){
            stride[ordering[i]]=stride[ordering[i-1]]*dimensions[ordering[i-1]];

        }
     }



    int AddRef(){
        return ++this->refCount;
        }

    int DelRef(){
        //deallocate
        return --this->refCount;
        }


    void allocate(){

        //i don't know the rank so can't metaprogram here

        long gsize=dimensions[0];

        for(int i=1;i<rank;i++)
                 gsize*=dimensions[i];

        data=new Type[gsize];
        calcStride();
        calcDimensions();

        }


    Type& operator()(int* idx){

		int pos=0;
        for (int i=0;i<rank-1;i++){
			pos+= idx[i]*stride[i];
        }
        int temp=pos + idx[rank-1];
        delete[] idx;

            return data[temp];

        }

	Type operator()(int* idx)const {

		int pos=0;
        for (int i=0;i<rank-1;i++){
			pos+= idx[i]*stride[i];
        }
        int temp=pos + idx[rank-1];
        delete[] idx;
            return data[temp];

        }





    };




#undef RELEASE
#endif // STORAGE_H_INCLUDED


