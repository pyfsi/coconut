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
#ifndef TINYARRAY_BASE_H_INCLUDED
#define TINYARRAY_BASE_H_INCLUDED

#include "../meta/Metaprograms.h"
#include <cstdarg>
#include <string>
#include <iostream>
#include <cassert>
//base class for array/tensor implementation

#define CHECK assert(data!=NULL)
#define RELEASE(x) if(x!=NULL) delete[] x;


template <class Type,int rank>
class TinyArray_base{

friend class literator;
    private:



public:

    Type* __restrict__ data;
	int stride[rank];
	unsigned int size[rank];
	int dataCount;

    //View Data

    public:

    class literator{
        private:
		Type* pos;
		int count;
		TinyArray_base<Type,rank> *base;

        public:
			literator(TinyArray_base<Type,rank> *base){
				this->base=base;
				count=0;
                }

			literator(){
				count=0;

			}

            Type& operator*(){
				return *pos;
                }

            void operator++(int){
				if(count<base->dataCount-1){
					pos++;
					count++;
				}
				else{
					pos=NULL;
				}
			}
			void setEnd(){
				pos=NULL;
				count=base->dataCount-1;
			}
			void setBegin(){
				pos=base->data;
			}

			bool operator==(const literator& it){
				return this->pos==it.pos;
			}
			bool operator!=(const literator& it){
				return this->pos!=it.pos;
			}


        };






    /////////////////////////////////////////////////////////////
	// Constructors / Destructors
	/////////////////////////////////////////////////////////////

//TODO overload for rank


    TinyArray_base(long dimension){
        //default CSTYLE
        dataCount=dimension;



        size[0]=dimension;


		if(dataCount==0)
			data=NULL;
		else
			data = new Type[dataCount];

		calcStride();
    }

    TinyArray_base(long dimension1, long dimension2){
        //default CSTYLE
        dataCount=dimension1*dimension2;


        size[0]=dimension1;
        size[1]=dimension2;


		if(dataCount==0)
			data=NULL;
		else
			data = new Type[dataCount];

		calcStride();
    }

    TinyArray_base(long dimension1, long dimension2, long dimension3){
        //default CSTYLE
        dataCount=dimension1*dimension2*dimension3;


        size[0]=dimension1;
        size[1]=dimension2;
        size[2]=dimension3;


		if(dataCount==0)
			data=NULL;
		else
			data = new Type[dataCount];

		calcStride();
    }


    TinyArray_base(long dimension1, long dimension2, long dimension3,long dimension4){
        //default CSTYLE
        dataCount=dimension1*dimension2*dimension3*dimension4;


        size[0]=dimension1;
        size[1]=dimension2;
        size[2]=dimension3;
        size[3]=dimension4;


		if(dataCount==0)
			data=NULL;
		else
			data = new Type[dataCount];

		calcStride();
    }

//    TinyArray_base(long dimensions, ...){
//        //default CSTYLE
//        dataCount=dimensions;
//
//
//		va_list l;
//        va_start(l,dimensions);
//        size[0]=dimensions;
//
//
//
//		for(int i=1;i<rank;i++){
//			size[i]=va_arg(l,long);
//			dataCount*=size[i];
//
//		}
//		/*#ifdef CHECK_OOB
//            assert(dataCount!=0);
//        #endif*/
//
//		if(dataCount==0)
//			data=NULL;
//		else
//			data = new Type[dataCount];
//
//		calcStride();
//    }

	inline void calcStride(){
		stride[0]=1;
        for(int i=1;i<rank;i++){
            stride[i]=stride[i-1]*size[rank-i];
        }

	};
    TinyArray_base(long* dimensions){

		dataCount=dimensions[0];
        #ifdef CHECK_OOB
            assert(dataCount!=0);
        #endif

        size[0]=dimensions[0];
		dataCount=dimensions[0];


		for(int i=1;i<rank;i++){
			size[i]=dimensions[i];
			dataCount*=size[i];
		}
        if(dataCount==0)
		{
		data=NULL;
		} else
		{
		data = new Type[dataCount];
		}
		calcStride();



        }



	//default constructor
    //doesn't initialize anything
    TinyArray_base(){
		data=NULL;
		dataCount=0;
		for(int i=0;i<rank;i++)
            size[i]=0;

    }



	void clear(void)
	{
	    if(data!=NULL)
            delete [] data;
		data=NULL;
		dataCount=0;
		for(int i=0;i<rank;i++)
            size[i]=0;

	}


    virtual ~TinyArray_base(){
        RELEASE(data);

        }

    //what to do with copy constructors
	TinyArray_base(const TinyArray_base<Type,rank> &rterm){

		this->dataCount=rterm->dataCount;

		if (rterm.data!=NULL){
            data = new Type[dataCount];
            copyVector<Type>(this->data,rterm.data,dataCount);
		}else {data=NULL;}

		copyVector<unsigned int>(this->size,rterm->size,rank);
		copyVector<int>(this->stride,rterm->stride,rank);

	}

	literator begin(){

		literator it(this);
		it.setBegin();
		return it;
    }


    literator end(){

		literator it(this);
		it.setEnd();
		return it;
	}


	template <class fType>
    TinyArray_base<Type,rank>& operator=(const TinyArray_base<fType,rank> &rterm){

	     	if(rterm.data!=NULL){

			if (data!=NULL ){
				if( this->dataCount!=rterm.dataCount){
					clear();
					data=new Type[rterm.dataCount];
                    }
                }
			else{
                data = new Type[rterm.dataCount];
                }

			this->dataCount=rterm.dataCount;
			copyVector<unsigned int>(this->size,rterm.size,rank);
			copyVector<Type>(this->data,rterm.data,this->dataCount);
			copyVector<int>(this->stride,rterm.stride,rank);
		}else{
			clear();
		}

       return *this;



    }


    TinyArray_base<Type,rank>& operator=(const TinyArray_base<Type,rank> &rterm){


	        	if(rterm.data!=NULL){

			if (data!=NULL ){
				if( this->dataCount!=rterm.dataCount){
					clear();
					data=new Type[rterm.dataCount];
                    }
                }
			else{
                data = new Type[rterm.dataCount];
                }

			this->dataCount=rterm.dataCount;
			copyVector<unsigned int>(this->size,rterm.size,rank);
			copyVector<Type>(this->data,rterm.data,this->dataCount);
			copyVector<int>(this->stride,rterm.stride,rank);
		}else{
			clear();
		}

       return *this;


    }

void resize(unsigned long* dimensions){



		bool redim=false;

		for(int i=0;i<rank;i++)
			if(size[i]!=dimensions[i])
				redim=true;

		if (!redim) return;

		clear();

        size[0]=dimensions[0];
		dataCount=dimensions[0];
		for(int i=1;i<rank;i++){
			size[i]=dimensions[i];
			dataCount*=size[i];
		}
		if (dataCount>0)
		{
			data = new Type[dataCount];
		} else
		{
			data = NULL;
		}
		calcStride();

	};



    //CAN'T HAVE VIEW OF SAME STORAGE WITH DIFF MATRIX TYPES

    //set Dimensions
    //set Dimension: assert no overflowed

	void setDimension(long d1,...){

        va_list l;
        va_start(l,d1);
        size[0]=d1;
        for (int i=1;i<rank;i++)
            size[i]=va_arg(l,long);
		resize(size);


        }

       inline  Type& getAt(int* n){

			int finalPos=0;
			for (int i=0;i<rank;i++){
				finalPos+= n[i] * stride[rank-1-i];
			}
            return data[finalPos];

        }



///////////////////////////////////////////////////////////////////////////////////////////
			inline Type operator()(const int n1)const {
				return data[n1];
			}
			inline Type& operator()(const int n1) {
				return data[n1];
			}
			inline Type operator()(const int n1,const int n2)const {
				return data[n1*stride[1]+n2*stride[0]];
			}
			inline Type& operator()(const int n1,const int n2) {
				return data[n1*stride[1]+n2*stride[0]];
			}

			inline Type operator()(const int n1,const int n2,const int n3)const {
				return data[n1*stride[2]+n2*stride[1]+n3*stride[0]];
			}

			inline Type& operator()(const int n1,const int n2,const int n3){

				return data[n1*stride[2]+n2*stride[1]+n3*stride[0]];
			}

			inline Type& operator()(const int n1,const int n2,const int n3,const int n4){
				return data[n1*stride[3]+n2*stride[2]+n3*stride[1]+n4*stride[0]];
			}

			inline Type operator()(const int n1,const int n2,const int n3,const int n4)const {
				return data[n1*stride[3]+n2*stride[2]+n3*stride[1]+n4*stride[0]];
			}



///////////////////////////////////////////////////////////////////////////////////////////

    };

#undef RELEASE
#undef CHECK
#endif // TYNYARRAY_BASE_H_INCLUDED
