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

#ifndef ARRAY_BASE_H_INCLUDED
#define ARRAY_BASE_H_INCLUDED
#include "../storage/storage.h"
#include "../meta/Metaprograms.h"
#include "Constants.h"
#include <cassert>
#include <cstdarg>
#include <string>
#include <iostream>
//base class for array/tensor implementation


#define CHECK assert(store!=NULL);
#define RELEASE(x) if(x!=NULL) delete x;



/*! \brief
*   This class represents a slice/view of an Array_base
*
*	The class provides the description of a slice or view of an Array_base
*	It's used to pass view info between Array_bases
*	\tparam Type The type of the data the ArrayBase of which this is a description
*	\tparam rank The rank of the Array_base of which this is a description
*/
template <class Type, int rank>
class SliceDesc{
public:
	/*!	\brief Creates a new SliceDesc

	Builds a new SliceDesc object, allocating the required memory
	\param rnk The rank of the Storage that holds the data that the ArrayBase references

	*/
	SliceDesc(int rnk){
		size=new int[rnk];
		beg=new int[rnk];
		stride=new int[rnk];
		fixed=new bool[rnk];
		dimensions=rnk;
	}

	/*!	\brief Destroys a SliceDesc

	Releases the allocated memory from the heap
	*/

	~SliceDesc(){
		delete []size;
		delete []beg;
		delete []stride;
		delete []fixed;

	}
	/*!	\brief The rank of the ArrayBase*/
	int dimensions;
	/*!	\brief Pointer to an array of size dimensions holding the size of each index*/
	int *size;
	/*!	\brief Pointer to an array of size dimensinons holding the index of the first element*/
	int *beg;
	/*!	\brief Pointer to an array of size dimensinons holding the stride of each index*/
	int *stride;
	/*!	\brief Pointer to an array of size dimensinons that indicates if the dimension is fixed
	*
	*	The dimension is fixed if it was indexed with a constant integer
	*/
	bool *fixed;

	/*!	\brief Pointer to the storage the Array_base uses*/

	GenericStorage<Type> *sto;

};


/*!	\brief This class is used to store information regarding the previous view modifications
*
*	This class is similar in functionality to SliceDesc. Both keep information of the view that an Arraybase
*	represents. But this class is used in an Array_base to store data of previous view modifications that the
*	current Array_base must inherit. It is necesary when a view from a view is created.
*/
class FrameReference{

public:
	/*!	\brief The rank of the Array_base*/
	int dimensions;
	/*!	\brief The constant int for the fixed dimensions*/
	int* value;
	/*!	\brief Pointer to an array of size dimensinons that indicates if the dimension is fixed
	*
	*	The dimension is fixed if it was indexed with a constant integer
	*/
	bool* fixed;

	/*!	\brief Creates a new FrameReference

	Builds a new SliceDesc object
	*/

	FrameReference(){
		value=NULL;
		fixed=NULL;
	}

	/*!	\brief Initializes the FrameReference

	Allocates memory, and configures the object, leaving it ready to use

	\param size The size of the Storage that can be different than the actual Array_base rank
	*/

	void init(int size){
		dimensions=size;
		if(value!=NULL) delete []value;
		if(fixed!=NULL) delete[]fixed;
		value=new int[dimensions];
		fixed=new bool[dimensions];
		for(int i=0;i<dimensions;i++){
			fixed[i]=false;
		}

	}
	/*!	\brief Copy operator overload*/
	FrameReference& operator=(const FrameReference &fr){
		dimensions=fr.dimensions;
		for (int i=0;i<dimensions;i++){
			value[i]=fr.value[i];
			fixed[i]=fr.fixed[i];

		}

		return *this;
	}

	/*!	\brief Destroys the FrameReference

	Deallocates the memory allocated on the constructor
	*/
	~FrameReference(){
		if (value!=NULL) delete[] value;
		if (fixed!=NULL) delete[] fixed;

	}

};
/*!	\brief Array class used as a parent for Marray
*
*	This class offers a view/data architecture. All the Array_bases are simple view of a GenericStorage.
*	It works by reference copy always. All the Marrays that use this as the base, share the same memory
*	space. It is a good choice when working with large Arrays, and memory copy issues are important.
*	It also gives the posibility to have arbitrary dimension ordering in order to have C-Style, FORTRAN-Style
*	and arbitrary-Style arrays.
*	It uses a separated Storage to keep the data, which can be specialized for different needs.
*	Despite the great functionality this class presents, speacial care must be taken when performance
*	is crucial. This is no a performance oriented class. If you want speed check TinyArray_base.
*
*	\tparam Type The type of the data that will be using this class
*	\tparam rank The rank of the tensor this class represents
*
*/

template <class Type,int rank>
class Array_base{

#define ORDER(i) base->store->ordering[i]
	friend class literator;
	friend class Array_base<int,rank>;
	friend class Array_base<double,rank>;
	friend class Array_base<float,rank>;

protected:

	/*!	\brief Pointer to the class that actually stores de data*/
	GenericStorage<Type> * __restrict__ store;

public:
	/*!	\brief Iterator of Array_base data
	*
	*	This class implements an STL like iterator, allowing to iterate an Array_base of
	*	any rank in a linear way.
	*/
	class literator{
	private:
		/*! \brief Array describing the current position in the Array_base. This takes into account the stride*/
		int  pos[rank];
		/*! \brief Array describing the number of iterated elements per dimension*/
		int  count[rank];
		/*! \brief Pointer to the ArrayBase iterated*/
		Array_base<Type,rank>* base;
		/*! \brief Indicates if the iterator has reached the end of the Array_base*/
		bool ended;


	public:
		/*! \brief Creates a new Array_base
		*	\param sto A pointer to the Array_base to iterate
		*/
		literator(Array_base<Type,rank>* sto){
			this->base=sto;
			ended=false;
		}

		/*! \brief Empty Constructor*/

		literator (){
			ended=false;
		}

		/*! \brief Assignation operator */
		literator& operator=(const literator &it){
			this->base=it.base;
			this->ended=false;
			assignLoop<int,const int,rank-1>::VloopAssign(pos,it.pos);
			assignLoop<int,const int,rank-1>::VloopAssign(count,it.count);

			return *this;
		}

		/*! \brief STL compliant operator*/
		Type& operator*(){

			return base->getAt(pos);

		}
		/*! \brief Debug function.
		\return String indicating the current position
		*/
		std::string debugPos(){

			std::string out;
			char temp[15];

			for(int i=rank-1;i>0;i--){
				sprintf(temp,"%d",count[ i ]);
				out+=temp ;
				out+=" ";
			}
			sprintf(temp,"%d",count[0]);
			out+=temp ;
			return out;

		}
		/*! \brief STL compliant operator. Incrementes the position*/
		void operator++(int){

			if(count[ORDER(0)]>=base->size[ORDER(0)]-1){

				for (int i=1;i<rank;i++){

					if (count[ORDER(i)]<base->size[ORDER(i)]-1){ //puedo incrementar
						count[ORDER(i)]++;
						for(int j=i-1;j>=0;j--){
							count[ORDER(j)]=0;
							pos[ORDER(j)]=0;
						}
						pos[ORDER(i)]+=base->stride[ORDER(i)];

						return;
					}

					//end of array

				}
				ended=true;
				pos[0]+=1;



			}
			else{
				count[ORDER(0)]++;
				pos[ORDER(0)]+=base->stride[ORDER(0)];

			}
		}


		/*! \brief Comparison operator*/
		bool operator==(const literator &it){
			return compareVec<int,rank-1>::comp(pos,it.pos);
		}
		/*! \brief Comparison operator*/
		bool operator!=(const literator &it){
			return !compareVec<const int,rank-1>::comp(pos,it.pos);
		}
	private:
		friend class Array_base;

		/*! \brief Positions the iterator at the end*/
		void setEnd(){
			assignLoop<int,long,rank-1>::VloopAssign(this->pos,base->size);
			vecProd<int,rank-1>::prod(this->pos,base->stride);
			ended=true;

		}
		/*! \brief Positions the iterator at the beggining*/
		void setBegin(){
			int b=0;
			assignLoopConst<int,int,rank-1>::loopAssign(pos,b);
			assignLoopConst<int,int,rank-1>::loopAssign(count,b);
			ended=false;

		}

	};

	/*! Indicates the beggining position of each dimension. This accounts the local view nature of the class*/
	long beg[rank];
	/*! Indicates the size of each dimension. This accounts the local view nature of the class*/
	unsigned long size[rank];
	/*! Indicates the stride of each dimension. This accounts the local view nature of the class*/
	int stride[rank];
	/*! Indicates the ordering of each dimension. This accounts the local view nature of the class*/
	bool reversed[rank];
	/*! FrameReference object that holds the inherited view modifications*/
	FrameReference ref;




	/*! \brief
	*	Creates a new Array_base Object
	*
	*	Inits all the variables and arrays, and allocates a new Storage with C-Style ordering
	*
	*	\param dimensions size of the first dimension
	*	\param ... size of the following dimensions
	*
	*	\remarks
	*	This constructor allows the construction of arrays up to rank 4.
	*/

	Array_base(unsigned long dimensions, ...){
		//default CSTYLE
		ref.init(rank);
		va_list l;
		va_start(l,dimensions);
		short ordering[rank];
		short *ran=Range<short>(rank-1,-1,0).getVector();

#ifdef CHECK_OOB
		assert(dimensions>0);
#endif

		copyVector(ordering,ran,rank);
		delete [] ran;

		beg[0]=0;
		size[0]=dimensions;
		stride[0]=1;



		for(int i=1;i<rank;i++){

			beg[i]=0;
			size[i]=va_arg(l,long);

#ifdef CHECK_OOB
			assert(size[i]>0);
#endif

			stride[i]=1;

		}
		//size is ordered
		store=new GenericStorage<Type>(size,rank,ordering);
		ref.dimensions=rank;




	}

	/*! \brief Returns an iterator to the beggining of the Array*/
	literator begin(){
		CHECK;
		literator it(this);
		it.setBegin();
		return it;
	}
	/*! \brief Returns an iterator to the end of the Array*/
	literator end(){
		CHECK;
		literator it(this);
		it.setEnd();
		return it;

	}


	/*! \brief
	*	Creates a new Array_base Object
	*
	*	Inits all the variables and arrays, and allocates a new Storage
	*
	*	\param dimensions Array with the dimensions of each dimension
	*	\param sto An enum describing the type of storage ordering
	*
	*	\remarks
	*	This constructor allows the construction of arrays up to rank 4.
	*/

	Array_base(unsigned long* dimensions,int sto=CSTYLE){


		short* ordering;
		ref.init(rank);
		if (sto==CSTYLE){
			ordering=new short[rank];
			short * ran=Range<short>(rank-1,-1,0).getVector();
			copyVector(ordering,ran,rank);
			delete[] ran;
		}
		else{
			//fortran order
			ordering=new short[rank];
			short * ran=Range<short>(rank-1,-1,0).getVector();
			copyVector(ordering,Range<short>(0,1,rank-1).getVector(),rank);
			delete[] ran;
		}

#ifdef CHECK_OOB
		for (int i=0;i<rank;i++)
			assert(dimensions[i]>0);
#endif

		store=new GenericStorage<Type>(dimensions,rank,ordering);
		delete[] ordering;

		//full matrix view
		for (int i=0;i<rank;i++){
			//begin[i]=store->getBeginAddress(0,i);
			beg[i]=0;
			size[i]=dimensions[i];


			stride[i]=1;


		}
		ref.dimensions=rank;



	}
	//default constructor
	//doesn't initialize anything


	/*! \brief
	*	Creates a new Array_base Object
	*
	*	No initialization or storage allocation performed here. It is left
	*	to the user.
	*
	*	\remarks
	*	The Array must be initialized by copy somewhere
	*/


	Array_base(){
		store=NULL;
		for(int i=0;i<rank;i++){
			size[i]=0;
		}

	}


	/*! \brief
	*	Creates a new Array_base Object by copy
	*
	*	\param rterm Array_base to copy data from
	*
	*	\remarks Both Objects will share de same data
	*/

	Array_base(const Array_base<Type,rank> &rterm){
		copy(rterm);
	}

	/*! \brief
	*	Copies another Array_base
	*
	*	\param rterm Array_base to copy data from
	*
	*	\remarks Both Objects will share de same data
	*/
	void copy(const Array_base<Type,rank> &rterm){
		this->store=rterm.store;
		assignLoop<long,const long,rank-1>::VloopAssign(this->beg,rterm.beg);
		assignLoop<int,const int,rank-1>::VloopAssign(this->stride,rterm.stride);
		assignLoop<unsigned long,const unsigned long,rank-1>::VloopAssign(this->size,rterm.size);
		this->ref.init(rterm.ref.dimensions);
		this->ref=rterm.ref;

		if(store!=NULL)
			store->AddRef();
	}


	/*! \brief
	*	Copies an Array_base
	*
	*	\return A Copy of the current Array_base
	*
	*	\remarks Both Objects will share de same data
	*/
	Array_base& getCopy(){

		Array_base temp;

		if(store!=NULL)
			*temp.store=*store;

		assignLoop<long,long,rank-1>::VloopAssign(temp.beg,this->beg);
		assignLoop<int,int,rank-1>::VloopAssign(temp.stride, this->stride);
		assignLoop<long,long,rank-1>::VloopAssign(temp.size,this->size);
		temp->ref.init(this->ref.dimensions);
		temp->ref=this->ref;

		return temp;

	}



	/*! \brief
	*	Creates a new Array_base Object
	*
	*	Inits all the variables and arrays, and allocates a new Storage
	*
	*	\param dimensions Array with the dimensions of each dimension
	*	\param ordering An array describing the storage ordering
	*
	*	\remarks
	*	This constructor allows the construction of arrays up to rank 4.
	*/

	Array_base(unsigned long* dimensions,short* ordering){

		//custom build
		ref.init(rank);

#ifdef CHECK_OOB
		for (int i=0;i<rank;i++)
			assert(dimensions[i]>0);
#endif

		store=new GenericStorage<Type>(dimensions,rank,ordering);
		for (int i=0;i<rank;i++){
			//begin[i]=store->getBeginAddress(0,i);
			beg[i]=0;
			size[i]=dimensions[i];
			stride[i]=1;
		}
		ref.dimensions=rank;
	}


	/*! \brief
	*	Destroys the Array_base Object
	*
	*	Destroyrs all the variables and arrays, and deallocates Storage if it the last
	*	one pointing it
	*/

	virtual ~Array_base(){

		if(store!=NULL){
			if(store->DelRef()==0)
				delete store;
			store=NULL;
		}
	}



	/////////////////////////////////////////////////////////////
	// Assignation Operations from SCALAR object
	/////////////////////////////////////////////////////////////


	/*! \brief Assignation operator from scalar*/
	template<class U>
	inline const Array_base& operator=(const U &u){

#ifdef CHECK_OOB
		CHECK
#endif
			typename Array_base<Type,rank>::literator it;
		it=begin();
		while(it!=end()){
			*it=u;
			it++;
		}
		return *this;
	}





	/*! \brief Assignation operator from another Array_base
	*
	*	\remarks The copy is by references. All the instances share the same data, EXCEPT:
	*	if the storages are of different Type, then a new storage is created and the data casted
	*	to the new Type. This is the case of the EXCEPTION.
	*/
	template <class fType>
	Array_base<Type,rank>& operator=(const Array_base<fType,rank> rterm){


		assignLoop<unsigned long,const unsigned long, rank-1>::VloopAssign(this->size,rterm.size);
		if(rterm.store!=NULL){
			//different types provides copy
			(*this->store)=(*rterm.store);

			assignLoop<long,const long,rank-1>::VloopAssign(this->beg,rterm.beg);
			assignLoop<int,const int,rank-1>::VloopAssign(this->stride,rterm.stride);
			assignLoop<unsigned long,const unsigned long,rank-1>::VloopAssign(this->size,rterm.size);
			this->ref.init(rterm.ref.dimensions);
			this->ref=rterm.ref;

		}
		return *this;

	}


	/*! \brief Assignation operator from another Array_base
	*
	*	\remarks The copy is by references. No data casting problems here. A storage is only created if
	*	there isn't one initialized yet
	*/
	virtual Array_base<Type,rank>& operator=(const Array_base<Type,rank> &rterm){

		if(store!=NULL){
			if(store->DelRef()==0) delete store;
			store=NULL;
		}
		assignLoop<unsigned long,const unsigned long, rank-1>::VloopAssign(this->size,rterm.size);
		if(rterm.store!=NULL){


			this->store=rterm.store;
			assignLoop<long,const long,rank-1>::VloopAssign(this->beg,rterm.beg);
			assignLoop<int,const int,rank-1>::VloopAssign(this->stride,rterm.stride);
			this->ref.init(rterm.ref.dimensions);
			this->ref=rterm.ref;
			store->AddRef();
		}



		return *this;
	}



	/*! \brief Gets the data from a Slice Description generated by another Array_base
	*
	*	\remarks The two Array_base involved in this operation can have different ranks.
	*	SliceDesc acts as the glue between them
	*/
	template <int rank2>
	Array_base<Type,rank>& operator=(const SliceDesc<Type,rank2> rterm){

		ref.init(rterm.dimensions);


		int curr=0;
		for (int i=0;i<rterm.dimensions;i++){

			if(!rterm.fixed[i]){
				size[curr]=rterm.size[i];
				beg[curr]=rterm.beg[i];
				stride[curr]=rterm.stride[i];
				ref.fixed[i]=false;

				curr++;

			}
			else{
				ref.value[i]=rterm.beg[i];
				ref.fixed[i]=true;

			}

		}
		if(store!=NULL)
			if(store->DelRef()==0) delete store;
		this->store=rterm.sto;
		return *this;



	}




	/*! \brief Resizes if necessary the Storage
	*
	*	\param dimensions An array of size rank containing the new size of each dimension.
	*/

	void resize(unsigned long *dimensions){


		if (store==NULL){
			short *ordering;
			ordering=new short[rank];
			short * ran=Range<short>(rank-1,-1,0).getVector();
			copyVector(ordering,ran,rank);
			delete[] ran;
			delete[] ordering;
			store=new GenericStorage<Type>(dimensions,rank,ordering);
			return;
		}

		bool redim=false;
		for(int i=0;i<rank;i++)

			if (dimensions[i]!=store->dimensions[0]){
				redim=true;
				size[i]=dimensions[i];
			}

			if(redim){
				store->setDimension(dimensions);
			}

	}



	/*! \brief Gets the data at a specific position. This is an alternative to the operator() overload
	*
	*	Returns the data contained in the desired position
	*	\param  n An array of size rank indicating the desired position
	*/
	inline
		Type& getAt(int* n){


#ifdef CHECK_OOB
			CHECK
#endif

				int *idxs=new int[ref.dimensions];
			int curr=0;
			for (int i=0;i<ref.dimensions;i++)
				//mapping function
				if(!ref.fixed[i]){
					idxs[i]=(n[curr]*stride[curr]+beg[curr]);
					curr++;
				}
				else{
					idxs[i]=ref.value[i];
				}
				//store deletes pointer
				return (*store)(idxs);
	}




	/*! \brief Indexing operator for rank 4 tensors	*/
	Type operator()(const int n1,const int n2,const int n3,const int n4)const {

#ifdef CHECK_OOB
		CHECK
#endif
			int* idxs=new int[ref.dimensions];
		int n[4]={n1,n2,n3,n4};
		int curr=0;

		for (int i=0;i<ref.dimensions;i++){
			//mapping function

			if(!ref.fixed[i]){
				idxs[i]=n[curr]*stride[curr]+beg[curr];
				curr++;

			}
			else{
				idxs[i]=ref.value[i];
			}
		}



		return (*store)(idxs);

	}


	/*! \brief Indexing operator for rank 4 tensors	*/
	Type& operator()(const int n1,const int n2,const int n3,const int n4) {

#ifdef CHECK_OOB
		CHECK
#endif
			int* idxs=new int[ref.dimensions];
		int n[4]={n1,n2,n3,n4};
		int curr=0;

		for (int i=0;i<ref.dimensions;i++){
			//mapping function

			if(!ref.fixed[i]){
				idxs[i]=n[curr]*stride[curr]+beg[curr];
				curr++;

			}
			else{
				idxs[i]=ref.value[i];
			}
		}



		return (*store)(idxs);

	}


	/*! \brief Indexing operator for rank 3 tensors	*/
	Type operator()(const int n1,const int n2,const int n3)const {

#ifdef CHECK_OOB
		CHECK
#endif
			int* idxs=new int[ref.dimensions];
		int n[3]={n1,n2,n3};
		int curr=0;

		for (int i=0;i<ref.dimensions;i++){
			//mapping function

			if(!ref.fixed[i]){
				idxs[i]=n[curr]*stride[curr]+beg[curr];
				curr++;

			}
			else{
				idxs[i]=ref.value[i];
			}
		}



		return (*store)(idxs);

	}


	/*! \brief Indexing operator for rank 3 tensors	*/
	Type& operator()(const int n1,const int n2,const int n3) {

#ifdef CHECK_OOB
		CHECK
#endif
			int* idxs=new int[ref.dimensions];
		int n[3]={n1,n2,n3};
		int curr=0;

		for (int i=0;i<ref.dimensions;i++){
			//mapping function

			if(!ref.fixed[i]){
				idxs[i]=n[curr]*stride[curr]+beg[curr];
				curr++;

			}
			else{
				idxs[i]=ref.value[i];
			}
		}



		return (*store)(idxs);

	}



	/*! \brief Indexing operator for rank 2 tensors	*/
	Type operator()(const int n1,const int n2)const {


#ifdef CHECK_OOB
		CHECK
#endif
			int* idxs=new int[ref.dimensions];
		int n[2]={n1,n2};
		int curr=0;

		for (int i=0;i<ref.dimensions;i++){
			//mapping function

			if(!ref.fixed[i]){
				idxs[i]=n[curr]*stride[curr]+beg[curr];
				curr++;

			}
			else{
				idxs[i]=ref.value[i];
			}
		}



		return (*store)(idxs);

	}

	/*! \brief Indexing operator for rank 2 tensors	*/
	Type& operator()(const int n1,const int n2) {

#ifdef CHECK_OOB
		CHECK
#endif
			int* idxs=new int[ref.dimensions];
		int n[2]={n1,n2};
		int curr=0;

		for (int i=0;i<ref.dimensions;i++){
			//mapping function

			if(!ref.fixed[i]){
				idxs[i]=n[curr]*stride[curr]+beg[curr];
				curr++;

			}
			else{
				idxs[i]=ref.value[i];
			}
		}



		return (*store)(idxs);

	}

	/*! \brief Indexing operator for rank 1 tensors	*/
	Type operator()(const int n1)const {

#ifdef CHECK_OOB
		CHECK
#endif
			int* idxs=new int[ref.dimensions];
		int curr=0;

		for (int i=0;i<ref.dimensions;i++){
			//mapping function

			if(!ref.fixed[i]){
				idxs[i]=n1*stride[curr]+beg[curr];
				curr++;

			}
			else{
				idxs[i]=ref.value[i];
			}
		}



		return (*store)(idxs);

	}

	/*! \brief Indexing operator for rank 1 tensors	*/
	Type& operator()(const int n1){

#ifdef CHECK_OOB
		CHECK
#endif
			int* idxs=new int[ref.dimensions];
		int curr=0;

		for (int i=0;i<ref.dimensions;i++){
			//mapping function

			if(!ref.fixed[i]){
				idxs[i]=n1*stride[curr]+beg[curr];
				curr++;

			}
			else{
				idxs[i]=ref.value[i];
			}
		}



		return (*store)(idxs);
	}

	//------------------------CHECKED UP TO HERE----------------------------------------

	/*! \brief Gets a new Slice from the current ArrayBase
	*
	*	Returns a new view description of the current Array_base
	*	\returns A SliceDesc that the assigned class will use to set itself
	*	\param  ranges A vector of range objects indicating the slice
	*/
	SliceDesc<Type, rank>& getSlice(sVector< Range<int>,rank > ranges){

		CHECK;
		//count how many freedom left

		SliceDesc<Type, rank>* __restrict__ desc=new SliceDesc<Type, rank>(ref.dimensions);



		desc->sto=this->store;
		this->store->AddRef();

		for(int i=0;i<ref.dimensions;i++){
			if (ref.fixed[i]==true){
				desc->fixed[i]=true;
				desc->beg[i]=ref.value[i];
				continue;
			}
			switch (ranges(i).pred){
					case Range<int>::FIXED:
						desc->fixed[i]=true;
						desc->beg[i]=ranges(i).size; //fixed number

						break;

					case Range<int>::ALL:

						desc->beg[i]=beg[i];
						desc->size[i]=size[i];
						desc->stride[i]=stride[i];
						desc->fixed[i]=false;
						break;

					default:

						desc->beg[i]=ranges(i).init;
						desc->size[i]=ranges(i).size;
						desc->stride[i]=ranges(i).stride;
						desc->fixed[i]=false;
						break;




			}//end switch
		}//end for

		return *desc;

	}


	/*! \brief Gets a new Slice from the current ArrayBase
	*
	*	Returns a new view description of the current Array_base
	*	\returns A SliceDesc that the assigned class will use to set itself
	*	\param  r1 A Range object indicating the slice of the first dimension
	*	\param  r2 A Range object indicating the slice of the second dimension
	*/
	SliceDesc<Type, rank>& getSlice(Range<int> r1, Range<int> r2){

		sVector <Range<int> ,rank> ranges;
		ranges(0)=r1;
		ranges(1)=r2;
		return getSlice(ranges);
	}


	/*! \brief Gets a new Slice from the current ArrayBase
	*
	*	 Returns a new view description of the current Array_base
	*	\returns A SliceDesc that the assigned class will use to set itself
	*	\param  r1 A Range object indicating the slice of the first dimension
	*/
	SliceDesc<Type, rank>& getSlice(Range<int> r1){
		sVector <Range<int> ,rank> ranges;
		ranges(0)=r1;
		return getSlice(ranges);
	}




};

#undef RELEASE
#undef CHECK

#endif // ARRAY_BASE_H_INCLUDED
