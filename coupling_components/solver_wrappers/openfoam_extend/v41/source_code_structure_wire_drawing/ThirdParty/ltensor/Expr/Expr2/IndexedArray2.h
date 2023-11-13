#ifndef INDEXEDARRAY2
#define INDEXEDARRAY2

template < class A, class T,char i , char j,int indexed>
class IndexedArray2
{


	private:

	A  &a;

      public:

    const IndexF<i> &index1;
    const IndexF<j> &index2;
      ~IndexedArray2 ()
	{
	}

	  IndexedArray2 (const A & arhs,const IndexF<i> &idx1,const IndexF<j> &idx2): a(arhs), index1(idx1), index2(idx2)
	{

	}


	inline T operator() (const int N1, const int N2) const
	{

		return a(index1(N1),index2(N2));
	}

	inline T & operator()(const int N1, const int N2)
	{

		return a(index1(N1),index2(N2));
	}

	int get_dim1 () const
	{
		return this->index1.get_size();
	}

	int get_dim2 () const
	{
		return this->index2.get_size();
	}

	
};




template < class A, class T,char i, char j>
class IndexedArray2<A,T,i,j,2>
{


	private:

	A  &a;

      public:

    const IndexF<i> &index2;
      ~IndexedArray2 ()
	{
	}

	  IndexedArray2 (const A & arhs,const IndexF<i> &idx2): a(arhs), index2(idx2)
	{

	}


	T operator() (const int N1, const int N2) const
	{

		return a(N1,index2(N2));
	}

	T & operator()(const int N1, const int N2)
	{

		return a(N1,index2(N2));
	}

	int get_dim1 () const
	{
		return this->A.get_dim1();
	}

	int get_dim2 () const
	{
		return this->index2.get_size();
	}

	
};




template < class A, class T,char i, char j>
class IndexedArray2<A,T,i,j,1>
{


	private:

	A  &a;

      public:

    const IndexF<i> &index1;
      ~IndexedArray2 ()
	{
	}

	  IndexedArray2 (const A & arhs,const IndexF<i> &idx1): a(arhs), index1(idx1)
	{

	}


	T operator() (const int N1, const int N2) const
	{

		return a(index1(N1),N2);
	}

	T & operator()(const int N1, const int N2)
	{

		return a(index1(N1),N2);
	}

	int get_dim1 () const
	{
		return this->index2.get_size();
	}

	int get_dim2 () const
	{
		return this->A.get_dim1();
	}

	
};



template < class A, class T,char i, char j>
class IndexedArray2<A,T,i,j,12>
{


	private:

	A  a;

      public:

    const IndexF<i> &index1;
    const IndexF<j> &index2;
      ~IndexedArray2 ()
	{
	}

	  IndexedArray2 (const A & arhs,const IndexF<i> &idx1,const IndexF<j> &idx2): a(arhs), index1(idx1), index2(idx2)
	{

	}


	T operator() (const int N1, const int N2) const
	{

		return a(index1(N1),index2(N2));
	}

	T & operator()(const int N1, const int N2)
	{

		return a(index1(N1),index2(N2));
	}

	int get_dim1 () const
	{
		return this->index1.get_size();
	}

	int get_dim2 () const
	{
		return this->index2.get_size();
	}


};





#endif

