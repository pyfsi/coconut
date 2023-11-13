#ifndef INDEXEDARRAY1
#define INDEXEDARRAY1

template < class A, class T,char i,int indexed>
class IndexedArray1
{
	
      private:
	  A  a;


      public:
      const IndexF<i> &index;

      ~IndexedArray1 ()
	{

	}

	IndexedArray1(const A & arhs, const IndexF<i> &idx): a(arhs),index(idx)

	{

	}

	T operator() (const int N) const
	{
		return a(index(N));
	}

	T & operator()(const int N)
	{
		return a(index(N));
	}

	int get_dim1 () const
	{
		return this->index.get_size();
	}

	

	

};


template < class A, class T,char i>
class IndexedArray1<A,T,i,1>
{
	
      private:
	  A  &a;


      public:
      const IndexF<i> &index;

      ~IndexedArray1 ()
	{

	}

	IndexedArray1(A & arhs, const IndexF<i> &idx): a(arhs),index(idx)

	{

	}

	inline T operator() (const int N) const
	{
		return a(index(N));
	}

	inline T & operator()(const int N)
	{
		return a(index(N));
	}

	int get_dim1 () const
	{
		return this->index.get_size();
	}

	

	

};

#endif
