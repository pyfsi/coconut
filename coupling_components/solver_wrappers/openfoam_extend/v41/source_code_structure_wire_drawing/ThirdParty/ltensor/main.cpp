#include "LTensor.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include<iomanip>
#include <string>

using namespace std;

int main(void)
  {

        Index<'i'> i;
        Index<'j'> j;
        Marray<double,1> q(2);
		Marray<double,1> p(2);
		Marray<double,2> A(2,2);
		q(i)=A(i,j)*p(j);
		cout<<q;
		//p(i)=q(i)+p(i);

  }
