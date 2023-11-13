#ifndef ALGORITHM_H
#define ALGORITHM_H

template <typename T>
T Sign(T t)
{
    if( t == 0 )
        return T(0);
    else
        return ( t < 0 ? T(-1) : T(1) );
}


#endif
