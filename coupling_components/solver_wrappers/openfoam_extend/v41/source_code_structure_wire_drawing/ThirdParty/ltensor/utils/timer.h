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
// very simple, high-precision (at least in theory) timer for timing API calls

#ifndef MSVC

#include <sys/time.h>

class HRTimer {

    public:

    struct timeval beg;
    struct timeval end;
    struct timezone tz;

    void Start(void) {
        gettimeofday(&beg, &tz);
    };
    void Stop(void) {
        gettimeofday(&end, &tz);
    };
    double GetDurationInSecs(void)
    {

        double duration = (end.tv_sec+end.tv_usec*1e-6)- (beg.tv_sec+beg.tv_usec*1e-6);
        return duration;
    }
};
#else

#include <Windows.h>
// very simple, high-precision (at least in theory) timer for timing API calls
struct HRTimer {

	LARGE_INTEGER freq;
	void Init(){
		QueryPerformanceFrequency(&freq);
	}
    void Start(void) {
        QueryPerformanceCounter(&mTimeStart);
    };
    void Stop(void) {
        QueryPerformanceCounter(&mTimeStop);
    };
    double GetDurationInSecs(void)
    {      
        
        double duration = (double)(mTimeStop.QuadPart-mTimeStart.QuadPart)/(double)freq.QuadPart;
        return duration;
    }

    LARGE_INTEGER mTimeStart;
    LARGE_INTEGER mTimeStop;
};



#endif

