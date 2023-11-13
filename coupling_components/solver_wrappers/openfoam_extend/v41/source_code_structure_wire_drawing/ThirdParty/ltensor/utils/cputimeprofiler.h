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
#ifndef cputimeprofiler_H
#define cputimeprofiler_H

#ifndef  MSVC

#include <stdio.h>
#include <sys/param.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>



class cputimeprofiler
{
public:
	struct tms t;
	struct tms u;
	long r1;
	long r2;
	double realtime_elapsed;
	double usertime_elapsed;
	double systemtime_elapsed;
	cputimeprofiler():r1(0),r2(0),realtime_elapsed(0),usertime_elapsed(0),systemtime_elapsed(0)
	{ }
	void resettimer()
	{
		realtime_elapsed=0;
		usertime_elapsed=0;
		systemtime_elapsed=0;
	}
	void Start()
	{
		r1 = times(&t);
	}
	void Stop()
	{
		r2 = times(&u);
		usertime_elapsed += (double)(u.tms_utime-t.tms_utime)/(HZ);
		systemtime_elapsed += (double)(u.tms_stime-t.tms_stime)/(HZ);
		realtime_elapsed += (double)(r2-r1)/(HZ);
	}
	void show_timetaken()
	{
		std::cout << "user time= " << usertime_elapsed <<std::endl;
		std::cout << "system time= " <<systemtime_elapsed<<std::endl;
		std::cout << "real time= " << realtime_elapsed <<std::endl;
		std::cout << std::flush;
	}
	double get_usertimetaken()
	{
		return usertime_elapsed;
	}
	double get_realtimetaken()
	{
		return realtime_elapsed;
	}
	double get_systimetaken()
	{
		return systemtime_elapsed;
	}
};

#else
#include <Windows.h>
#include <iostream>
#include <cassert>
class cputimeprofiler
{
public:

	FILETIME CreationTime ;
	FILETIME ExitTime     ;
	FILETIME KernelTime1   ;
	FILETIME KernelTime2  ;
	FILETIME UserTime1    ;
	FILETIME UserTime2   ;
	SYSTEMTIME SystemTime ;
	ULARGE_INTEGER aux1;
	ULARGE_INTEGER aux2;
	ULARGE_INTEGER aux3;
	double userFinal;
	double kernelFinal;
	double totalFinal;


	cputimeprofiler(){resettimer();	  }

	void resettimer()
	{
		userFinal=0.0;
		kernelFinal=0.0;
		totalFinal=0.0;

	}
	void Start()
	{
		//SetThreadPriority(GetCurrentThread(),2);
		int ret = GetProcessTimes( GetCurrentProcess(),
			&CreationTime,
			&ExitTime,
			&KernelTime1,
			&UserTime1 );
		assert(ret!=0);

	}
	void Stop()
	{
		int ret =GetProcessTimes( GetCurrentProcess(),
			&CreationTime,
			&ExitTime,
			&KernelTime2,
			&UserTime2 );

		assert(ret!=0);

		aux1.LowPart = UserTime1.dwLowDateTime;
		aux1.HighPart = UserTime1.dwHighDateTime;
		aux2.LowPart = UserTime2.dwLowDateTime;
		aux2.HighPart = UserTime2.dwHighDateTime;
		aux3.QuadPart = aux2.QuadPart-aux1.QuadPart;

		userFinal += aux3.QuadPart / 10000000.0;


		aux1.LowPart = KernelTime1.dwLowDateTime;
		aux1.HighPart = KernelTime1.dwHighDateTime;
		aux2.LowPart = KernelTime2.dwLowDateTime;
		aux2.HighPart = KernelTime2.dwHighDateTime;
		aux3.QuadPart = aux2.QuadPart-aux1.QuadPart;

		kernelFinal += aux3.QuadPart / 10000000.0;



		totalFinal = userFinal+kernelFinal;



	}
	void show_timetaken()
	{

		std::cout << "user time= " << userFinal <<std::endl;
		std::cout << "kernel time= " <<kernelFinal<<std::endl;
		std::cout << "elapsed time= " << totalFinal<<std::endl;


	}
	double get_usertimetaken()
	{
		return userFinal;
	}
	double get_realtimetaken()
	{
		return totalFinal;
	}
	double get_systimetaken()

	{
		return kernelFinal;
	}
};



#endif



#endif
