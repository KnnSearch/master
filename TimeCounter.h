#ifndef _TIMECOUNTER_H
#define _TIMECOUNTER_H

#include <windows.h> 

#include <stdlib.h> 
#include <iostream>

class TimeCounter
{
	LARGE_INTEGER Time;
public:
	void StartTime(void);
	double EndTime();
};

#endif
