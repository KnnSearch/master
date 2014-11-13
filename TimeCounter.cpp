#include "TimeCounter.h"

void TimeCounter::StartTime(void) 
{
	QueryPerformanceCounter(&Time); 
}  

double TimeCounter::EndTime()
{     
	LARGE_INTEGER liDiff;     
	LARGE_INTEGER liFreq;      
	QueryPerformanceCounter(&liDiff);      
	liDiff.QuadPart -= Time.QuadPart;     
	liDiff.QuadPart *= 1000000;
	QueryPerformanceFrequency(&liFreq);
	return liDiff.QuadPart/ 1000.0 / liFreq.QuadPart;
}