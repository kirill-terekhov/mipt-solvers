#ifndef _GET_TIME_H
#define _GET_TIME_H

#if defined(_WIN32)
#include <windows.h>
static double get_time()
{
	LARGE_INTEGER ticksPerSecond;
	LARGE_INTEGER tick;
	QueryPerformanceFrequency(&ticksPerSecond);
	QueryPerformanceCounter(&tick);
	return (double)tick.QuadPart / (double)ticksPerSecond.QuadPart;
}
#else
#include <sys/time.h>
#include <stddef.h>
static double get_time()
{
	double t1;
	struct timeval time;
	gettimeofday(&(time), NULL);
	t1 = (double)time.tv_sec + (double)time.tv_usec / (1000.0 * 1000.0);
	return t1;
}
#endif

/* requires cxx-11
#include <chrono>
static double get_time()
{
    return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count() / 1.0e6;
}
*/

#endif //_GET_TIME_H
