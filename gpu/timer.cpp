#include <cuda_runtime.h>
#include <cstdio>
#include <sys/time.h>
#include "timer.h"
using namespace G2G;
using namespace std;

Timer::Timer(void) : started(false) {
}

void Timer::start(void) {
#ifdef DO_TIMINGS
	gettimeofday(&t0, NULL);
	started = true;
#endif
}

void Timer::stop(void) {
#ifdef DO_TIMINGS
	gettimeofday(&t1, NULL);
	timersub(&t1, &t0, &res);
	started = false;
#endif
}

bool Timer::isStarted(void) const {
	return started;	
}

unsigned long Timer::getMicrosec(void) const {
	return res.tv_usec;
}

unsigned long Timer::getSec(void) const {
	return res.tv_sec;
}

bool Timer::operator<(const Timer& other) const {
	return (res.tv_sec < other.res.tv_sec ||
					(res.tv_sec == other.res.tv_sec && res.tv_usec < other.res.tv_usec));
}

void Timer::sync(void) {
#ifdef DO_TIMINGS
	cudaThreadSynchronize();
#endif
}

std::ostream& G2G::operator<<(std::ostream& o, const Timer& t) {
#ifdef DO_TIMINGS
	if (t.getSec() != 0)
		o << t.getSec() << "s. " << t.getMicrosec() << "us.";
	else
		o << t.getMicrosec() << "us.";
	
#else
	o << "[TIMINGS NOT ENABLED]";
#endif
	return o;
}

void Timer::print(void) {
#ifdef DO_TIMINGS
	if (getSec() != 0)
		printf("%lus. %luus.", getSec(), getMicrosec());
	else
		printf("%luus.", getMicrosec());
#else
	printf("%s","[TIMINGS NOT ENABLED]");
#endif
}

/**** to be used by fortran ****/
Timer global_timer;
extern "C" void timer_start_(void) {
#ifdef DO_TIMINGS
	global_timer.start();
#endif
}

extern "C" void timer_stop_(const char* what) {
#ifdef DO_TIMINGS
	Timer::sync();
	global_timer.stop();
	printf("TIMER (%s): ", what);
	global_timer.print();
	printf("\n");
#endif
}
