#include <cuda_runtime.h>
#include <cstdio>
#include <sys/time.h>
#include <map>
#include <string>
#include "timer.h"
using namespace G2G;
using namespace std;

Timer::Timer(void) : started(false) {
  timerclear(&res);
}

Timer::Timer(const timeval& t) : started(false) {
  res = t;
}

void Timer::start(void) {
#ifdef TIMINGS
	gettimeofday(&t0, NULL);
	started = true;
#endif
}

void Timer::pause(void) {
#ifdef TIMINGS
  gettimeofday(&t1, NULL);
  timeval partial_res;
	timersub(&t1, &t0, &partial_res);
  timeradd(&res, &partial_res, &res);
  timerclear(&t0);
#endif
}

void Timer::stop(void) {
#ifdef TIMINGS
	gettimeofday(&t1, NULL);
	timersub(&t1, &t0, &res);
  timerclear(&t0);
	started = false;
#endif
}

void Timer::start_and_sync(void) {
  sync();
  start();
}

void Timer::stop_and_sync(void) {
  sync();
  stop();
}

void Timer::pause_and_sync(void) {
  sync();
  pause();
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
#ifdef TIMINGS
  #if !CPU_KERNELS
	cudaThreadSynchronize();
  #endif
#endif
}

std::ostream& G2G::operator<<(std::ostream& o, const Timer& t) {
#ifdef TIMINGS
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
#ifdef TIMINGS
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
map<string, Timer> fortran_timers;

extern "C" void timer_start_(const char* timer_name) {
#ifdef TIMINGS
  if (fortran_timers.find(timer_name) == fortran_timers.end()) fortran_timers[timer_name] = Timer();
  Timer::sync();
	fortran_timers[timer_name].start();
#endif
}

extern "C" void timer_stop_(const char* timer_name) {
#ifdef TIMINGS
	Timer::sync();
  if (fortran_timers.find(timer_name) == fortran_timers.end()) cout << "no existe timer!" << endl;
	fortran_timers[timer_name].stop();
  cout << "TIMER [" << timer_name << "]: " << fortran_timers[timer_name] << endl;
#endif
}

extern "C" void timer_pause_(const char* timer_name) {
#ifdef TIMINGS
	Timer::sync();
	if (fortran_timers.find(timer_name) == fortran_timers.end()) cout << "no existe timer!" << endl;
  fortran_timers[timer_name].pause();
  cout << "TIMER [" << timer_name << "]: " << fortran_timers[timer_name] << "(so far)" << endl;
#endif
}
