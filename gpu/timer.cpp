#include <sys/time.h>
#include "timer.h"
using namespace G2G;
using namespace std;

Timer::Timer(void) : started(false) {
}

void Timer::start(void) {
	gettimeofday(&t0, NULL);
	started = true;
}

void Timer::stop(void) {
	gettimeofday(&t1, NULL);
	timersub(&t1, &t0, &res);
	started = false;
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

std::ostream& G2G::operator<<(std::ostream& o, const Timer& t) {
	if (t.getSec() != 0)
		o << t.getSec() << "s. " << t.getMicrosec() << "us.";
	else
		o << t.getMicrosec() << "us.";
	
	return o;
}

