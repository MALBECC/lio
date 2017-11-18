#ifndef __TIMER_H__
#define __TIMER_H__

#include <iostream>

namespace G2G {

class Timer {
 public:
  Timer(void);
  Timer(const timespec& t);

  void start(void);
  void stop(void);
  void pause(void);
  void start_and_sync(void);
  void stop_and_sync(void);
  void pause_and_sync(void);

  unsigned long getMicrosec(void) const;
  unsigned long getSec(void) const;
  double getTotal(void) const;

  bool isStarted(void) const;

  friend std::ostream& operator<<(std::ostream& o, const Timer& t);
  // to compare stopped timers
  bool operator<(const Timer& other) const;
  void print(void);
  static void sync(void);

 private:
  timespec t0, t1, res;
  bool started;
};

std::ostream& operator<<(std::ostream& o, const Timer& t);
}

#endif
