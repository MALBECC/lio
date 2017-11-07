#include <cstdio>
#include <time.h>
#include <sys/time.h>
#include <map>
#include <set>
#include <string>
#include "init.h"
#include "timer.h"

// Only include it to sync timings in cuda threads.
#include "cuda_includes.h"

using namespace G2G;
using namespace std;

/* adapt macros to timespec */
# define timerspecclear(tvp)  ((tvp)->tv_sec = (tvp)->tv_nsec = 0)
# define timerspecadd(a, b, result)                 \
  do {                        \
    (result)->tv_sec = (a)->tv_sec + (b)->tv_sec;           \
    (result)->tv_nsec = (a)->tv_nsec + (b)->tv_nsec;            \
    if ((result)->tv_nsec >= 1000000000)               \
      {                       \
  ++(result)->tv_sec;                 \
  (result)->tv_nsec -= 1000000000;               \
      }                       \
  } while (0)
# define timerspecsub(a, b, result)                 \
  do {                        \
    (result)->tv_sec = (a)->tv_sec - (b)->tv_sec;           \
    (result)->tv_nsec = (a)->tv_nsec - (b)->tv_nsec;            \
    if ((result)->tv_nsec < 0) {                \
      --(result)->tv_sec;                 \
      (result)->tv_nsec += 1000000000;               \
    }                       \
  } while (0)


Timer::Timer(void) : started(false)
{
   timerspecclear(&res);
}

Timer::Timer(const timespec& t) : started(false) 
{
   res = t;
}

void Timer::start(void) 
{
   started = true;
   clock_gettime(CLOCK_MONOTONIC, &t0);
}

void Timer::pause(void)
{
   clock_gettime(CLOCK_MONOTONIC, &t1);
   timespec partial_res;
   timerspecsub(&t1, &t0, &partial_res);
   timerspecadd(&res, &partial_res, &res);
   timerspecclear(&t0);
}

void Timer::stop(void)
{
   clock_gettime(CLOCK_MONOTONIC, &t1);
   timerspecsub(&t1, &t0, &res);
   timerspecclear(&t0);
   started = false;
}

void Timer::start_and_sync(void)
{
   sync();
   start();
}

void Timer::stop_and_sync(void)
{
   sync();
   stop();
}

void Timer::pause_and_sync(void) 
{
   sync();
   pause();
}

bool Timer::isStarted(void) const
{
   return started;
}

unsigned long Timer::getMicrosec(void) const 
{
   return res.tv_nsec / 1000;
}

unsigned long Timer::getSec(void) const
{
   return res.tv_sec;
}

double Timer::getTotal(void) const
{
   return res.tv_nsec + res.tv_sec * 1000.0 * 1000.0 * 1000.0;
}

bool Timer::operator<(const Timer& other) const 
{
   return (res.tv_sec < other.res.tv_sec || (res.tv_sec == other.res.tv_sec && res.tv_nsec < other.res.tv_nsec));
}

void Timer::sync(void) 
{
#if GPU_KERNELS
   cudaThreadSynchronize();
#endif
}

std::ostream& G2G::operator<<(std::ostream& o, const Timer& t) 
{
   if (t.getSec() != 0)
   {
      o << t.getSec() << "s. " << t.getMicrosec() << "us.";
   } else {
      o << t.getMicrosec() << "us.";
   }
   return o;
}

void Timer::print(void)
{
   if (getSec() != 0)
   {
      printf("%lus. %luus.", getSec(), getMicrosec());
   } else {
      printf("%luus.", getMicrosec());
   }
}

/**** to be used by fortran ****/
Timer global_timer;
string current_timer;
map<string, set<string> > timer_children;
map<string, string> timer_parents;
map<string, Timer*> all_timers;
map<string, Timer*> top_timers;
map<string, Timer> fortran_timers;

// One-time Fortran timer calls - these time sections and report timings immediately
extern "C" void g2g_timer_start_(const char* timer_name, unsigned int length_arg)
{
   if (G2G::timer_single)
   {
      string tname(timer_name, length_arg);
      tname.append("\0");

      if (fortran_timers.find(tname) == fortran_timers.end())
      {
         fortran_timers[tname] = Timer();
      }
      Timer::sync();
      fortran_timers[tname].start();
   }
}

extern "C" void g2g_timer_stop_(const char* timer_name, unsigned int length_arg)
{
   if (G2G::timer_single)
   {
      string tname(timer_name, length_arg);
      tname.append("\0");

      Timer::sync();
      if (fortran_timers.find(tname) == fortran_timers.end())
      {
         cout << "no existe timer! (" << tname << ")" << endl;
      }
      fortran_timers[tname].stop();

      cout << "TIMER [" << tname << "]: " << fortran_timers[tname] << endl;
   }
}

extern "C" void g2g_timer_pause_(const char* timer_name, unsigned int length_arg)
{
   if (G2G::timer_single)
   {
      string tname(timer_name, length_arg);
      tname.append("\0");

      Timer::sync();
      if (fortran_timers.find(tname) == fortran_timers.end())
      {
         cout << "no existe timer! (" << tname << ")" << endl;
      }
      fortran_timers[tname].pause();

      cout << "TIMER [" << tname << "]: " << fortran_timers[tname] << "(so far)" << endl;
   }
}

// Summary timer calls - these store timings until g2g_timer_summary is called, when a summary
// of all timings up to that point is given in a tree-sorted display
extern "C" void g2g_timer_sum_start_(const char* timer_name, unsigned int length_arg)
{
   if (G2G::timer_sum)
   {
      string tname(timer_name,length_arg);
      tname.append("\0");
      if (timer_children.find(tname) == timer_children.end())
      {
         timer_children[tname] = set<string>();
      }
      if (current_timer.length() == 0)
      {
         if (top_timers.find(tname) == top_timers.end())
         {
            top_timers[tname] = new Timer();
         }
         Timer::sync();
         top_timers[tname]->start();
         if (all_timers.find(tname) == all_timers.end())
         {
            all_timers[tname] = top_timers[tname];
         }
        current_timer = tname;
      } else {
         if ((timer_children[current_timer]).find(tname) == (timer_children[current_timer]).end())
         {
            timer_children[current_timer].insert(tname);
            all_timers[tname] = new Timer();
            timer_parents[tname] = current_timer;
         }
         Timer::sync();
         all_timers[tname]->start();
         current_timer = tname;
      }
   }
}

extern "C" void g2g_timer_sum_stop_(const char* timer_name, unsigned int length_arg)
{
   if (G2G::timer_sum)
   {
      string tname(timer_name, length_arg);
      tname.append("\0");
      Timer::sync();
      if (current_timer.compare(tname) != 0)
      {
         cout << "Error: not the current timer: (" << tname << "," << current_timer << ")" << endl;
      } else {
         all_timers[current_timer]->stop();
         if (timer_parents.find(current_timer) == timer_parents.end())
         {
            current_timer = "";
         } else {
            current_timer = timer_parents[current_timer];
         }
      }
   }
}

extern "C" void g2g_timer_sum_pause_(const char* timer_name, unsigned int length_arg)
{
   if (G2G::timer_sum)
   {
      string tname(timer_name, length_arg);
      tname.append("\0");
      Timer::sync();
      if (current_timer.compare(tname) != 0)
      {
         cout << "Error: not the current timer: (" << tname << ")" << endl;
      } else {
         all_timers[current_timer]->pause();
         if (timer_parents.find(current_timer) == timer_parents.end())
         {
            current_timer = "";
         } else { 
            current_timer = timer_parents[current_timer];
         }
      }
   }
}


extern "C" void g2g_timer_clear_( void )
{
   if (G2G::timer_sum)
   {
      for (map<string,Timer*>::iterator it = all_timers.begin(); it != all_timers.end(); ++it)
      {
         delete it->second;
      }
      top_timers.clear();
      all_timers.clear();
      timer_parents.clear();
      timer_children.clear();
   }
}

void print_timer(string indent, string timer_name, Timer& timer, float total, string parent) 
{
   float time = timer.getSec() + (float)(timer.getMicrosec()) / 1000000.0f;
   printf("%s%-35s%12.6fs (%6.2f%% of %s)\n",indent.c_str(),timer_name.c_str(),time,(100.0f * time / total),parent.c_str());
   indent.append("  ");
   for (set<string>::iterator it = timer_children[timer_name].begin(); it != timer_children[timer_name].end(); ++it)
   {
      print_timer(indent, *it, *all_timers[(*it)], time, timer_name);
   }
}

extern "C" void g2g_timer_summary_( void )
{
   if (G2G::timer_sum)
   {
      Timer total_timer = *all_timers["Total"];
      float total_time = 0.0f;
      total_time = total_timer.getSec() + (float)(total_timer.getMicrosec()) / 1000000.0f;

      string indent = "";
      cout << "-------------------------------------------------------------------------------------------" << endl;
      cout << "                   LIO TIMING INFORMATION:" << endl;
      cout << "-------------------------------------------------------------------------------------------" << endl << endl;
      cout << "Total time: " << total_time << "s" << endl;
      for (map<string,Timer*>::iterator it = top_timers.begin(); it != top_timers.end(); ++it) {
         print_timer(indent, it->first, *(it->second), total_time, "Total");
      }
      cout << "-------------------------------------------------------------------------------------------" << endl;
   }
}
