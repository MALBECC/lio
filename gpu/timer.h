#ifndef __TIMER_H__
#define __TIMER_H__

 #include <iostream>

namespace G2G {
	class Timer {
		public:
			Timer(void);
      Timer(const timeval& t);
		
			void start(void);
			void stop(void);
      void pause(void);
		
			unsigned long getMicrosec(void) const;
			unsigned long getSec(void) const;
		
			bool isStarted(void) const;
			
			friend std::ostream& operator<<(std::ostream& o, const Timer& t);
		
			/* to compare stopped timers */
			bool operator<(const Timer& other) const;
		
			void print(void);
		
			static void sync(void);
		
		private:
			timeval t0, t1, res;
			bool started;
	};
	
	std::ostream& operator<<(std::ostream& o, const Timer& t);	
}



#endif
