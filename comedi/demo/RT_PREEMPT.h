#ifndef RT_PREEMPT
#define RT_PREEMPT

#include <time.h>
#include <sched.h>
#include <sys/mman.h>

class RT_preempt{
 public:
  int interval;///< real time resolution
  struct timespec t; ///< real time variable
  /*
  //! constructor
  void initialization();
  //! functions
  void nanowait();
  void next_time_interval();
  };*/

void initialization()
{
  interval= DEMIPERIOD; /* 50us*/
  //int interval = 50000; /* 50us*/
  set_RTpriority();  
  lock_memory_pagination();  
  stack_prefault();  
  gettime();  
}

void set_RTpriority() 
{
#define MY_PRIORITY (49) /* we use 49 as the PRREMPT_RT use 50
                            as the priority of kernel tasklets
                            and interrupt handler by default */
  
  struct sched_param param;
  /* Declare ourself as a real time task */
  
  param.sched_priority = MY_PRIORITY;
  if(sched_setscheduler(0, SCHED_FIFO, &param) == -1) 
    {
      perror("sched_setscheduler failed");
      exit(-1);
    }
}
  
void lock_memory_pagination() 
{
  /* Lock memory */
  
  if(mlockall(MCL_CURRENT|MCL_FUTURE) == -1) {
    perror("mlockall failed");
    exit(-2);
  }
  
}
  
void stack_prefault() {
#define MAX_SAFE_STACK (8*1024) /* The maximum stack size which is
                                   guranteed safe to access without
                                   faulting */
  
  unsigned char dummy[MAX_SAFE_STACK];
  
  memset(dummy, 0, MAX_SAFE_STACK);
  return;
}

void gettime()
{
  clock_gettime(CLOCK_MONOTONIC ,&t);
  /* start after one second */
  t.tv_sec++;
  
} 
  
void nanowait()
{
  /* Like nanosleep(2), clock_nanosleep() allows the calling thread to sleep for an
     interval specified with nanosecond precision.*/
  
  clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &t, NULL);
  
}

void next_time_interval()
{
#define NSEC_PER_SEC    (1000000000) /* The number of nsecs per sec. */
  t.tv_nsec += interval;
  
  while (t.tv_nsec >= NSEC_PER_SEC) {
    t.tv_nsec -= NSEC_PER_SEC;
    t.tv_sec++;}
  
}

}; //CLASS RT_PREEMPT

#endif// RT_PREEMPT
