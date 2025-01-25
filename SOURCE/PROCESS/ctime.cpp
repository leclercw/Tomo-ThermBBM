#include <time.h> 
#include <sys/time.h> 

#include "ctime.h"

R CPUtime(){
#ifdef SYSTIMES
  struct tms buf;
  if (times(&buf)!=-1)
    return ((R)buf.tms_utime+(R)buf.tms_stime)/(long) sysconf(_SC_CLK_TCK);
  else
#endif
    return ((R) clock())/CLOCKS_PER_SEC;
}
