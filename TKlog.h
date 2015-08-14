
#ifdef __cplusplus
#include <cstdio>
#include <ctime>
extern "C" {
#else
#include <stdio.h>
#include <time.h>
#endif
   extern int debugFlag;   /* Global = 1 if debug mode is running */
   extern int writeResiduals;   /* Global = 1 if we are writing out post-fit residuals */
   extern int tcheck;   /* Global = 1 if time check message should be printed is running */
   extern clock_t timer_clk;

#ifdef __cplusplus
}
#endif

/* define some functions for log message 
 * M.Keith 2012 - let me know if this fails to compile anywhere.
 * mkeith@pulsarastronomy.net
 **/
#ifndef LOG_OUTFILE
#define LOG_OUTFILE stdout
#endif
#define WHERESTR  "[%s:%d] "
#define WHEREARG  __FILE__, __LINE__
#define ENDL "\n"
#define WHEREERR "******\nERROR [%s:%d] "
#define WHERETCHK "[%s:%d] T=%.2f s: "
#define _LOG(...) fprintf(LOG_OUTFILE,__VA_ARGS__)
#define logmsg(_fmt, ...) _LOG(WHERESTR _fmt ENDL, WHEREARG,##__VA_ARGS__)
#define logdbg(_fmt, ...)  if(debugFlag)logmsg(_fmt,##__VA_ARGS__)
#define logerr(_fmt, ...) _LOG(WHEREERR _fmt ENDL, WHEREARG,##__VA_ARGS__)
#define logtchk(_fmt, ...) if(tcheck)_LOG(WHERETCHK _fmt ENDL, WHEREARG,(clock()-timer_clk)/(float)CLOCKS_PER_SEC,##__VA_ARGS__)

