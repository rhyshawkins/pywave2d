
#include <stdio.h>
#include <stdarg.h>

#include "pywave2dexception.hpp"

extern "C" {
  #include "slog.h"
};

pywave2dexception::pywave2dexception(const char *srcfile,
				     const char *function,
				     int lineno,
				     const char *fmt, ...)
{
  va_list ap;
  
  va_start(ap, fmt);
  vslog(SLOG_ERROR,
	srcfile,
	function,
	lineno,
	fmt,
	ap);
  va_end(ap);
  
}

pywave2dexception::~pywave2dexception() throw()
{
}
