#pragma once
#ifndef pywave2dexception_hpp
#define pywave2dexception_hpp

#include <exception>

#define PYWAVE2DEXCEPTION(fmt, ...) pywave2dexception(__FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)

class pywave2dexception : public std::exception {
public:

  
  pywave2dexception(const char *srcfile,
		      const char *function,
		      int lineno,
		      const char *fmt, ...);
  virtual ~pywave2dexception() throw();
  
};

#endif // pywave2dexception_hpp
