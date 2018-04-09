//
//    pyWave2D : A python 2 extension for inversion of 2D datasets using trans-dimensional
//    trees. For reference see
//
//      R Hawkins and M Sambridge, "Geophysical imaging using trans-dimensional trees",
//      Geophysical Journal International, 2015, 203:2, 972 - 1000,
//      https://doi.org/10.1093/gji/ggv326
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
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
