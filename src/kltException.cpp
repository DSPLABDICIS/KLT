/*
 * kltException.cpp
 *
 *  Created on: Feb 9, 2013
 *      Author: luzdora
 */

/* $Id: kltException.cpp 1871 2006-07-17 13:17:18Z croussil $ */

#include "kltException.hpp"

#include <sstream>
typedef int ExceptionId; //special structure*********************************************?

using std::string;
using namespace klt;

/*
 *
 * class KltException
 */


KltException::KltException(ExceptionId id_, const string& message_, const string& file_, int line_) throw()
{
//	error handling code
	}

KltException::~KltException() throw() {}

KltException::ExceptionId KltException::getExceptionId() const throw() {
  return id;
}

string KltException::exceptionIdToString(ExceptionId id_) throw() {
  switch(id_) {
//   case MY_ERROR:
//     return "MY_ERROR";

  default:
    std::stringstream s;
    s << id_;
    return s.str();
  }
}

