/*
 * error.h
 *
 *  Created on: Feb 9, 2013
 *      Author: luzdora
 */

#ifndef ERROR_H_
#define ERROR_H_

#include <stdio.h>
#include <stdarg.h>

#define PRECOND(predicate, string)\
	if(!(predicate)){\
		printf(string);\
		printf(#predicate);\
		throw;\
}

#ifdef __cplusplus
extern "C" {
#endif

void KLTError(char *fmt, ...);
void KLTWarning(const char *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif /* ERROR_H_ */
