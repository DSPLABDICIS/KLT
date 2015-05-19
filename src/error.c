/*
 * error.c
 *
 *  Created on: Feb 9, 2013
 *      Author: luzdora
 */

/*********************************************************************
 * error.c
 *
 * Error and warning messages, and system commands.
 *********************************************************************/


/* Standard includes */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "error.h"

/*********************************************************************
 * KLTError
 *
 * Prints an error message and dies.
 *
 * INPUTS
 * exactly like printf
 */

void KLTError(char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
	char buffer[150];
	printf(buffer, fmt, args);
  va_end(args);
//  JFR_ERROR(KltException, KltException::KLT_ERROR, buffer);
 }


/*********************************************************************
 * KLTWarning
 *
 * Prints a warning message.
 *
 * INPUTS
 * exactly like printf
 */

void KLTWarning(const char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
	char buffer[150];
	printf("Warning : ");
	printf(buffer, fmt, args);
  va_end(args);
}

/*void assert(int predicate)
{
	JFR_ASSERT(predicate, "Standard KLT lib assert");
}*/
