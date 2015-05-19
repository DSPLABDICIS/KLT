/*
 * base.h
 *
 *  Created on: Feb 9, 2013
 *      Author: luzdora
 */

#ifndef BASE_H_
#define BASE_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef uchar
#define uchar unsigned char
#endif
//
//#ifndef schar
//#define schar signed char
//#endif
//
//#ifndef uint
//#define uint unsigned int
//#endif
//
//#ifndef ushort
//#define ushort unsigned short
//#endif
//
//#ifndef ulong
//#define ulong unsigned long
//#endif


int max(int a,int b);
int min(int a,int b);

#define max3(a,b,c)	((a) > (b) ? max((a),(c)) : max((b),(c)))
#define min3(a,b,c)	((a) < (b) ? min((a),(c)) : min((b),(c)))


#ifdef __cplusplus
}
#endif

#endif /* BASE_H_ */
