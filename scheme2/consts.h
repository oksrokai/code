#ifndef CONSTS_H
#define CONSTS_H

#include "params.h"

#define _16Q 0
#define _16QINV 16
#define _16V 32
#define _16ZERO 48
#define _ZETAS_EXP 64
/* The C ABI on MacOS exports all symbols with a leading
 * underscore. This means that any symbols we refer to from
 * C files (functions) can't be found, and all symbols we
 * refer to from ASM also can't be found.
 *
 * This define helps us get around this
 */
#ifdef __ASSEMBLER__
#if defined(__WIN32__) || defined(__APPLE__)
#define cdecl(s) _##s
#else
#define cdecl(s) s
#endif
#endif

#ifndef __ASSEMBLER__
#include "align.h"
typedef ALIGNED_INT16(528) qdata_t;
extern const qdata_t qdata;
#endif

#endif
