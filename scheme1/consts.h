#ifndef CONSTS_H
#define CONSTS_H

#include "params.h"

#define _16XQ            0
#define _16XQINV        16
#define _16XV           32
#define _16XFLO         48
#define _16XFHI         64
#define _16XMONTSQLO    80
#define _16XMONTSQHI    96
#define _16XMASK       112
#define _REVIDXB       128
#define _REVIDXD       144
#define _ZETAS_EXP     160
#define	_16XSHIFT      624
#define _16XHALFQF     640
#define _16XHALFQC     656
#define _16XONE        672
#define _16XTWO        688
#define _16XSEVEN      704
#define _8XONE         720
#define _8XFOUR        736
#define _8XEIGHT       752
#define _8XFIFTEEN     768
#define _8XTHIRTYONE   784
#define _8X512         800
#define _8X1024        816
#define _8XQ           832
#define _8XHALFQF      848
#define _8XHALFQC      864
#define _16XMAGIC      880
#define _8XMAGIC       896
#define _8XMAGICG      912
#define _32XMASKSHUF   928
#define _16XMASK10     944
#define _16XMASK11     960
#define _16XZERO       976
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
typedef ALIGNED_INT16(992) qdata_t;
extern const qdata_t qdata;
#endif

#endif
