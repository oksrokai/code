#include "fq.inc"

#if defined(__WIN32__) || defined(__APPLE__)
#define cdecl(s) _##s
#else
#define cdecl(s) s
#endif

.macro schoolbook off,sign
#load
vmovdqa		\off+32(%rsi),%ymm7	# b
vmovdqa		\off+32(%rdx),%ymm8	# d
vmovdqa		\off(%rsi),%ymm9	# a
vmovdqa		\off(%rdx),%ymm10	# c

#mul
vpmullw		%ymm7,%ymm8,%ymm11	# bd.lo
vpmulhw		%ymm7,%ymm8,%ymm12	# bd.hi
vpmullw		%ymm7,%ymm10,%ymm13	# bc.lo
vpmulhw		%ymm7,%ymm10,%ymm7	# bc.hi
vpmullw		%ymm9,%ymm8,%ymm14	# ad.lo
vpmulhw		%ymm9,%ymm8,%ymm8	# ad.hi
vpmullw		%ymm9,%ymm10,%ymm15	# ac.lo
vpmulhw		%ymm9,%ymm10,%ymm9	# ac.hi

#reduce
vpmullw		%ymm1,%ymm11,%ymm11
vpmulhw		%ymm0,%ymm11,%ymm11
vpsubw		%ymm11,%ymm12,%ymm11	# bd

#mul
vpmullw		%ymm2,%ymm11,%ymm10	# rbd.lo
vpmulhw		%ymm2,%ymm11,%ymm11	# rbd.hi

#unpack
vpunpcklwd	%ymm7,%ymm13,%ymm12	# bc0
vpunpckhwd	%ymm7,%ymm13,%ymm7	# bc1
vpunpcklwd	%ymm8,%ymm14,%ymm13	# ad0
vpunpckhwd	%ymm8,%ymm14,%ymm8	# ad1
vpunpcklwd	%ymm9,%ymm15,%ymm14	# ac0
vpunpckhwd	%ymm9,%ymm15,%ymm9	# ac1
vpunpcklwd	%ymm11,%ymm10,%ymm15	# rbd0
vpunpckhwd	%ymm11,%ymm10,%ymm10	# rbd1

#add
.ifeq \sign
vpaddd		%ymm14,%ymm15,%ymm14	# x0
vpaddd		%ymm9,%ymm10,%ymm9	# x1
.else
vpsubd		%ymm15,%ymm14,%ymm14	# x0
vpsubd		%ymm10,%ymm9,%ymm9	# x1
.endif
vpaddd		%ymm12,%ymm13,%ymm12	# y0
vpaddd		%ymm7,%ymm8,%ymm7	# y1
.endm

.macro karatsuba off,sign
#load
vmovdqa		\off+32(%rsi),%ymm7	# b
vmovdqa		\off+32(%rdx),%ymm8	# d
vmovdqa		\off(%rsi),%ymm9    # a
vmovdqa		\off(%rdx),%ymm10   # c

vpaddw      %ymm7,%ymm9,%ymm11  # a+b
vpaddw      %ymm8,%ymm10,%ymm12 # c+d

#mul
vpmullw		%ymm8,%ymm7,%ymm13	# bd.lo
vpmulhw		%ymm8,%ymm7,%ymm14	# bd.hi
vpmullw		%ymm10,%ymm9,%ymm15	# ac.lo
vpmulhw		%ymm10,%ymm9,%ymm9	# ac.hi
vpmullw		%ymm12,%ymm11,%ymm7 # (a+b)(c+d).lo
vpmulhw		%ymm12,%ymm11,%ymm11 # (a+b)(c+d).hi

#reduce
vpmullw		%ymm1,%ymm13,%ymm8
vpmulhw		%ymm0,%ymm8,%ymm8
vpsubw		%ymm8,%ymm14,%ymm8 # bd'

#mul
vpmullw		%ymm2,%ymm8,%ymm10	# rbd.lo
vpmulhw		%ymm2,%ymm8,%ymm8	# rbd.hi

#unpack
vpunpcklwd	%ymm14,%ymm13,%ymm12 # bd0
vpunpckhwd	%ymm14,%ymm13,%ymm13 # bd1
vpunpcklwd	%ymm9,%ymm15,%ymm14	 # ac0
vpunpckhwd	%ymm9,%ymm15,%ymm9	 # ac1
vpunpcklwd	%ymm11,%ymm7,%ymm15   # (a+b)(c+d)0
vpunpckhwd	%ymm11,%ymm7,%ymm7    # (a+b)(c+d)1
vpunpcklwd	%ymm8,%ymm10,%ymm11  # rbd0
vpunpckhwd	%ymm8,%ymm10,%ymm8 # rbd1

#add
vpaddd      %ymm12,%ymm14,%ymm10 # (ac+bd)0=ac0+bd0
vpaddd      %ymm13,%ymm9,%ymm12  # (ac+bd)1=ac1+bd1
.ifeq \sign
vpaddd		%ymm11,%ymm14,%ymm13 # x0=ac0+rbd0
vpaddd		%ymm8,%ymm9,%ymm8  # x1=ac1+rbd1
.else
vpsubd		%ymm11,%ymm14,%ymm13 # x0=ac0-rbd0
vpsubd		%ymm8,%ymm9,%ymm8  # x1=ac1-rbd1
.endif
vpsubd      %ymm10,%ymm15,%ymm10 # y0=(a+b)(c+d)0-(ac+bd)0
vpsubd      %ymm12,%ymm7,%ymm7  # y1=(a+b)(c+d)1-(ac+bd)1
.endm

.macro red a0,a1,b0,b1,x,y,z
vpxor		%ymm\x,%ymm\x,%ymm\x
#a
#pack
vpblendw	$0xAA,%ymm\x,%ymm\a0,%ymm\y
vpblendw	$0xAA,%ymm\x,%ymm\a1,%ymm\z
vpsrld		$16,%ymm\a0,%ymm\a0
vpsrld		$16,%ymm\a1,%ymm\a1
vpackusdw	%ymm\z,%ymm\y,%ymm\z
vpackusdw	%ymm\a1,%ymm\a0,%ymm\a0
#reduce
vpmullw		%ymm1,%ymm\z,%ymm\z
vpmulhw		%ymm0,%ymm\z,%ymm\z
vpsubw		%ymm\z,%ymm\a0,%ymm\a0
#b
#pack
vpblendw	$0xAA,%ymm\x,%ymm\b0,%ymm\y
vpblendw	$0xAA,%ymm\x,%ymm\b1,%ymm\x
vpsrld		$16,%ymm\b0,%ymm\b0
vpsrld		$16,%ymm\b1,%ymm\b1
vpackusdw	%ymm\x,%ymm\y,%ymm\y
vpackusdw	%ymm\b1,%ymm\b0,%ymm\b0
#reduce
vpmullw		%ymm1,%ymm\y,%ymm\y
vpmulhw		%ymm0,%ymm\y,%ymm\y
vpsubw		%ymm\y,%ymm\b0,%ymm\b0
.endm

.global cdecl(basemul_acc_avx)
cdecl(basemul_acc_avx):
#consts
vmovdqa		cdecl(_16xq)(%rip),%ymm0
vmovdqa		cdecl(_16xqinv)(%rip),%ymm1
vmovdqu		(%rcx),%ymm2

poly0.0:
karatsuba	0,0

#mov
vmovdqa		%ymm13,%ymm3 # x0
vmovdqa		%ymm8,%ymm4 # x1
vmovdqa		%ymm10,%ymm5 # y0
vmovdqa		%ymm7,%ymm6	 # y1

poly1.0:
karatsuba	512,0

#add
vpaddd		%ymm13,%ymm3,%ymm3
vpaddd		%ymm8,%ymm4,%ymm4
vpaddd		%ymm10,%ymm5,%ymm5
vpaddd		%ymm7,%ymm6,%ymm6

#reduce
red		3,4,5,6,7,8,9

#if OKAI_K > 2
poly2.0:
karatsuba	1024,0

#reduce
red		13,8,10,7,9,11,12

vpaddw		%ymm13,%ymm3,%ymm3
vpaddw		%ymm10,%ymm5,%ymm5
#endif

#store
vmovdqa		%ymm3,(%rdi)
vmovdqa		%ymm5,32(%rdi)

poly0.1:
karatsuba	64,1

#mov
vmovdqa		%ymm13,%ymm3
vmovdqa		%ymm8,%ymm4
vmovdqa		%ymm10,%ymm5
vmovdqa		%ymm7,%ymm6

poly1.1:
karatsuba	576,1

#add
vpaddd		%ymm13,%ymm3,%ymm3
vpaddd		%ymm8,%ymm4,%ymm4
vpaddd		%ymm10,%ymm5,%ymm5
vpaddd		%ymm7,%ymm6,%ymm6

#reduce
red		3,4,5,6,7,8,9

#if OKAI_K > 2
poly2.1:
karatsuba	1088,1

#reduce
red		13,8,10,7,9,11,12

vpaddw		%ymm13,%ymm3,%ymm3
vpaddw		%ymm10,%ymm5,%ymm5
#endif

#store
vmovdqa		%ymm3,64(%rdi)
vmovdqa		%ymm5,96(%rdi)

ret

.global cdecl(basemul_avx)
cdecl(basemul_avx):
#consts
vmovdqa		cdecl(_16xq)(%rip),%ymm0
vmovdqa		cdecl(_16xqinv)(%rip),%ymm1
vmovdqu		(%rcx),%ymm2

schoolbook	0,0

#reduce
red		14,9,12,7,8,10,11

#store
vmovdqa		%ymm14,(%rdi)
vmovdqa		%ymm12,32(%rdi)

schoolbook	64,1

#reduce
red		14,9,12,7,8,10,11

#store
vmovdqa		%ymm14,64(%rdi)
vmovdqa		%ymm12,96(%rdi)

ret

