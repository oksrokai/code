#include "consts.h"

.macro schoolbook_multi_pos off
vmovdqa		_16XQINV*2(%rcx),%ymm0
vmovdqa		(32*\off+ 0)*2(%rdi),%ymm1		# f0[0,8,...,120]
vmovdqa		(32*\off+16)*2(%rdi),%ymm2		# f0[1,9,...,121]
vmovdqa		(32*\off+ 0)*2(%rsi),%ymm3		# g0[0,8,...,120]
vmovdqa		(32*\off+16)*2(%rsi),%ymm4		# g0[1,9,...,121]

vmovdqa		(32*\off+256)*2(%rdi),%ymm5		# f1[0,8,...,120]
vmovdqa		(32*\off+272)*2(%rdi),%ymm6		# f1[1,9,...,121]
vmovdqa		(32*\off+256)*2(%rsi),%ymm7		# g1[0,8,...,120]
vmovdqa		(32*\off+272)*2(%rsi),%ymm8		# g1[1,9,...,121]

vpmullw		%ymm0,%ymm1,%ymm9			# f0*QINV 0,...,120
vpmullw		%ymm0,%ymm2,%ymm10			# f0*QINV 1,...,121
vpmullw		%ymm0,%ymm5,%ymm11			# f1*QINV 0,...,120
vpmullw		%ymm0,%ymm6,%ymm12			# f1*QINV 1,...,121

vpmulhw		%ymm1,%ymm3,%ymm13			# f0[0]*g0[0]
vpmulhw		%ymm1,%ymm4,%ymm1			# f0[0]*g0[1]
vpmulhw		%ymm2,%ymm3,%ymm14			# f0[1]*g0[0]
vpmulhw		%ymm2,%ymm4,%ymm2			# f0[1]*g0[1]

vpmulhw		%ymm5,%ymm7,%ymm15			# f1[0]*g1[0]
vpmulhw		%ymm5,%ymm8,%ymm5			# f1[0]*g1[1]
vpmulhw		%ymm6,%ymm7,%ymm0			# f1[1]*g1[0]
vpmulhw		%ymm6,%ymm8,%ymm6			# f1[1]*g1[1]

vmovdqa		%ymm13,(%rsp)

vpmullw		%ymm3,%ymm9,%ymm13			# f0[0]*g0[0]*QINV
vpmullw		%ymm4,%ymm9,%ymm9			# f0[0]*g0[1]*QINV
vpmullw		%ymm3,%ymm10,%ymm3			# f0[1]*g0[0]*QINV
vpmullw		%ymm4,%ymm10,%ymm10         # f0[1]*g0[1]*QINV

vpmullw		%ymm7,%ymm11,%ymm4			# f1[0]*g1[0]*QINV
vpmullw		%ymm8,%ymm11,%ymm11			# f1[0]*g1[1]*QINV
vpmullw		%ymm7,%ymm12,%ymm7			# f1[1]*g1[0]*QINV
vpmullw		%ymm8,%ymm12,%ymm12			# f1[1]*g1[1]*QINV

vmovdqa		_16XQ*2(%rcx),%ymm8
vpmulhw		%ymm8,%ymm13,%ymm13
vpmulhw		%ymm8,%ymm9,%ymm9
vpmulhw		%ymm8,%ymm3,%ymm3
vpmulhw		%ymm8,%ymm10,%ymm10
vpmulhw		%ymm8,%ymm4,%ymm4
vpmulhw		%ymm8,%ymm11,%ymm11
vpmulhw		%ymm8,%ymm7,%ymm7
vpmulhw		%ymm8,%ymm12,%ymm12

vpsubw		(%rsp),%ymm13,%ymm13		# -f0[0]*g0[0]
vpsubw		%ymm9,%ymm1,%ymm9			# f0[0]*g0[1]
vpsubw		%ymm3,%ymm14,%ymm3			# f0[1]*g0[0]
vpsubw		%ymm10,%ymm2,%ymm10			# f0[1]*g0[1]

vpsubw		%ymm4,%ymm15,%ymm4			# f1[0]*g1[0]
vpsubw		%ymm11,%ymm5,%ymm11			# f1[0]*g1[1]
vpsubw		%ymm7,%ymm0,%ymm7			# f1[1]*g1[0]
vpsubw		%ymm12,%ymm6,%ymm12			# f1[1]*g1[1]

vmovdqa		(%r9),%ymm0
vmovdqa		32(%r9),%ymm1
vpmullw		%ymm0,%ymm10,%ymm2
vpmullw		%ymm0,%ymm12,%ymm6
vpmulhw		%ymm1,%ymm10,%ymm10
vpmulhw		%ymm1,%ymm12,%ymm12
vpmulhw		%ymm8,%ymm2,%ymm2
vpmulhw		%ymm8,%ymm6,%ymm6
vpsubw		%ymm2,%ymm10,%ymm10			# rb0d0
vpsubw		%ymm6,%ymm12,%ymm12

vpsubw		%ymm13,%ymm10,%ymm13        # p0hat[0]
vpaddw		%ymm3,%ymm9,%ymm9           # p0hat[1]
vpaddw		%ymm12,%ymm4,%ymm15         # p1hat[0]
vpaddw		%ymm7,%ymm11,%ymm11         # p1hat[1]

vmovdqa		(32*\off+ 0)*2(%rdi),%ymm1		# f0[0,8,...,120]
vmovdqa		(32*\off+16)*2(%rdi),%ymm2		# f0[1,9,...,121]
vmovdqa		(32*\off+ 0)*2(%rsi),%ymm3		# g0[0,8,...,120]
vmovdqa		(32*\off+16)*2(%rsi),%ymm4		# g0[1,9,...,121]

vmovdqa		(32*\off+256)*2(%rdi),%ymm5		# f1[0,8,...,120]
vmovdqa		(32*\off+272)*2(%rdi),%ymm6		# f1[1,9,...,121]
vmovdqa		(32*\off+256)*2(%rsi),%ymm7		# g1[0,8,...,120]
vmovdqa		(32*\off+272)*2(%rsi),%ymm8		# g1[1,9,...,121]

vpaddw      %ymm1,%ymm5,%ymm1               # F[0]
vpaddw      %ymm2,%ymm6,%ymm2               # F[1]
vpaddw      %ymm3,%ymm7,%ymm3               # G[0]
vpaddw      %ymm4,%ymm8,%ymm4               # G[1]

vmovdqa		_16XQINV*2(%rcx),%ymm0
vpmullw		%ymm0,%ymm1,%ymm5			# F[0]*QINV
vpmullw		%ymm0,%ymm2,%ymm6			# F[1]*QINV

vpmulhw		%ymm1,%ymm3,%ymm7			# F[0]*G[0]
vpmulhw		%ymm1,%ymm4,%ymm8			# F[0]*G[1]
vpmulhw		%ymm2,%ymm3,%ymm10			# F[1]*G[0]
vpmulhw		%ymm2,%ymm4,%ymm12			# F[1]*G[1]

vpmullw		%ymm5,%ymm3,%ymm1			# F[0]*G[0]*QINV
vpmullw		%ymm5,%ymm4,%ymm2			# F[0]*G[1]*QINV
vpmullw		%ymm6,%ymm3,%ymm3			# F[1]*G[0]*QINV
vpmullw		%ymm6,%ymm4,%ymm4           # F[1]*G[1]*QINV

vmovdqa		_16XQ*2(%rcx),%ymm14
vpmulhw		%ymm14,%ymm1,%ymm1
vpmulhw		%ymm14,%ymm2,%ymm2
vpmulhw		%ymm14,%ymm3,%ymm3
vpmulhw		%ymm14,%ymm4,%ymm4

vpsubw		%ymm1,%ymm7,%ymm1		    # F[0]*G[0]
vpsubw		%ymm2,%ymm8,%ymm2			# F[0]*G[1]
vpsubw		%ymm3,%ymm10,%ymm3			# F[1]*G[0]
vpsubw		%ymm4,%ymm12,%ymm4			# F[1]*G[1]

vmovdqa		(%r9),%ymm0
vmovdqa		32(%r9),%ymm7
vpmullw		%ymm0,%ymm4,%ymm5
vpmulhw		%ymm7,%ymm4,%ymm6
vpmulhw		%ymm14,%ymm5,%ymm5
vpsubw		%ymm5,%ymm6,%ymm6

vpaddw		%ymm1,%ymm6,%ymm1        # phat[0]
vpaddw		%ymm2,%ymm3,%ymm2        # phat[1]

vpmullw		%ymm0,%ymm11,%ymm3
vpmulhw		%ymm7,%ymm11,%ymm4
vpmulhw		%ymm14,%ymm3,%ymm3
vpsubw		%ymm3,%ymm4,%ymm3

vpaddw      %ymm13,%ymm3,%ymm5    # h0hat[0]
vpaddw      %ymm9,%ymm15,%ymm6    # h0hat[1]

vpsubw      %ymm13,%ymm1,%ymm1
vpsubw      %ymm15,%ymm1,%ymm1    # h1hat[0]
vpsubw      %ymm9,%ymm2,%ymm2
vpsubw      %ymm11,%ymm2,%ymm2    # h1hat[1]

vmovdqa		%ymm5,(32*\off+ 0)*2(%rdx)
vmovdqa		%ymm6,(32*\off+16)*2(%rdx)
vmovdqa		%ymm1,(32*\off+256)*2(%rdx)
vmovdqa		%ymm2,(32*\off+272)*2(%rdx)
.endm

.macro schoolbook_multi_neg off
vmovdqa		_16XQINV*2(%rcx),%ymm0
vmovdqa		(32*\off+ 0)*2(%rdi),%ymm1		# f0[0,8,...,120]
vmovdqa		(32*\off+16)*2(%rdi),%ymm2		# f0[1,9,...,121]
vmovdqa		(32*\off+ 0)*2(%rsi),%ymm3		# g0[0,8,...,120]
vmovdqa		(32*\off+16)*2(%rsi),%ymm4		# g0[1,9,...,121]

vmovdqa		(32*\off+256)*2(%rdi),%ymm5		# f1[0,8,...,120]
vmovdqa		(32*\off+272)*2(%rdi),%ymm6		# f1[1,9,...,121]
vmovdqa		(32*\off+256)*2(%rsi),%ymm7		# g1[0,8,...,120]
vmovdqa		(32*\off+272)*2(%rsi),%ymm8		# g1[1,9,...,121]

vpmullw		%ymm0,%ymm1,%ymm9			# f0*QINV 0,...,120
vpmullw		%ymm0,%ymm2,%ymm10			# f0*QINV 1,...,121
vpmullw		%ymm0,%ymm5,%ymm11			# f1*QINV 0,...,120
vpmullw		%ymm0,%ymm6,%ymm12			# f1*QINV 1,...,121

vpmulhw		%ymm1,%ymm3,%ymm13			# f0[0]*g0[0]
vpmulhw		%ymm1,%ymm4,%ymm1			# f0[0]*g0[1]
vpmulhw		%ymm2,%ymm3,%ymm14			# f0[1]*g0[0]
vpmulhw		%ymm2,%ymm4,%ymm2			# f0[1]*g0[1]

vpmulhw		%ymm5,%ymm7,%ymm15			# f1[0]*g1[0]
vpmulhw		%ymm5,%ymm8,%ymm5			# f1[0]*g1[1]
vpmulhw		%ymm6,%ymm7,%ymm0			# f1[1]*g1[0]
vpmulhw		%ymm6,%ymm8,%ymm6			# f1[1]*g1[1]

vmovdqa		%ymm13,(%rsp)

vpmullw		%ymm3,%ymm9,%ymm13			# f0[0]*g0[0]*QINV
vpmullw		%ymm4,%ymm9,%ymm9			# f0[0]*g0[1]*QINV
vpmullw		%ymm3,%ymm10,%ymm3			# f0[1]*g0[0]*QINV
vpmullw		%ymm4,%ymm10,%ymm10         # f0[1]*g0[1]*QINV

vpmullw		%ymm7,%ymm11,%ymm4			# f1[0]*g1[0]*QINV
vpmullw		%ymm8,%ymm11,%ymm11			# f1[0]*g1[1]*QINV
vpmullw		%ymm7,%ymm12,%ymm7			# f1[1]*g1[0]*QINV
vpmullw		%ymm8,%ymm12,%ymm12			# f1[1]*g1[1]*QINV

vmovdqa		_16XQ*2(%rcx),%ymm8
vpmulhw		%ymm8,%ymm13,%ymm13
vpmulhw		%ymm8,%ymm9,%ymm9
vpmulhw		%ymm8,%ymm3,%ymm3
vpmulhw		%ymm8,%ymm10,%ymm10
vpmulhw		%ymm8,%ymm4,%ymm4
vpmulhw		%ymm8,%ymm11,%ymm11
vpmulhw		%ymm8,%ymm7,%ymm7
vpmulhw		%ymm8,%ymm12,%ymm12

vpsubw		(%rsp),%ymm13,%ymm13		# -f0[0]*g0[0]
vpsubw		%ymm9,%ymm1,%ymm9			# f0[0]*g0[1]
vpsubw		%ymm3,%ymm14,%ymm3			# f0[1]*g0[0]
vpsubw		%ymm10,%ymm2,%ymm10			# f0[1]*g0[1]

vpsubw		%ymm4,%ymm15,%ymm4			# f1[0]*g1[0]
vpsubw		%ymm11,%ymm5,%ymm11			# f1[0]*g1[1]
vpsubw		%ymm7,%ymm0,%ymm7			# f1[1]*g1[0]
vpsubw		%ymm12,%ymm6,%ymm12			# f1[1]*g1[1]

vmovdqa		(%r9),%ymm0
vmovdqa		32(%r9),%ymm1
vpmullw		%ymm0,%ymm10,%ymm2
vpmullw		%ymm0,%ymm12,%ymm6
vpmulhw		%ymm1,%ymm10,%ymm10
vpmulhw		%ymm1,%ymm12,%ymm12
vpmulhw		%ymm8,%ymm2,%ymm2
vpmulhw		%ymm8,%ymm6,%ymm6
vpsubw		%ymm2,%ymm10,%ymm10			# rb0d0
vpsubw		%ymm6,%ymm12,%ymm12

vmovdqa		_16XZERO*2(%rcx),%ymm1
vpsubw      %ymm13,%ymm1,%ymm13
vpsubw		%ymm10,%ymm13,%ymm13        # p0hat[0]
vpaddw		%ymm3,%ymm9,%ymm9           # p0hat[1]
vpsubw		%ymm12,%ymm4,%ymm15         # p1hat[0]
vpaddw		%ymm7,%ymm11,%ymm11         # p1hat[1]

vmovdqa		(32*\off+ 0)*2(%rdi),%ymm1		# f0[0,8,...,120]
vmovdqa		(32*\off+16)*2(%rdi),%ymm2		# f0[1,9,...,121]
vmovdqa		(32*\off+ 0)*2(%rsi),%ymm3		# g0[0,8,...,120]
vmovdqa		(32*\off+16)*2(%rsi),%ymm4		# g0[1,9,...,121]

vmovdqa		(32*\off+256)*2(%rdi),%ymm5		# f1[0,8,...,120]
vmovdqa		(32*\off+272)*2(%rdi),%ymm6		# f1[1,9,...,121]
vmovdqa		(32*\off+256)*2(%rsi),%ymm7		# g1[0,8,...,120]
vmovdqa		(32*\off+272)*2(%rsi),%ymm8		# g1[1,9,...,121]

vpaddw      %ymm1,%ymm5,%ymm1               # F[0]
vpaddw      %ymm2,%ymm6,%ymm2               # F[1]
vpaddw      %ymm3,%ymm7,%ymm3               # G[0]
vpaddw      %ymm4,%ymm8,%ymm4               # G[1]

vmovdqa		_16XQINV*2(%rcx),%ymm0
vpmullw		%ymm0,%ymm1,%ymm5			# F[0]*QINV
vpmullw		%ymm0,%ymm2,%ymm6			# F[1]*QINV

vpmulhw		%ymm1,%ymm3,%ymm7			# F[0]*G[0]
vpmulhw		%ymm1,%ymm4,%ymm8			# F[0]*G[1]
vpmulhw		%ymm2,%ymm3,%ymm10			# F[1]*G[0]
vpmulhw		%ymm2,%ymm4,%ymm12			# F[1]*G[1]

vpmullw		%ymm5,%ymm3,%ymm1			# F[0]*G[0]*QINV
vpmullw		%ymm5,%ymm4,%ymm2			# F[0]*G[1]*QINV
vpmullw		%ymm6,%ymm3,%ymm3			# F[1]*G[0]*QINV
vpmullw		%ymm6,%ymm4,%ymm4           # F[1]*G[1]*QINV

vmovdqa		_16XQ*2(%rcx),%ymm14
vpmulhw		%ymm14,%ymm1,%ymm1
vpmulhw		%ymm14,%ymm2,%ymm2
vpmulhw		%ymm14,%ymm3,%ymm3
vpmulhw		%ymm14,%ymm4,%ymm4

vpsubw		%ymm1,%ymm7,%ymm1		    # F[0]*G[0]
vpsubw		%ymm2,%ymm8,%ymm2			# F[0]*G[1]
vpsubw		%ymm3,%ymm10,%ymm3			# F[1]*G[0]
vpsubw		%ymm4,%ymm12,%ymm4			# F[1]*G[1]

vmovdqa		(%r9),%ymm0
vmovdqa		32(%r9),%ymm7
vpmullw		%ymm0,%ymm4,%ymm5
vpmulhw		%ymm7,%ymm4,%ymm6
vpmulhw		%ymm14,%ymm5,%ymm5
vpsubw		%ymm5,%ymm6,%ymm6

vpsubw		%ymm6,%ymm1,%ymm1        # phat[0]
vpaddw		%ymm2,%ymm3,%ymm2        # phat[1]

vpmullw		%ymm0,%ymm11,%ymm3
vpmulhw		%ymm7,%ymm11,%ymm4
vpmulhw		%ymm14,%ymm3,%ymm3
vpsubw		%ymm3,%ymm4,%ymm3

vpsubw      %ymm3,%ymm13,%ymm5    # h0hat[0]
vpaddw      %ymm9,%ymm15,%ymm6    # h0hat[1]

vpsubw      %ymm13,%ymm1,%ymm1
vpsubw      %ymm15,%ymm1,%ymm1    # h1hat[0]
vpsubw      %ymm9,%ymm2,%ymm2
vpsubw      %ymm11,%ymm2,%ymm2    # h1hat[1]

vmovdqa		%ymm5,(32*\off+ 0)*2(%rdx)
vmovdqa		%ymm6,(32*\off+16)*2(%rdx)
vmovdqa		%ymm1,(32*\off+256)*2(%rdx)
vmovdqa		%ymm2,(32*\off+272)*2(%rdx)
.endm

.text
.global cdecl(basemul_multi_avx)
cdecl(basemul_multi_avx):
mov		%rsp,%r8
and		$-32,%rsp
sub		$32,%rsp

lea		(_ZETAS_EXP+176)*2(%rcx),%r9
schoolbook_multi_pos 0
schoolbook_multi_neg 1

add		$32*2,%r9
schoolbook_multi_pos 2
schoolbook_multi_neg 3

add		$192*2,%r9
schoolbook_multi_pos 4
schoolbook_multi_neg 5

add		$32*2,%r9
schoolbook_multi_pos 6
schoolbook_multi_neg 7

mov		%r8,%rsp
ret