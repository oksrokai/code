#include "consts.h"

.macro dec_from_10  a,b=11,c=12
vpmullw		%ymm14,%ymm\a,%ymm\b
vpmulhw		%ymm14,%ymm\a,%ymm\a
vpunpcklwd	%ymm\a,%ymm\b,%ymm\c
vpunpckhwd	%ymm\a,%ymm\b,%ymm\b
vpaddd      %ymm13,%ymm\c,%ymm\c
vpaddd      %ymm13,%ymm\b,%ymm\b
vpsrld      $10,%ymm\c,%ymm\c
vpsrld      $10,%ymm\b,%ymm\b
vpackusdw   %ymm\b,%ymm\c,%ymm\a
.endm

.macro dec_from_11  a,b=11,c=12
vpmullw		%ymm14,%ymm\a,%ymm\b
vpmulhw		%ymm14,%ymm\a,%ymm\a
vpunpcklwd	%ymm\a,%ymm\b,%ymm\c
vpunpckhwd	%ymm\a,%ymm\b,%ymm\b
vpaddd      %ymm13,%ymm\c,%ymm\c
vpaddd      %ymm13,%ymm\b,%ymm\b
vpsrld      $11,%ymm\c,%ymm\c
vpsrld      $11,%ymm\b,%ymm\b
vpackusdw   %ymm\b,%ymm\c,%ymm\a
.endm

.text
.global cdecl(poly_decompress_from_10_avx2)
cdecl(poly_decompress_from_10_avx2):
#const
vmovdqa		_16XMASK10*2(%rdx),%ymm15
vmovdqa		_16XQ*2(%rdx),%ymm14
vmovdqa		_8X512*2(%rdx),%ymm13

#load
vmovdqu		(%rsi),%ymm0
vmovdqu		32(%rsi),%ymm1
vmovdqu		64(%rsi),%ymm2
vmovdqu		96(%rsi),%ymm3
vmovdqu		128(%rsi),%ymm4
vmovdqu		160(%rsi),%ymm5
vmovdqu		192(%rsi),%ymm6
vmovdqu		224(%rsi),%ymm7
vmovdqu		256(%rsi),%ymm8
vmovdqu		288(%rsi),%ymm9

vpand		%ymm15,%ymm0,%ymm10
dec_from_10 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi

vpsrlw		$10,%ymm0,%ymm0
vpsllw		$6,%ymm1,%ymm10
vpaddw		%ymm10,%ymm0,%ymm0
vpand		%ymm15,%ymm0,%ymm0
dec_from_10 0
vmovdqa		%ymm0,(%rdi)
add		    $32,%rdi

vpsrlw		$4,%ymm1,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_10 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi

vpsrlw		$14,%ymm1,%ymm1
vpsllw		$2,%ymm2,%ymm10
vpaddw		%ymm10,%ymm1,%ymm1
vpand		%ymm15,%ymm1,%ymm1
dec_from_10 1
vmovdqa		%ymm1,(%rdi)
add		    $32,%rdi

vpsrlw		$8,%ymm2,%ymm2
vpsllw		$8,%ymm3,%ymm10
vpaddw		%ymm10,%ymm2,%ymm2
vpand		%ymm15,%ymm2,%ymm2
dec_from_10 2
vmovdqa		%ymm2,(%rdi)
add		    $32,%rdi

vpsrlw		$2,%ymm3,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_10 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi

vpsrlw		$12,%ymm3,%ymm3
vpsllw		$4,%ymm4,%ymm10
vpaddw		%ymm10,%ymm3,%ymm3
vpand		%ymm15,%ymm3,%ymm3
dec_from_10 3
vmovdqa		%ymm3,(%rdi)
add		    $32,%rdi

vpsrlw		$6,%ymm4,%ymm4
dec_from_10 4
vmovdqa		%ymm4,(%rdi)
add		    $32,%rdi

vpand		%ymm15,%ymm5,%ymm10
dec_from_10 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi

vpsrlw		$10,%ymm5,%ymm5
vpsllw		$6,%ymm6,%ymm10
vpaddw		%ymm10,%ymm5,%ymm5
vpand		%ymm15,%ymm5,%ymm5
dec_from_10 5
vmovdqa		%ymm5,(%rdi)
add		    $32,%rdi

vpsrlw		$4,%ymm6,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_10 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi

vpsrlw		$14,%ymm6,%ymm6
vpsllw		$2,%ymm7,%ymm10
vpaddw		%ymm10,%ymm6,%ymm6
vpand		%ymm15,%ymm6,%ymm6
dec_from_10 6
vmovdqa		%ymm6,(%rdi)
add		    $32,%rdi

vpsrlw		$8,%ymm7,%ymm7
vpsllw		$8,%ymm8,%ymm10
vpaddw		%ymm10,%ymm7,%ymm7
vpand		%ymm15,%ymm7,%ymm7
dec_from_10 7
vmovdqa		%ymm7,(%rdi)
add		    $32,%rdi

vpsrlw		$2,%ymm8,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_10 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi

vpsrlw		$12,%ymm8,%ymm8
vpsllw		$4,%ymm9,%ymm10
vpaddw		%ymm10,%ymm8,%ymm8
vpand		%ymm15,%ymm8,%ymm8
dec_from_10 8
vmovdqa		%ymm8,(%rdi)
add		    $32,%rdi

vpsrlw		$6,%ymm9,%ymm9
dec_from_10 9
vmovdqa		%ymm9,(%rdi)

ret

.global cdecl(poly_decompress_from_11_avx2)
cdecl(poly_decompress_from_11_avx2):
#const
vmovdqa		_16XMASK11*2(%rdx),%ymm15
vmovdqa		_16XQ*2(%rdx),%ymm14
vmovdqa		_8X1024*2(%rdx),%ymm13

#load
vmovdqu		(%rsi),%ymm0
vmovdqu		32(%rsi),%ymm1
vmovdqu		64(%rsi),%ymm2
vmovdqu		96(%rsi),%ymm3
vmovdqu		128(%rsi),%ymm4
vmovdqu		160(%rsi),%ymm5
vmovdqu		192(%rsi),%ymm6
vmovdqu		224(%rsi),%ymm7
vmovdqu		256(%rsi),%ymm8
vmovdqu		288(%rsi),%ymm9

vpand		%ymm15,%ymm0,%ymm10
dec_from_11 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi  #0

vpsrlw		$11,%ymm0,%ymm0
vpsllw		$5,%ymm1,%ymm10
vpaddw		%ymm10,%ymm0,%ymm0
vpand		%ymm15,%ymm0,%ymm0
dec_from_11 0
vmovdqa		%ymm0,(%rdi)
add		    $32,%rdi   #1

vmovdqu		320(%rsi),%ymm0

vpsrlw		$6,%ymm1,%ymm1
vpsllw		$10,%ymm2,%ymm10
vpaddw		%ymm10,%ymm1,%ymm1
vpand		%ymm15,%ymm1,%ymm1
dec_from_11 1
vmovdqa		%ymm1,(%rdi)
add		    $32,%rdi   #2

vpsrlw		$1,%ymm2,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_11 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi   #3

vpsrlw		$12,%ymm2,%ymm2
vpsllw		$4,%ymm3,%ymm10
vpaddw		%ymm10,%ymm2,%ymm2
vpand		%ymm15,%ymm2,%ymm2
dec_from_11 2
vmovdqa		%ymm2,(%rdi)
add		    $32,%rdi     #4

vpsrlw		$7,%ymm3,%ymm3
vpsllw		$9,%ymm4,%ymm10
vpaddw		%ymm10,%ymm3,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_11 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi   #5

vpsrlw		$2,%ymm4,%ymm4
vpand		%ymm15,%ymm4,%ymm4
dec_from_11 4
vmovdqa		%ymm4,(%rdi)
add		    $32,%rdi   #6

vpsrlw		$13,%ymm4,%ymm4
vpsllw		$3,%ymm5,%ymm10
vpaddw		%ymm10,%ymm4,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_11 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi  #7

vpsrlw		$8,%ymm5,%ymm5
vpsllw		$8,%ymm6,%ymm10
vpaddw		%ymm10,%ymm5,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_11 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi  #8

vpsrlw		$3,%ymm6,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_11 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi   #9

vpsrlw		$14,%ymm6,%ymm6
vpsllw		$2,%ymm7,%ymm10
vpaddw		%ymm10,%ymm6,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_11 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi   #10

vpsrlw		$9,%ymm7,%ymm7
vpsllw		$7,%ymm8,%ymm10
vpaddw		%ymm10,%ymm7,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_11 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi   #11

vpsrlw		$4,%ymm8,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_11 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi   #12

vpsrlw		$15,%ymm8,%ymm8
vpsllw		$1,%ymm9,%ymm10
vpaddw		%ymm10,%ymm8,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_11 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi   #13

vpsrlw		$10,%ymm9,%ymm9
vpsllw		$6,%ymm0,%ymm10
vpaddw		%ymm10,%ymm9,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_11 10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi  #14

vpsrlw		$5,%ymm0,%ymm0
dec_from_11 0
vmovdqa		%ymm0,(%rdi)

ret
