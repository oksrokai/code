#if defined(__WIN32__) || defined(__APPLE__)
#define cdecl(s) _##s
#else
#define cdecl(s) s
#endif

.macro dec_from_8  a,b=11,c=12
vpmullw		%ymm14,%ymm\a,%ymm\b
vpmulhw		%ymm14,%ymm\a,%ymm\a
vpunpcklwd	%ymm\a,%ymm\b,%ymm\c
vpunpckhwd	%ymm\a,%ymm\b,%ymm\b
vpaddd      %ymm13,%ymm\c,%ymm\c
vpaddd      %ymm13,%ymm\b,%ymm\b
vpsrld      $8,%ymm\c,%ymm\c
vpsrld      $8,%ymm\b,%ymm\b
vpackusdw   %ymm\b,%ymm\c,%ymm\a
.endm

.macro dec_from_9  a,b=11,c=12
vpmullw		%ymm14,%ymm\a,%ymm\b
vpmulhw		%ymm14,%ymm\a,%ymm\a
vpunpcklwd	%ymm\a,%ymm\b,%ymm\c
vpunpckhwd	%ymm\a,%ymm\b,%ymm\b
vpaddd      %ymm13,%ymm\c,%ymm\c
vpaddd      %ymm13,%ymm\b,%ymm\b
vpsrld      $9,%ymm\c,%ymm\c
vpsrld      $9,%ymm\b,%ymm\b
vpackusdw   %ymm\b,%ymm\c,%ymm\a
.endm

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

.global cdecl(poly_decompress_from_8_avx2)
cdecl(poly_decompress_from_8_avx2):
.p2align 5
#const
vmovdqa		cdecl(_16xmask8)(%rip),%ymm15
vmovdqa		cdecl(_16xq)(%rip),%ymm14
vmovdqa		cdecl(_8x128)(%rip),%ymm13
#load
vmovdqu		(%rsi),%ymm0
vmovdqu		32(%rsi),%ymm1
vmovdqu		64(%rsi),%ymm2
vmovdqu		96(%rsi),%ymm3
vmovdqu		128(%rsi),%ymm4
vmovdqu		160(%rsi),%ymm5
vmovdqu		192(%rsi),%ymm6
vmovdqu		224(%rsi),%ymm7

vpand		%ymm15,%ymm0,%ymm10
dec_from_8  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi

vpsrlw		$8,%ymm0,%ymm0
vpand		%ymm15,%ymm0,%ymm0
dec_from_8  0
vmovdqa		%ymm0,(%rdi)
add		    $32,%rdi

vpand		%ymm15,%ymm1,%ymm10
dec_from_8  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi

vpsrlw		$8,%ymm1,%ymm1
vpand		%ymm15,%ymm1,%ymm1
dec_from_8  1
vmovdqa		%ymm1,(%rdi)
add		    $32,%rdi

vpand		%ymm15,%ymm2,%ymm10
dec_from_8  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi

vpsrlw		$8,%ymm2,%ymm2
vpand		%ymm15,%ymm2,%ymm2
dec_from_8  2
vmovdqa		%ymm2,(%rdi)
add		    $32,%rdi

vpand		%ymm15,%ymm3,%ymm10
dec_from_8  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi

vpsrlw		$8,%ymm3,%ymm3
vpand		%ymm15,%ymm3,%ymm3
dec_from_8  3
vmovdqa		%ymm3,(%rdi)
add		    $32,%rdi

vpand		%ymm15,%ymm4,%ymm10
dec_from_8  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi

vpsrlw		$8,%ymm4,%ymm4
vpand		%ymm15,%ymm4,%ymm4
dec_from_8  4
vmovdqa		%ymm4,(%rdi)
add		    $32,%rdi

vpand		%ymm15,%ymm5,%ymm10
dec_from_8  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi

vpsrlw		$8,%ymm5,%ymm5
vpand		%ymm15,%ymm5,%ymm5
dec_from_8  5
vmovdqa		%ymm5,(%rdi)
add		    $32,%rdi

vpand		%ymm15,%ymm6,%ymm10
dec_from_8  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi

vpsrlw		$8,%ymm6,%ymm6
vpand		%ymm15,%ymm6,%ymm6
dec_from_8  6
vmovdqa		%ymm6,(%rdi)
add		    $32,%rdi

vpand		%ymm15,%ymm7,%ymm10
dec_from_8  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi

vpsrlw		$8,%ymm7,%ymm7
vpand		%ymm15,%ymm7,%ymm7
dec_from_8  7
vmovdqa		%ymm7,(%rdi)
add		    $32,%rdi

ret


.global cdecl(poly_decompress_from_9_avx2)
cdecl(poly_decompress_from_9_avx2):
.p2align 5

#const
vmovdqa		cdecl(_16xmask9)(%rip),%ymm15
vmovdqa		cdecl(_16xq)(%rip),%ymm14
vmovdqa		cdecl(_8x256)(%rip),%ymm13

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

vpand		%ymm15,%ymm0,%ymm10
dec_from_9  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi  #0

vpsrlw		$9,%ymm0,%ymm0
vpsllw		$7,%ymm1,%ymm10
vpaddw		%ymm10,%ymm0,%ymm0
vpand		%ymm15,%ymm0,%ymm0
dec_from_9  0
vmovdqa		%ymm0,(%rdi)
add		    $32,%rdi  #1

vpsrlw		$2,%ymm1,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_9  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi  #2

vpsrlw		$11,%ymm1,%ymm1
vpsllw		$5,%ymm2,%ymm10
vpaddw		%ymm10,%ymm1,%ymm1
vpand		%ymm15,%ymm1,%ymm1
dec_from_9  1
vmovdqa		%ymm1,(%rdi)
add		    $32,%rdi  #3

vpsrlw		$4,%ymm2,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_9  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi  #4

vpsrlw		$13,%ymm2,%ymm2
vpsllw		$3,%ymm3,%ymm10
vpaddw		%ymm10,%ymm2,%ymm2
vpand		%ymm15,%ymm2,%ymm2
dec_from_9  2
vmovdqa		%ymm2,(%rdi)
add		    $32,%rdi  #5

vpsrlw		$6,%ymm3,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_9  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi  #6

vpsrlw		$15,%ymm3,%ymm3
vpsllw		$1,%ymm4,%ymm10
vpaddw		%ymm10,%ymm3,%ymm3
vpand		%ymm15,%ymm3,%ymm3
dec_from_9  3
vmovdqa		%ymm3,(%rdi)
add		    $32,%rdi  #7

vpsrlw		$8,%ymm4,%ymm4
vpsllw		$8,%ymm5,%ymm10
vpaddw		%ymm10,%ymm4,%ymm4
vpand		%ymm15,%ymm4,%ymm4
dec_from_9  4
vmovdqa		%ymm4,(%rdi)
add		    $32,%rdi  #8

vpsrlw		$1,%ymm5,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_9  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi  #9

vpsrlw		$10,%ymm5,%ymm5
vpsllw		$6,%ymm6,%ymm10
vpaddw		%ymm10,%ymm5,%ymm5
vpand		%ymm15,%ymm5,%ymm5
dec_from_9  5
vmovdqa		%ymm5,(%rdi)
add		    $32,%rdi  #10

vpsrlw		$3,%ymm6,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_9  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi  #11

vpsrlw		$12,%ymm6,%ymm6
vpsllw		$4,%ymm7,%ymm10
vpaddw		%ymm10,%ymm6,%ymm6
vpand		%ymm15,%ymm6,%ymm6
dec_from_9  6
vmovdqa		%ymm6,(%rdi)
add		    $32,%rdi  #12

vpsrlw		$5,%ymm7,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_9  10
vmovdqa		%ymm10,(%rdi)
add		    $32,%rdi  #13

vpsrlw		$14,%ymm7,%ymm7
vpsllw		$2,%ymm8,%ymm10
vpaddw		%ymm10,%ymm7,%ymm7
vpand		%ymm15,%ymm7,%ymm7
dec_from_9  7
vmovdqa		%ymm7,(%rdi)
add		    $32,%rdi  #14

vpsrlw		$7,%ymm8,%ymm10
vpand		%ymm15,%ymm10,%ymm10
dec_from_9  10
vmovdqa		%ymm10,(%rdi)  #15

ret

.global cdecl(poly_decompress_from_10_avx2)
cdecl(poly_decompress_from_10_avx2):
.p2align 5

#const
vmovdqa		cdecl(_16xmask10)(%rip),%ymm15
vmovdqa		cdecl(_16xq)(%rip),%ymm14
vmovdqa		cdecl(_8x512)(%rip),%ymm13

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


.global cdecl(poly_frombytes_avx2)
cdecl(poly_frombytes_avx2):
.p2align 5

#const
vmovdqa		cdecl(_16xmask)(%rip),%ymm15

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
vmovdqu		320(%rsi),%ymm10
vmovdqu		352(%rsi),%ymm11
vmovdqu		384(%rsi),%ymm12

vpand		%ymm15,%ymm0,%ymm13
vmovdqa		%ymm13,(%rdi)
add		$32,%rdi

vpsrlw		$13,%ymm0,%ymm0
vpsllw		$3,%ymm1,%ymm13
vpaddw		%ymm13,%ymm0,%ymm0
vpand		%ymm15,%ymm0,%ymm0
vmovdqa		%ymm0,(%rdi)
add		$32,%rdi

vpsrlw		$10,%ymm1,%ymm1
vpsllw		$6,%ymm2,%ymm13
vpaddw		%ymm13,%ymm1,%ymm1
vpand		%ymm15,%ymm1,%ymm1
vmovdqa		%ymm1,(%rdi)
add		$32,%rdi

vpsrlw		$7,%ymm2,%ymm2
vpsllw		$9,%ymm3,%ymm13
vpaddw		%ymm13,%ymm2,%ymm2
vpand		%ymm15,%ymm2,%ymm2
vmovdqa		%ymm2,(%rdi)
add		$32,%rdi

vpsrlw		$4,%ymm3,%ymm3
vpsllw		$12,%ymm4,%ymm13
vpaddw		%ymm13,%ymm3,%ymm3
vpand		%ymm15,%ymm3,%ymm3
vmovdqa		%ymm3,(%rdi)
add		$32,%rdi

vpsrlw		$1,%ymm4,%ymm4
vpand		%ymm15,%ymm4,%ymm13
vmovdqa		%ymm13,(%rdi)
add		$32,%rdi

vpsrlw		$13,%ymm4,%ymm4
vpsllw		$2,%ymm5,%ymm13
vpaddw		%ymm13,%ymm4,%ymm4
vpand		%ymm15,%ymm4,%ymm4
vmovdqa		%ymm4,(%rdi)
add		$32,%rdi

vpsrlw		$11,%ymm5,%ymm5
vpsllw		$5,%ymm6,%ymm13
vpaddw		%ymm13,%ymm5,%ymm5
vpand		%ymm15,%ymm5,%ymm5
vmovdqa		%ymm5,(%rdi)
add		$32,%rdi

vpsrlw		$8,%ymm6,%ymm6
vpsllw		$8,%ymm7,%ymm13
vpaddw		%ymm13,%ymm6,%ymm6
vpand		%ymm15,%ymm6,%ymm6
vmovdqa		%ymm6,(%rdi)
add		$32,%rdi

vpsrlw		$5,%ymm7,%ymm7
vpsllw		$11,%ymm8,%ymm13
vpaddw		%ymm13,%ymm7,%ymm7
vpand		%ymm15,%ymm7,%ymm7
vmovdqa		%ymm7,(%rdi)
add		$32,%rdi

vpsrlw		$2,%ymm8,%ymm8
vpand		%ymm15,%ymm8,%ymm13
vmovdqa		%ymm13,(%rdi)
add		$32,%rdi

vpsrlw		$13,%ymm8,%ymm8
vpsllw		$1,%ymm9,%ymm13
vpaddw		%ymm13,%ymm8,%ymm8
vpand		%ymm15,%ymm8,%ymm8
vmovdqa		%ymm8,(%rdi)
add		$32,%rdi

vpsrlw		$12,%ymm9,%ymm9
vpsllw		$4,%ymm10,%ymm13
vpaddw		%ymm13,%ymm9,%ymm9
vpand		%ymm15,%ymm9,%ymm9
vmovdqa		%ymm9,(%rdi)
add		$32,%rdi

vpsrlw		$9,%ymm10,%ymm10
vpsllw		$7,%ymm11,%ymm13
vpaddw		%ymm13,%ymm10,%ymm10
vpand		%ymm15,%ymm10,%ymm10
vmovdqa		%ymm10,(%rdi)
add		$32,%rdi

vpsrlw		$6,%ymm11,%ymm11
vpsllw		$10,%ymm12,%ymm13
vpaddw		%ymm13,%ymm11,%ymm11
vpand		%ymm15,%ymm11,%ymm11
vmovdqa		%ymm11,(%rdi)
add		$32,%rdi

vpsrlw		$3,%ymm12,%ymm12
vmovdqa		%ymm12,(%rdi)

ret
