#include "consts.h"

.macro compress_to_10		a,b,a1,b1

vextracti128	$1,%ymm\a,%xmm\b

vpmovzxwd		%xmm\a,%ymm\a
vpslld			$0xa,%ymm\a,%ymm\a
vpaddd			%ymm0,%ymm\a,%ymm\a
vpshufd			$0xf5,%ymm\a,%ymm\a1
vpmuludq		%ymm1,%ymm\a,%ymm\a
vpmuludq		%ymm1,%ymm\a1,%ymm\a1
vpshufd			$0xf5,%ymm\a,%ymm\a
vpblendd		$0xaa,%ymm\a1,%ymm\a,%ymm\a
vpshufb			%ymm2,%ymm\a,%ymm\a
vpermq			$0xf8,%ymm\a,%ymm\a

vpmovzxwd		%xmm\b,%ymm\b
vpslld			$0xa,%ymm\b,%ymm\b
vpaddd			%ymm0,%ymm\b,%ymm\b
vpshufd			$0xf5,%ymm\b,%ymm\b1
vpmuludq		%ymm1,%ymm\b,%ymm\b
vpmuludq		%ymm1,%ymm\b1,%ymm\b1
vpshufd			$0xf5,%ymm\b,%ymm\b
vpblendd		$0xaa,%ymm\b1,%ymm\b,%ymm\b
vpshufb			%ymm2,%ymm\b,%ymm\b
vpermq			$0x8f,%ymm\b,%ymm\b

vpblendd		$0xf0,%ymm\b,%ymm\a,%ymm\a
vpsrlw			$6,%ymm\a,%ymm\a

.endm

.macro compress_to_11		a,b,a1,b1

vextracti128	$1,%ymm\a,%xmm\b

vpmovzxwd		%xmm\a,%ymm\a
vpslld			$0xb,%ymm\a,%ymm\a
vpaddd			%ymm0,%ymm\a,%ymm\a
vpshufd			$0xf5,%ymm\a,%ymm\a1
vpmuludq		%ymm1,%ymm\a,%ymm\a
vpmuludq		%ymm1,%ymm\a1,%ymm\a1
vpshufd			$0xf5,%ymm\a,%ymm\a
vpblendd		$0xaa,%ymm\a1,%ymm\a,%ymm\a
vpslld          $15,%ymm\a,%ymm\a
vpsrld          $21,%ymm\a,%ymm\a
vpshufb			%ymm2,%ymm\a,%ymm\a
vpermq			$0xf8,%ymm\a,%ymm\a

vpmovzxwd		%xmm\b,%ymm\b
vpslld			$0xb,%ymm\b,%ymm\b
vpaddd			%ymm0,%ymm\b,%ymm\b
vpshufd			$0xf5,%ymm\b,%ymm\b1
vpmuludq		%ymm1,%ymm\b,%ymm\b
vpmuludq		%ymm1,%ymm\b1,%ymm\b1
vpshufd			$0xf5,%ymm\b,%ymm\b
vpblendd		$0xaa,%ymm\b1,%ymm\b,%ymm\b
vpslld          $15,%ymm\b,%ymm\b
vpsrld          $21,%ymm\b,%ymm\b
vpshufb			%ymm2,%ymm\b,%ymm\b
vpermq			$0x8f,%ymm\b,%ymm\b

vpblendd		$0xf0,%ymm\b,%ymm\a,%ymm\a

.endm

.text

.global cdecl(poly_compress_to_10_avx2)
cdecl(poly_compress_to_10_avx2):
#consts
vmovdqa		_8XHALFQF*2(%rdx),%ymm0
vmovdqa		_8XMAGIC*2(%rdx),%ymm1
vmovdqa		_32XMASKSHUF*2(%rdx),%ymm2

#load
vmovdqa		(%rsi),%ymm3
vmovdqa		32(%rsi),%ymm4
vmovdqa		64(%rsi),%ymm5
vmovdqa		96(%rsi),%ymm6
vmovdqa		128(%rsi),%ymm7
vmovdqa		160(%rsi),%ymm8
vmovdqa		192(%rsi),%ymm9
vmovdqa		224(%rsi),%ymm10

compress_to_10		3,11,12,13
compress_to_10		4,14,15,11
compress_to_10		5,12,13,14
compress_to_10		6,15,11,12
compress_to_10		7,13,14,15
compress_to_10		8,11,12,13
compress_to_10		9,14,15,11
compress_to_10		10,12,13,14

# use ymm15 as tmp, result store in ymm3
vpsllw		$10,%ymm4,%ymm15
vpaddw		%ymm15,%ymm3,%ymm3
vmovdqu		%ymm3,(%rdi)
add		$32,%rdi

# use ymm15 as tmp, result store in ymm4
vpsrlw		$6,%ymm4,%ymm4
vpsllw		$4,%ymm5,%ymm15
vpaddw		%ymm15,%ymm4,%ymm4
vpsllw		$14,%ymm6,%ymm15
vpaddw		%ymm15,%ymm4,%ymm4
vmovdqu		%ymm4,(%rdi)
add		$32,%rdi

# use ymm15 as tmp, result store in ymm6
vpsrlw		$2,%ymm6,%ymm6
vpsllw		$8,%ymm7,%ymm15
vpaddw		%ymm15,%ymm6,%ymm6
vmovdqu		%ymm6,(%rdi)
add		$32,%rdi

# use ymm15 as tmp, result store in ymm7
vpsrlw		$8,%ymm7,%ymm7
vpsllw		$2,%ymm8,%ymm15
vpaddw		%ymm15,%ymm7,%ymm7
vpsllw		$12,%ymm9,%ymm15
vpaddw		%ymm15,%ymm7,%ymm7
vmovdqu		%ymm7,(%rdi)
add		$32,%rdi

# use ymm15 as tmp, result store in ymm9
vpsrlw		$4,%ymm9,%ymm9
vpsllw		$6,%ymm10,%ymm15
vpaddw		%ymm15,%ymm9,%ymm9
vmovdqu		%ymm9,(%rdi)
add		$32,%rdi

vmovdqa		256(%rsi),%ymm11
vmovdqa		288(%rsi),%ymm12
vmovdqa		320(%rsi),%ymm13
vmovdqa		352(%rsi),%ymm14
vmovdqa		384(%rsi),%ymm15
vmovdqa		416(%rsi),%ymm3
vmovdqa		448(%rsi),%ymm4
vmovdqa		480(%rsi),%ymm5

compress_to_10		11,6,7,8
compress_to_10		12,9,10,6
compress_to_10		13,7,8,9
compress_to_10		14,10,6,7
compress_to_10		15,8,9,10
compress_to_10		3,6,7,8
compress_to_10		4,9,10,6
compress_to_10		5,7,8,9

# use ymm6 as tmp, result store in ymm11
vpsllw		$10,%ymm12,%ymm6
vpaddw		%ymm6,%ymm11,%ymm11
vmovdqu		%ymm11,(%rdi)
add		$32,%rdi

# use ymm6 as tmp, result store in ymm12
vpsrlw		$6,%ymm12,%ymm12
vpsllw		$4,%ymm13,%ymm6
vpaddw		%ymm6,%ymm12,%ymm12
vpsllw		$14,%ymm14,%ymm6
vpaddw		%ymm6,%ymm12,%ymm12
vmovdqu		%ymm12,(%rdi)
add		$32,%rdi

# use ymm6 as tmp, result store in ymm14
vpsrlw		$2,%ymm14,%ymm14
vpsllw		$8,%ymm15,%ymm6
vpaddw		%ymm6,%ymm14,%ymm14
vmovdqu		%ymm14,(%rdi)
add		$32,%rdi

# use ymm6 as tmp, result store in ymm15
vpsrlw		$8,%ymm15,%ymm15
vpsllw		$2,%ymm3,%ymm6
vpaddw		%ymm6,%ymm15,%ymm15
vpsllw		$12,%ymm4,%ymm6
vpaddw		%ymm6,%ymm15,%ymm15
vmovdqu		%ymm15,(%rdi)
add		$32,%rdi

# use ymm6 as tmp, result store in ymm4
vpsrlw		$4,%ymm4,%ymm4
vpsllw		$6,%ymm5,%ymm6
vpaddw		%ymm6,%ymm4,%ymm4
vmovdqu		%ymm4,(%rdi)

ret

.global cdecl(poly_compress_to_11_avx2)
cdecl(poly_compress_to_11_avx2):
#consts
vmovdqa		_8XHALFQF*2(%rdx),%ymm0
vmovdqa		_8XMAGIC*2(%rdx),%ymm1
vmovdqa		_32XMASKSHUF*2(%rdx),%ymm2

#load
vmovdqa		(%rsi),%ymm3
vmovdqa		32(%rsi),%ymm4
vmovdqa		64(%rsi),%ymm5
vmovdqa		96(%rsi),%ymm6
vmovdqa		128(%rsi),%ymm7
vmovdqa		160(%rsi),%ymm8
vmovdqa		192(%rsi),%ymm9
vmovdqa		224(%rsi),%ymm10

compress_to_11		3,11,12,13
compress_to_11		4,14,15,11
compress_to_11		5,12,13,14
compress_to_11		6,15,11,12
compress_to_11		7,13,14,15
compress_to_11		8,11,12,13
compress_to_11		9,14,15,11
compress_to_11		10,12,13,14

# use ymm15 as tmp, result store in ymm3
vpsllw		$11,%ymm4,%ymm15
vpaddw		%ymm15,%ymm3,%ymm3
vmovdqu		%ymm3,(%rdi)
add		$32,%rdi

# use ymm15 as tmp, result store in ymm4
vpsrlw		$5,%ymm4,%ymm4
vpsllw		$6,%ymm5,%ymm15
vpaddw		%ymm15,%ymm4,%ymm4
vmovdqu		%ymm4,(%rdi)
add		$32,%rdi


# use ymm15 as tmp, result store in ymm5
vpsrlw		$10,%ymm5,%ymm5
vpsllw		$1,%ymm6,%ymm15
vpaddw		%ymm15,%ymm5,%ymm5
vpsllw		$12,%ymm7,%ymm15
vpaddw		%ymm15,%ymm5,%ymm5
vmovdqu		%ymm5,(%rdi)
add		$32,%rdi

# use ymm15 as tmp, result store in ymm7
vpsrlw		$4,%ymm7,%ymm7
vpsllw		$7,%ymm8,%ymm15
vpaddw		%ymm15,%ymm7,%ymm7
vmovdqu		%ymm7,(%rdi)
add		$32,%rdi

# use ymm15 as tmp, result store in ymm8
vpsrlw		$9,%ymm8,%ymm8
vpsllw		$2,%ymm9,%ymm15
vpaddw		%ymm15,%ymm8,%ymm8
vpsllw		$13,%ymm10,%ymm15
vpaddw		%ymm15,%ymm8,%ymm8
vmovdqu		%ymm8,(%rdi)
add		$32,%rdi

vmovdqa		256(%rsi),%ymm11
vmovdqa		288(%rsi),%ymm12
vmovdqa		320(%rsi),%ymm13
vmovdqa		352(%rsi),%ymm14
vmovdqa		384(%rsi),%ymm15
vmovdqa		416(%rsi),%ymm3
vmovdqa		448(%rsi),%ymm4
vmovdqa		480(%rsi),%ymm5

compress_to_11		11,6,7,8
compress_to_11		12,9,6,7
compress_to_11		13,8,9,6
compress_to_11		14,7,8,9
compress_to_11		15,6,7,8
compress_to_11		3,9,6,7
compress_to_11		4,8,9,6
compress_to_11		5,7,8,9

# use ymm6 as tmp, result store in ymm10
vpsrlw      $3,%ymm10,%ymm10
vpsllw		$8,%ymm11,%ymm6
vpaddw		%ymm6,%ymm10,%ymm10
vmovdqu		%ymm10,(%rdi)
add		$32,%rdi

# use ymm6 as tmp, result store in ymm11
vpsrlw		$8,%ymm11,%ymm11
vpsllw		$3,%ymm12,%ymm6
vpaddw		%ymm6,%ymm11,%ymm11
vpsllw		$14,%ymm13,%ymm6
vpaddw		%ymm6,%ymm11,%ymm11
vmovdqu		%ymm11,(%rdi)
add		$32,%rdi

# use ymm6 as tmp, result store in ymm13
vpsrlw		$2,%ymm13,%ymm13
vpsllw		$9,%ymm14,%ymm6
vpaddw		%ymm6,%ymm13,%ymm13
vmovdqu		%ymm13,(%rdi)
add		$32,%rdi

# use ymm6 as tmp, result store in ymm14
vpsrlw		$7,%ymm14,%ymm14
vpsllw		$4,%ymm15,%ymm6
vpaddw		%ymm6,%ymm14,%ymm14
vpsllw		$15,%ymm3,%ymm6
vpaddw		%ymm6,%ymm14,%ymm14
vmovdqu		%ymm14,(%rdi)
add		$32,%rdi

# use ymm6 as tmp, result store in ymm3
vpsrlw		$1,%ymm3,%ymm3
vpsllw		$10,%ymm4,%ymm6
vpaddw		%ymm6,%ymm3,%ymm3
vmovdqu		%ymm3,(%rdi)
add		$32,%rdi

# use ymm6 as tmp, result store in ymm4
vpsrlw		$6,%ymm4,%ymm4
vpsllw		$5,%ymm5,%ymm6
vpaddw		%ymm6,%ymm4,%ymm4
vmovdqu		%ymm4,(%rdi)

ret
