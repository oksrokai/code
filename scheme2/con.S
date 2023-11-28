#if defined(__WIN32__) || defined(__APPLE__)
#define cdecl(s) _##s
#else
#define cdecl(s) s
#endif

.macro con_3  x,t,n
vpsrlw      $\n,%ymm11,%ymm\t
vpand       %ymm14,%ymm\t,%ymm\t       #k
vpmullw     %ymm13,%ymm\t,%ymm\t
vpaddw      %ymm\t,%ymm\x,%ymm\x      #s += k
vpsllw      $3,%ymm\x,%ymm\x
vpaddw      %ymm12,%ymm\x,%ymm\x
vpmulhuw    %ymm15,%ymm\x,%ymm\x
vpsllw      $3,%ymm\x,%ymm\x
vpsrlw      $13,%ymm\x,%ymm\x
.endm

.macro con_4  x,t,n
vpsrld      $\n,%ymm11,%ymm\t
vpand       %ymm14,%ymm\t,%ymm\t       #k
vpmulld     %ymm13,%ymm\t,%ymm\t
vpaddd      %ymm\t,%ymm\x,%ymm\x      #s += k
vpslld      $4,%ymm\x,%ymm\x
vpaddd      %ymm12,%ymm\x,%ymm\x
vpmulld     %ymm15,%ymm\x,%ymm\x
vpslld      $2,%ymm\x,%ymm\x
vpsrld      $28,%ymm\x,%ymm\x
.endm

.global cdecl(con_avx_3)
cdecl(con_avx_3):
.p2align 5
# load consts
vmovdqa		cdecl(_16xhalfq)(%rip),%ymm12  #Q/2
vmovdqa		cdecl(_16xhalfqc)(%rip),%ymm13  #(Q+1)/2
vmovdqa		cdecl(_16x1)(%rip),%ymm14
vmovdqa		cdecl(_16xmagic)(%rip),%ymm15

# load m
vmovdqu     (%rsi),%ymm11

# load s
vmovdqa     (%rdx),%ymm2
vmovdqa     32(%rdx),%ymm3
vmovdqa     64(%rdx),%ymm4
vmovdqa     96(%rdx),%ymm5
vmovdqa		128(%rdx),%ymm6
vmovdqa		160(%rdx),%ymm7
con_3       2,8,0
con_3       3,9,1
con_3       4,10,2
con_3       5,8,3
con_3       6,9,4
con_3       7,10,5
# use ymm10 as tmp, result store in ymm2
vpsllw		$3,%ymm3,%ymm10
vpor		%ymm10,%ymm2,%ymm2
vpsllw      $6,%ymm4,%ymm10
vpor		%ymm10,%ymm2,%ymm2
vpsllw      $9,%ymm5,%ymm10
vpor		%ymm10,%ymm2,%ymm2
vpsllw      $12,%ymm6,%ymm10
vpor		%ymm10,%ymm2,%ymm2
vpsllw      $15,%ymm7,%ymm10
vpor		%ymm10,%ymm2,%ymm2
# store c
vmovdqu     %ymm2,(%rdi)

vmovdqa		192(%rdx),%ymm8
vmovdqa		224(%rdx),%ymm9
vmovdqa		256(%rdx),%ymm10
vmovdqa		288(%rdx),%ymm2
vmovdqa		320(%rdx),%ymm3
con_3       8,4,6
con_3       9,5,7
con_3       10,6,8
con_3       2,4,9
con_3       3,5,10
# use ymm4 as tmp, result store in ymm7
vpsrlw		$1,%ymm7,%ymm7
vpsllw		$2,%ymm8,%ymm4
vpor		%ymm4,%ymm7,%ymm7
vpsllw		$5,%ymm9,%ymm4
vpor		%ymm4,%ymm7,%ymm7
vpsllw		$8,%ymm10,%ymm4
vpor		%ymm4,%ymm7,%ymm7
vpsllw		$11,%ymm2,%ymm4
vpor		%ymm4,%ymm7,%ymm7
vpsllw		$14,%ymm3,%ymm4
vpor		%ymm4,%ymm7,%ymm7
vmovdqu     %ymm7,32(%rdi)

vmovdqa		352(%rdx),%ymm4
vmovdqa		384(%rdx),%ymm5
vmovdqa		416(%rdx),%ymm6
vmovdqa		448(%rdx),%ymm7
vmovdqa		480(%rdx),%ymm8
con_3       4,2,11
con_3       5,9,12
con_3       6,10,13
con_3       7,9,14
con_3       8,10,15
# use ymm2 as tmp, result store in ymm3
vpsrlw		$2,%ymm3,%ymm3
vpsllw		$1,%ymm4,%ymm2
vpor		%ymm2,%ymm3,%ymm3
vpsllw		$4,%ymm5,%ymm2
vpor		%ymm2,%ymm3,%ymm3
vpsllw		$7,%ymm6,%ymm2
vpor		%ymm2,%ymm3,%ymm3
vpsllw		$10,%ymm7,%ymm2
vpor		%ymm2,%ymm3,%ymm3
vpsllw		$13,%ymm8,%ymm2
vpor		%ymm2,%ymm3,%ymm3
vmovdqu     %ymm3,64(%rdi)

ret

.global cdecl(con_avx_4)
cdecl(con_avx_4):
.p2align 5
# load consts
vmovdqa		cdecl(_8xhalfq)(%rip),%ymm12  #Q/2
vmovdqa		cdecl(_8xhalfqc)(%rip),%ymm13  #(Q+1)/2
vmovdqa		cdecl(_8x1)(%rip),%ymm14
vmovdqa		cdecl(_8xmagic)(%rip),%ymm15

# load m
vmovdqu     (%rsi),%ymm11

# load s
vmovdqa     (%rdx),%ymm0
vmovdqa     32(%rdx),%ymm2
vmovdqa     64(%rdx),%ymm4
vmovdqa     96(%rdx),%ymm6
vextracti128    $1,%ymm0,%xmm1
vextracti128    $1,%ymm2,%xmm3
vextracti128    $1,%ymm4,%xmm5
vextracti128    $1,%ymm6,%xmm7
vpmovzxwd   %xmm0,%ymm0
vpmovzxwd   %xmm1,%ymm1
vpmovzxwd   %xmm2,%ymm2
vpmovzxwd   %xmm3,%ymm3
vpmovzxwd   %xmm4,%ymm4
vpmovzxwd   %xmm5,%ymm5
vpmovzxwd   %xmm6,%ymm6
vpmovzxwd   %xmm7,%ymm7

con_4       0,8,0
con_4       1,9,1
con_4       2,10,2
con_4       3,8,3
con_4       4,9,4
con_4       5,10,5
con_4       6,8,6
con_4       7,9,7

# use ymm10 as tmp, result store in ymm0
vpslld		$4,%ymm1,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $8,%ymm2,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $12,%ymm3,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $16,%ymm4,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $20,%ymm5,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $24,%ymm6,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $28,%ymm7,%ymm10
vpor		%ymm10,%ymm0,%ymm0

# store c
vmovdqu     %ymm0,(%rdi)

# load s
vmovdqa		128(%rdx),%ymm0
vmovdqa		160(%rdx),%ymm2
vmovdqa		192(%rdx),%ymm4
vmovdqa		224(%rdx),%ymm6
vextracti128    $1,%ymm0,%xmm1
vextracti128    $1,%ymm2,%xmm3
vextracti128    $1,%ymm4,%xmm5
vextracti128    $1,%ymm6,%xmm7
vpmovzxwd   %xmm0,%ymm0
vpmovzxwd   %xmm1,%ymm1
vpmovzxwd   %xmm2,%ymm2
vpmovzxwd   %xmm3,%ymm3
vpmovzxwd   %xmm4,%ymm4
vpmovzxwd   %xmm5,%ymm5
vpmovzxwd   %xmm6,%ymm6
vpmovzxwd   %xmm7,%ymm7

con_4       0,8,8
con_4       1,9,9
con_4       2,10,10
con_4       3,8,11
con_4       4,9,12
con_4       5,10,13
con_4       6,8,14
con_4       7,9,15

# use ymm10 as tmp, result store in ymm0
vpslld		$4,%ymm1,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $8,%ymm2,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $12,%ymm3,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $16,%ymm4,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $20,%ymm5,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $24,%ymm6,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $28,%ymm7,%ymm10
vpor		%ymm10,%ymm0,%ymm0

# store c
vmovdqu     %ymm0,32(%rdi)

# load s
vmovdqa		256(%rdx),%ymm0
vmovdqa		288(%rdx),%ymm2
vmovdqa		320(%rdx),%ymm4
vmovdqa		352(%rdx),%ymm6
vextracti128    $1,%ymm0,%xmm1
vextracti128    $1,%ymm2,%xmm3
vextracti128    $1,%ymm4,%xmm5
vextracti128    $1,%ymm6,%xmm7
vpmovzxwd   %xmm0,%ymm0
vpmovzxwd   %xmm1,%ymm1
vpmovzxwd   %xmm2,%ymm2
vpmovzxwd   %xmm3,%ymm3
vpmovzxwd   %xmm4,%ymm4
vpmovzxwd   %xmm5,%ymm5
vpmovzxwd   %xmm6,%ymm6
vpmovzxwd   %xmm7,%ymm7

con_4       0,8,16
con_4       1,9,17
con_4       2,10,18
con_4       3,8,19
con_4       4,9,20
con_4       5,10,21
con_4       6,8,22
con_4       7,9,23

# use ymm10 as tmp, result store in ymm0
vpslld		$4,%ymm1,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $8,%ymm2,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $12,%ymm3,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $16,%ymm4,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $20,%ymm5,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $24,%ymm6,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $28,%ymm7,%ymm10
vpor		%ymm10,%ymm0,%ymm0

# store c
vmovdqu     %ymm0,64(%rdi)

# load s
vmovdqa		384(%rdx),%ymm0
vmovdqa		416(%rdx),%ymm2
vmovdqa		448(%rdx),%ymm4
vmovdqa		480(%rdx),%ymm6
vextracti128    $1,%ymm0,%xmm1
vextracti128    $1,%ymm2,%xmm3
vextracti128    $1,%ymm4,%xmm5
vextracti128    $1,%ymm6,%xmm7
vpmovzxwd   %xmm0,%ymm0
vpmovzxwd   %xmm1,%ymm1
vpmovzxwd   %xmm2,%ymm2
vpmovzxwd   %xmm3,%ymm3
vpmovzxwd   %xmm4,%ymm4
vpmovzxwd   %xmm5,%ymm5
vpmovzxwd   %xmm6,%ymm6
vpmovzxwd   %xmm7,%ymm7

con_4       0,8,24
con_4       1,9,25
con_4       2,10,26
con_4       3,8,27
con_4       4,9,28
con_4       5,10,29
con_4       6,8,30
con_4       7,9,31

# use ymm10 as tmp, result store in ymm0
vpslld		$4,%ymm1,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $8,%ymm2,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $12,%ymm3,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $16,%ymm4,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $20,%ymm5,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $24,%ymm6,%ymm10
vpor		%ymm10,%ymm0,%ymm0
vpslld      $28,%ymm7,%ymm10
vpor		%ymm10,%ymm0,%ymm0

# store c
vmovdqu     %ymm0,96(%rdi)

ret