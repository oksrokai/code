#include "consts.h"
.include "fq.inc"

.macro rec_3  x,t,m,n
# load s
vmovdqa     \m(%rdx),%ymm\t
# ymm10 reserved for reduce tmp
vpsllw      $3,%ymm\t,%ymm\t
vpaddw      %ymm14,%ymm\x,%ymm\x
vpmullw     %ymm12,%ymm\x,%ymm\x
vpsubw      %ymm\t,%ymm\x,%ymm\x
vpmulhw     %ymm13,%ymm\x,%ymm\t
vpsraw      $10,%ymm\t,%ymm\x
vpsllw      $13,%ymm\x,%ymm\x
vpsrlw      $15,%ymm\x,%ymm\x
vpsllw      $\n,%ymm\x,%ymm\x
vpor        %ymm\x,%ymm11,%ymm11
.endm

.macro rec_4  x,t,m,n
# load s
vmovdqu     \m(%rdx),%ymm\t
vpmovzxwd   %xmm\t,%ymm\t
vpslld      $4,%ymm\t,%ymm\t
vpaddd      %ymm14,%ymm\x,%ymm\x
vpmulld     %ymm12,%ymm\x,%ymm\x
vpsubd      %ymm\t,%ymm\x,%ymm\x
vpmulld     %ymm13,%ymm\x,%ymm\t
vpsrld      $29,%ymm\t,%ymm\x
vpand       %ymm0,%ymm\x,%ymm\x
vpslld      $\n,%ymm\x,%ymm\x
vpor        %ymm\x,%ymm11,%ymm11
.endm

.macro rec_5  x,t,m,n
# load s
vmovdqu     \m(%rdx),%ymm\t
vpmovzxwd   %xmm\t,%ymm\t
vpslld      $5,%ymm\t,%ymm\t
vpaddd      %ymm14,%ymm\x,%ymm\x
vpmulld     %ymm12,%ymm\x,%ymm\x
vpsubd      %ymm\t,%ymm\x,%ymm\x
vpmulld     %ymm13,%ymm\x,%ymm\t
vpsrld      $30,%ymm\t,%ymm\x
vpand       %ymm0,%ymm\x,%ymm\x
vpslld      $\n,%ymm\x,%ymm\x
vpor        %ymm\x,%ymm11,%ymm11
.endm

.text

.global cdecl(rec_avx_3)
cdecl(rec_avx_3):
#consts
vmovdqa		_16XSEVEN*2(%rcx),%ymm15
vmovdqa		_16XTWO*2(%rcx),%ymm14
vmovdqa		_16XMAGIC*2(%rcx),%ymm13
vmovdqa		_16XQ*2(%rcx),%ymm12

vpxor       %ymm11,%ymm11,%ymm11

#load
vmovdqu		(%rsi),%ymm2
vmovdqu		32(%rsi),%ymm3
vmovdqu		64(%rsi),%ymm4

vpand		%ymm15,%ymm2,%ymm5
rec_3       5,9,0,0

vpsrlw		$3,%ymm2,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_3       6,5,32,1

vpsrlw		$6,%ymm2,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_3       7,6,64,2

vpsrlw		$9,%ymm2,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_3       8,7,96,3

vpsrlw		$12,%ymm2,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_3       9,8,128,4

vpsrlw		$15,%ymm2,%ymm2
vpsllw		$1,%ymm3,%ymm5
vpaddw		%ymm2,%ymm5,%ymm5
vpand		%ymm15,%ymm5,%ymm5
rec_3       5,9,160,5

vpsrlw		$2,%ymm3,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_3       6,5,192,6

vpsrlw		$5,%ymm3,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_3       7,6,224,7

vpsrlw		$8,%ymm3,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_3       8,7,256,8

vpsrlw		$11,%ymm3,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_3       9,8,288,9

vpsrlw		$14,%ymm3,%ymm3
vpsllw		$2,%ymm4,%ymm5
vpaddw		%ymm3,%ymm5,%ymm5
vpand		%ymm15,%ymm5,%ymm5
rec_3       5,9,320,10

vpsrlw		$1,%ymm4,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_3       6,5,352,11

vpsrlw		$4,%ymm4,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_3       7,6,384,12

vpsrlw		$7,%ymm4,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_3       8,7,416,13

vpsrlw		$10,%ymm4,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_3       9,8,448,14

vpsrlw		$13,%ymm4,%ymm4
rec_3       4,3,480,15

vmovdqu     %ymm11,(%rdi)

ret

.global cdecl(rec_avx_4)
cdecl(rec_avx_4):
#consts
vmovdqa		_8XFIFTEEN*2(%rcx),%ymm15  #0x1111
vmovdqa		_8XFOUR*2(%rcx),%ymm14   #4
vmovdqa		_8XMAGICG*2(%rcx),%ymm13
vmovdqa		_8XQ*2(%rcx),%ymm12     #Q
vmovdqa		_8XONE*2(%rcx),%ymm0

vpxor       %ymm11,%ymm11,%ymm11

#load c
vmovdqu		(%rsi),%ymm1
vmovdqu		32(%rsi),%ymm2
vmovdqu		64(%rsi),%ymm3
vmovdqu		96(%rsi),%ymm4

#rec 0-63
vpand		%ymm15,%ymm1,%ymm5
rec_4       5,9,0,0

vpsrld		$4,%ymm1,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_4       6,5,16,1

vpsrld		$8,%ymm1,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_4       7,6,32,2

vpsrld		$12,%ymm1,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_4       8,7,48,3

vpsrld		$16,%ymm1,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_4       9,8,64,4

vpsrld		$20,%ymm1,%ymm5
vpand		%ymm15,%ymm5,%ymm5
rec_4       5,9,80,5

vpsrld		$24,%ymm1,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_4       6,5,96,6

vpsrld		$28,%ymm1,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_4       7,6,112,7

#rec 64-127
vpand		%ymm15,%ymm2,%ymm5
rec_4       5,9,128,8

vpsrld		$4,%ymm2,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_4       6,5,144,9

vpsrld		$8,%ymm2,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_4       7,6,160,10

vpsrld		$12,%ymm2,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_4       8,7,176,11

vpsrld		$16,%ymm2,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_4       9,8,192,12

vpsrld		$20,%ymm2,%ymm5
vpand		%ymm15,%ymm5,%ymm5
rec_4       5,9,208,13

vpsrld		$24,%ymm2,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_4       6,5,224,14

vpsrld		$28,%ymm2,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_4       7,6,240,15

#rec 128-191
vpand		%ymm15,%ymm3,%ymm5
rec_4       5,9,256,16

vpsrld		$4,%ymm3,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_4       6,5,272,17

vpsrld		$8,%ymm3,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_4       7,6,288,18

vpsrld		$12,%ymm3,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_4       8,7,304,19

vpsrld		$16,%ymm3,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_4       9,8,320,20

vpsrld		$20,%ymm3,%ymm5
vpand		%ymm15,%ymm5,%ymm5
rec_4       5,9,336,21

vpsrld		$24,%ymm3,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_4       6,5,352,22

vpsrld		$28,%ymm3,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_4       7,6,368,23

#rec 192-256
vpand		%ymm15,%ymm4,%ymm5
rec_4       5,9,384,24

vpsrld		$4,%ymm4,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_4       6,5,400,25

vpsrld		$8,%ymm4,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_4       7,6,416,26

vpsrld		$12,%ymm4,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_4       8,7,432,27

vpsrld		$16,%ymm4,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_4       9,8,448,28

vpsrld		$20,%ymm4,%ymm5
vpand		%ymm15,%ymm5,%ymm5
rec_4       5,9,464,29

vpsrld		$24,%ymm4,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_4       6,5,480,30

vpsrld		$28,%ymm4,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_4       7,6,496,31

vmovdqu     %ymm11,(%rdi)

ret

.global cdecl(rec_avx_5)
cdecl(rec_avx_5):
#consts
vmovdqa		_8XTHIRTYONE*2(%rcx),%ymm15  #0b11111
vmovdqa		_8XEIGHT*2(%rcx),%ymm14   #8
vmovdqa		_8XMAGICG*2(%rcx),%ymm13
vmovdqa		_8XQ*2(%rcx),%ymm12     #Q
vmovdqa		_8XONE*2(%rcx),%ymm0

vpxor       %ymm11,%ymm11,%ymm11

#load c
vmovdqu		(%rsi),%ymm1
vmovdqu		32(%rsi),%ymm2
vmovdqu		64(%rsi),%ymm3
vmovdqu		96(%rsi),%ymm4
vmovdqu		128(%rsi),%ymm5

#rec 0-63
vpand		%ymm15,%ymm1,%ymm6
rec_5       6,9,0,0

vpsrld		$5,%ymm1,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_5       7,6,16,1

vpsrld		$10,%ymm1,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_5       8,7,32,2

vpsrld		$15,%ymm1,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_5       9,8,48,3

vpsrld		$20,%ymm1,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_5       6,9,64,4

vpsrld		$25,%ymm1,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_5       7,6,80,5

vpsrld		$30,%ymm1,%ymm8
vpslld      $2,%ymm2,%ymm10
vpor        %ymm10,%ymm8,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_5       8,7,96,6

vpsrld		$3,%ymm2,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_5       9,8,112,7

#rec 64-127
vpsrld		$8,%ymm2,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_5       6,9,128,8

vpsrld		$13,%ymm2,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_5       7,6,144,9

vpsrld		$18,%ymm2,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_5       8,7,160,10

vpsrld		$23,%ymm2,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_5       9,8,176,11

vpsrld		$28,%ymm2,%ymm6
vpslld      $4,%ymm3,%ymm10
vpor        %ymm10,%ymm6,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_5       6,9,192,12

vpsrld		$1,%ymm3,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_5       7,6,208,13

vpsrld		$6,%ymm3,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_5       8,7,224,14

vpsrld		$11,%ymm3,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_5       9,8,240,15

#rec 128-191
vpsrld		$16,%ymm3,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_5       6,9,256,16

vpsrld		$21,%ymm3,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_5       7,6,272,17

vpsrld		$26,%ymm3,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_5       8,7,288,18

vpsrld		$31,%ymm3,%ymm9
vpslld      $1,%ymm4,%ymm10
vpor        %ymm10,%ymm9,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_5       9,8,304,19

vpsrld		$4,%ymm4,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_5       6,9,320,20

vpsrld		$9,%ymm4,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_5       7,6,336,21

vpsrld		$14,%ymm4,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_5       8,7,352,22

vpsrld		$19,%ymm4,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_5       9,8,368,23

#rec 192-256
vpsrld		$24,%ymm4,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_5       6,9,384,24

vpsrld		$29,%ymm4,%ymm7
vpslld      $3,%ymm5,%ymm10
vpor        %ymm10,%ymm7,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_5       7,6,400,25

vpsrld		$2,%ymm5,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_5       8,7,416,26

vpsrld		$7,%ymm5,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_5       9,8,432,27

vpsrld		$12,%ymm5,%ymm6
vpand		%ymm15,%ymm6,%ymm6
rec_5       6,9,448,28

vpsrld		$17,%ymm5,%ymm7
vpand		%ymm15,%ymm7,%ymm7
rec_5       7,6,464,29

vpsrld		$22,%ymm5,%ymm8
vpand		%ymm15,%ymm8,%ymm8
rec_5       8,7,480,30

vpsrld		$27,%ymm5,%ymm9
vpand		%ymm15,%ymm9,%ymm9
rec_5       9,8,496,31

vmovdqu     %ymm11,(%rdi)

ret
