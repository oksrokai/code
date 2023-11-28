#include "macros.i"
.syntax unified
.cpu cortex-m4
.thumb

.macro plant_red q, qa, qinv, tmp
	mul \tmp, \tmp, \qinv     
	//tmp*qinv mod 2^2n/ 2^n; in high half
	smlatt \tmp, \tmp, \q, \qa
	// result in high half
.endm

.macro single_basemul poly0, poly1, tmp, tmp2, q, qa, qinv, zeta
    //z*g0[1]
    smulwt \tmp, \zeta, \poly1 
    //z*g0[1]*q+qa
    smlabt \tmp, \tmp, \q, \qa  
    //z*f0[1]*g0[1]
    smultt \tmp, \poly0, \tmp  
    //f0[0]*g0[0]+z*f0[1]*g0[1]
    smlabb \tmp, \poly0, \poly1, \tmp 
    //f0[0]*g0[0]+z*f0[1]*g0[1] -> r0
    plant_red \q, \qa, \qinv, \tmp
    //f0[1]*g0[0]+f0[0]*g0[1] -> r1
    smuadx \tmp2, \poly0, \poly1 
    plant_red \q, \qa, \qinv, \tmp2
    // r0 r1 all upper half
    pkhtb \tmp, \tmp2, \tmp, asr#16
.endm

.macro single_basemul_yhat poly0, tmp, tmp2, q, qa, qinv, zeta
    //z*g0[1]
    smulwt \tmp, \zeta, \poly0 
    //z*g0[1]*q+qa
    smlabt \tmp, \tmp, \q, \qa
    //r0
    //plant_red \q, \qa, \qinv, \tmp
    //r1
    lsl \tmp2, \poly0, #16
    // r0 r1 all upper half
    pkhtb \tmp, \tmp2, \tmp, asr#16
.endm

.global multi_basemul_asm
.type multi_basemul_asm, %function
.align 2
multi_basemul_asm:
    push {r4-r11, lr}
    rptr    .req r0
	aptr    .req r1
	bptr    .req r2
	zetaptr .req r3
	poly0   .req r4
	poly1   .req r6
	poly2   .req r5
	poly3   .req r7
	q       .req r8
	qa      .req r14
	qinv    .req r9
	tmp     .req r10
	tmp2    .req r11
	zeta    .req r12

    movt q, #7681
	### qinv=0x8b83e201
	movw qinv, #57857
	movt qinv, #35715
    movw qa, #61448
    .rept 64
    ldr poly0, [aptr], #4
    ldr poly1, [bptr], #4
    ldr.w poly2, [aptr, #508]
    ldr.w poly3, [bptr, #508]
    ldr zeta, [zetaptr], #4

    single_basemul poly0,poly1,tmp,tmp2,q,qa,qinv,zeta
    uadd16 poly0, poly0, poly2
    uadd16 poly1, poly1, poly3
    single_basemul poly2,poly3,tmp2,poly2,q,qa,qinv,zeta
    single_basemul poly0,poly1,poly2,poly3,q,qa,qinv,zeta
    single_basemul_yhat tmp2,poly0,poly1,q,qa,qinv,zeta
    uadd16 poly1, tmp, poly0
    usub16 poly2, poly2, tmp
    usub16 poly2, poly2, tmp2
    movw poly3, #8737
    doublebarrett poly2,tmp,tmp2,q,poly3
    doublebarrett poly1,tmp,tmp2,q,poly3
    str poly1, [rptr], #4
    str.w poly2, [rptr, #508]

    ldr poly0, [aptr], #4
    ldr poly1, [bptr], #4
    ldr.w poly2, [aptr, #508]
    ldr.w poly3, [bptr, #508]
    neg zeta, zeta
    single_basemul poly0,poly1,tmp,tmp2,q,qa,qinv,zeta
    uadd16 poly0, poly0, poly2
    uadd16 poly1, poly1, poly3
    single_basemul poly2,poly3,tmp2,poly2,q,qa,qinv,zeta
    single_basemul poly0,poly1,poly2,poly3,q,qa,qinv,zeta
    single_basemul_yhat tmp2,poly0,poly1,q,qa,qinv,zeta
    uadd16 poly1, tmp, poly0
    usub16 poly2, poly2, tmp
    usub16 poly2, poly2, tmp2
    movw poly3, #8737
    doublebarrett poly2,tmp,tmp2,q,poly3
    doublebarrett poly1,tmp,tmp2,q,poly3
    str poly1, [rptr], #4
    str.w poly2, [rptr, #508]
    .endr
    pop {r4-r11, pc}