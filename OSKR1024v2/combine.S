.syntax unified
.cpu cortex-m4
.thumb

.global combine_asm
.type combine_asm, %function
.align 2
combine_asm:
    push {r4-r11, lr}
    rptr    .req r0
	aptr    .req r1
	bptr    .req r2
	poly0   .req r3
	poly1   .req r4
	poly2   .req r5
	poly3   .req r6
	tmp0    .req r7
	tmp1    .req r8
	tmp2    .req r9
	tmp3    .req r10
    .rept 64
    ldr poly0, [aptr], #4
    ldr poly1, [bptr], #4
    ldr poly2, [aptr], #4
    ldr poly3, [bptr], #4

    pkhbt tmp0, poly0, poly1, lsl#16
    pkhtb tmp1, poly1, poly0, asr#16

    pkhbt tmp2, poly2, poly3, lsl#16
    pkhtb tmp3, poly3, poly2, asr#16

    str tmp0, [rptr], #4
    str tmp1, [rptr], #4
    str tmp2, [rptr], #4
    str tmp3, [rptr], #4
    .endr
    pop {r4-r11, pc}