.syntax unified
.cpu cortex-m4
.thumb

.macro shift_a tmp0,tmp1,tmp2,tmp3,restmp
  and \restmp,\tmp0,#0xff
  lsr \tmp0,#8
  bfi \restmp,\tmp1,#8,#4
  lsr \tmp1,#4
  bfi \restmp,\tmp2,#16,#4
  lsr \tmp2,#4
  bfi \restmp,\tmp3,#20,#8
  lsr \tmp3,#8
.endm

.align 2
.global from_bytes
.type from_bytes, %function
from_bytes:
  push {r2-r8, lr}

  rcoeffs      .req r0
  a_val        .req r1

  atmp0         .req r2
  atmp1         .req r3
  atmp2         .req r4
  res0         .req r5
  res1         .req r6
  res2         .req r7
  res3         .req r8

  ldr atmp0,[a_val],#4
  ldr atmp1,[a_val],#4
  ldr atmp2,[a_val],#4
  
  shift_a atmp0,atmp0,atmp0,atmp0,res0
  shift_a atmp0,atmp1,atmp1,atmp1,res1
  shift_a atmp1,atmp1,atmp1,atmp2,res2
  shift_a atmp2,atmp2,atmp2,atmp2,res3

  str res0,[rcoeffs],#4
  str res1,[rcoeffs],#4
  str res2,[rcoeffs],#4
  str res3,[rcoeffs],#4

  pop {r2-r8,pc}