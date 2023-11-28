.syntax unified
.cpu cortex-m4
.thumb

.macro doublecon shift,stmp
  lsr polyk, mtmp, \shift
  bfi polyk, polyk, #15, #2
  and polyk, #0X00010001

  pkhbt pack, polyk, \stmp, lsl#16
  smuad pack, pack, mult
  mla pack, pack, magic, ad
  lsr pack, pack, #26
  and result1, pack, #15

  pkhtb pack, \stmp, polyk, asr#16
  smuad pack, pack, mult
  mla pack, pack, magic, ad
  lsr pack, pack, #22
  and pack, pack, #0Xf0

  orr result1, result1, pack
.endm

.align 2
.global con
.type con, %function
con:
  push {r3-r12, lr}

  polyc       .req r0
  msg         .req r1
  polys       .req r2

  mtmp        .req r3
  stmp0       .req r4
  stmp1       .req r5
  polyk       .req r6
  mult        .req r7
  pack        .req r8
  magic       .req r9
  ad          .req r10
  result1     .req r11
  result2     .req r12

  movw magic, #20159
  movw result2, #0
  ldr mult, =#1075216
  ldr ad, =#33544576
  // load m[0-8] and store in c[0-3]
  ldr mtmp, [msg]
  ldr stmp0, [polys], #4 
  ldr stmp1, [polys], #4
  doublecon #0,stmp0
  orr result2, result1
  doublecon #2,stmp1
  lsl result1, result1, #8
  orr result2, result1

  ldr stmp0, [polys], #4
  ldr stmp1, [polys], #4
  doublecon #4,stmp0
  lsl result1, result1, #16
  orr result2, result1
  doublecon #6,stmp1
  lsl result1, result1, #24
  orr result2, result1

  str result2, [polyc]

  // load m[9-15] and store in c[32-35]
  movw result2, #0

  ldr stmp0, [polys], #4
  ldr stmp1, [polys], #4
  doublecon #8,stmp0
  orr result2, result1
  doublecon #10,stmp1
  lsl result1, result1, #8
  orr result2, result1

  ldr stmp0, [polys], #4
  ldr stmp1, [polys], #4
  doublecon #12,stmp0
  lsl result1, result1, #16
  orr result2, result1
  doublecon #14,stmp1
  lsl result1, result1, #24
  orr result2, result1

  str result2, [polyc, 32]

  // load m[16-23] and store in c[64-67]
  movw result2, #0

  ldr stmp0, [polys], #4
  ldr stmp1, [polys], #4
  doublecon #16,stmp0
  orr result2, result1
  doublecon #18,stmp1
  lsl result1, result1, #8
  orr result2, result1

  ldr stmp0, [polys], #4
  ldr stmp1, [polys], #4
  doublecon #20,stmp0
  lsl result1, result1, #16
  orr result2, result1
  doublecon #22,stmp1
  lsl result1, result1, #24
  orr result2, result1

  str result2, [polyc, 64]

  // load m[24-31] and store in c[96-99]
  movw result2, #0

  ldr stmp0, [polys], #4
  ldr stmp1, [polys], #4
  doublecon #24,stmp0
  orr result2, result1
  doublecon #26,stmp1
  lsl result1, result1, #8
  orr result2, result1

  ldr stmp0, [polys], #4
  ldr stmp1, [polys], #4
  doublecon #28,stmp0
  lsl result1, result1, #16
  orr result2, result1
  doublecon #30,stmp1
  lsl result1, result1, #24
  orr result2, result1

  str result2, [polyc, 96]

  pop {r3-r12, pc}