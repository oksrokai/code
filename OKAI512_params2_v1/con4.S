.syntax unified
.cpu cortex-m4
.thumb

.macro doublecon shift
  lsr polyk, mtmp, \shift
  bfi polyk, polyk, #15, #2
  and polyk, #0X00010001

  pkhbt pack, polyk, stmp0, lsl#16
  smuad pack, pack, mult
  mla pack, pack, magic, ad
  lsr pack, pack, #26
  and result1, pack, #15

  pkhtb pack, stmp0, polyk, asr#16
  smuad pack, pack, mult
  mla pack, pack, magic, ad
  lsr pack, pack, #22
  and pack, pack, #0Xf0

  orr result1, result1, pack
.endm

.align 2
.global con4
.type con4, %function
con4:
  push {r3-r11, lr}

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

  movw magic, #8737 @ a/q == a * m >> k, magic = m
  ldr mult, =#1110032 @ 
  ldr ad, =#33550080 @ ad = m * Math.floor(q/2)

  // load m[0-8] and store in c[0-3]
  ldr mtmp, [msg]
  movw result2, #0

  ldr stmp0, [polys], #4 
  doublecon #0
  orr result2, result1

  ldr stmp0, [polys], #4
  doublecon #2
  lsl result1, result1, #8
  orr result2, result1

  ldr stmp0, [polys], #4
  doublecon #4
  lsl result1, result1, #16
  orr result2, result1

  ldr stmp0, [polys], #4
  doublecon #6
  lsl result1, result1, #24
  orr result2, result1

  str result2, [polyc]

  // load m[9-15] and store in c[32-35]
  movw result2, #0

  ldr stmp0, [polys], #4
  doublecon #8
  orr result2, result1

  ldr stmp0, [polys], #4
  doublecon #10
  lsl result1, result1, #8
  orr result2, result1

  ldr stmp0, [polys], #4
  doublecon #12
  lsl result1, result1, #16
  orr result2, result1

  ldr stmp0, [polys], #4
  doublecon #14
  lsl result1, result1, #24
  orr result2, result1

  str result2, [polyc, 32]

  // load m[16-23] and store in c[64-67]
  movw result2, #0

  ldr stmp0, [polys], #4
  doublecon #16
  orr result2, result1

  ldr stmp0, [polys], #4
  doublecon #18
  lsl result1, result1, #8
  orr result2, result1

  ldr stmp0, [polys], #4
  doublecon #20
  lsl result1, result1, #16
  orr result2, result1

  ldr stmp0, [polys], #4
  doublecon #22
  lsl result1, result1, #24
  orr result2, result1

  str result2, [polyc, 64]

  // load m[24-31] and store in c[96-99]
  movw result2, #0

  ldr stmp0, [polys], #4
  doublecon #24
  orr result2, result1

  ldr stmp0, [polys], #4
  doublecon #26
  lsl result1, result1, #8
  orr result2, result1

  ldr stmp0, [polys], #4
  doublecon #28
  lsl result1, result1, #16
  orr result2, result1

  ldr stmp0, [polys], #4
  doublecon #30
  lsl result1, result1, #24
  orr result2, result1

  str result2, [polyc, 96]

  pop {r3-r11, pc}
