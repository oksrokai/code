.syntax unified
.cpu cortex-m4
.thumb

.macro compute1 loc,add_pos,mul_pos
  pkhbt pack, ttmp, stmp, lsl#16
  smlsd pack, pack, mult, \add_pos
  mul pack, \mul_pos
  lsr pack, #30
  bfi result, pack, \loc, #1

  pkhtb pack, stmp, ttmp, asr#16
  smlsd pack, pack, mult, \add_pos
  mul pack,\mul_pos
  lsr pack, #30
  bfi result, pack, \loc+1, #1
.endm

.macro load_c ctmp
  and ttmp, \ctmp, #0x1f
  lsr \ctmp, #5
  bfi ttmp, \ctmp, #16, #5
  lsr \ctmp, #5
.endm

.align 2
.global rec5
.type rec5, %function
rec5:
  push {r3-r12, lr}

  t           .req r0
  c           .req r1
  spp         .req r2

  ctmp1       .req r3
  ctmp2       .req r4
  ctmp3       .req r5
  ctmp4       .req r6
  ctmp5       .req r7
  ttmp        .req r8
  stmp        .req r9
  pack        .req r10
  mult        .req r11
  result      .req r12

  movw mult, #3329
  movt mult, #32
  movw result, #0
  movw ctmp5, #20159
  movw ctmp4, #26632

  // load c[0-3] [32-35] [64-67] [96-99] [128-131]
  ldr ctmp1, [c], #32
  ldr ctmp2, [c], #32
  ldr ctmp3, [c], #32

  //load s[0 1]
  ldr stmp, [spp], #4
  load_c ctmp1
  compute1 #0,ctmp4,ctmp5

  ldr stmp, [spp], #4
  load_c ctmp1
  compute1 #2,ctmp4,ctmp5

  ldr stmp, [spp], #4
  load_c ctmp1
  compute1 #4,ctmp4,ctmp5

  ldr stmp, [spp], #4
  and ttmp, ctmp1, #0x3
  lsr ctmp1, #2
  bfi ttmp, ctmp2, #2, #3
  lsr ctmp2, #3
  bfi ttmp, ctmp2, #16, #5
  lsr ctmp2, #5
  compute1 #6,ctmp4,ctmp5

  ldr stmp, [spp], #4
  load_c ctmp2
  compute1 #8,ctmp4,ctmp5

  ldr stmp, [spp], #4
  load_c ctmp2
  compute1 #10,ctmp4,ctmp5

  ldr stmp, [spp], #4
  and ttmp, ctmp2, #0xf
  lsr ctmp2, #4
  bfi ttmp, ctmp3, #4, #1
  lsr ctmp3, #1
  bfi ttmp, ctmp3, #16, #5
  lsr ctmp3, #5
  compute1 #12,ctmp4,ctmp5

  ldr stmp, [spp], #4
  load_c ctmp3
  compute1 #14,ctmp4,ctmp5

  ldr stmp, [spp], #4
  load_c ctmp3
  compute1 #16,ctmp4,ctmp5

  mov ctmp2,ctmp4
  ldr ctmp4, [c], #32

  ldr stmp, [spp], #4
  and ttmp, ctmp3, #0x1f
  lsr ctmp3, #5
  bfi ttmp, ctmp3, #16, #1
  lsr ctmp3, #1
  bfi ttmp, ctmp4, #17, #4
  lsr ctmp4, #4
  compute1 #18,ctmp2,ctmp5

  ldr stmp, [spp], #4
  load_c ctmp4
  compute1 #20,ctmp2,ctmp5

  ldr stmp, [spp], #4
  load_c ctmp4
  compute1 #22,ctmp2,ctmp5

  mov ctmp3,ctmp5
  ldr ctmp5, [c]
  
  ldr stmp, [spp], #4
  and ttmp, ctmp4, #0x1f
  lsr ctmp4, #5
  bfi ttmp, ctmp4, #16, #3
  lsr ctmp4, #3
  bfi ttmp, ctmp5, #19, #2
  lsr ctmp5, #2
  compute1 #24,ctmp2,ctmp3

  ldr stmp, [spp], #4
  load_c ctmp5
  compute1 #26,ctmp2,ctmp3

  ldr stmp, [spp], #4
  load_c ctmp5
  compute1 #28,ctmp2,ctmp3

  ldr stmp, [spp], #4
  load_c ctmp5
  compute1 #30,ctmp2,ctmp3

  str result, [t]

  pop {r3-r12, pc}