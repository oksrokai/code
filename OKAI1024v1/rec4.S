.syntax unified
.cpu cortex-m4
.thumb

.macro doublerec loc,stmp
  and ttmp, ctmp, #0xf
  lsr ctmp, #4
  bfi ttmp, ctmp, #16, #4
  lsr ctmp, #4

  pkhbt pack, ttmp, \stmp, lsl#16
  smlsd pack, pack, mult, ad
  mul pack, magic
  lsr pack, #29
  bfi result, pack, \loc, #1

  pkhtb pack, \stmp, ttmp, asr#16
  smlsd pack, pack, mult, ad
  mul pack, magic
  lsr pack, #29
  bfi result, pack, \loc+1, #1
.endm

.align 2
.global rec
.type rec, %function
rec:
  push {r3-r12, lr}

  t           .req r0
  c           .req r1
  spp         .req r2

  ctmp        .req r3
  ttmp        .req r4
  stmp1       .req r5
  pack        .req r6
  mult        .req r7
  ad          .req r8
  magic       .req r9
  result      .req r10
  stmp2       .req r11

  stmp3        .req r12

  movw mult, #3329
  movt mult, #16
  movw ad, #13316
  movw magic, #20159
  movw result, #0
  
  // load c[0-3]
  ldr ctmp, [c], #32
  ldr stmp1, [spp], #4
  ldr stmp2, [spp], #4
  ldr stmp3, [spp], #4
  doublerec #0,stmp1
  doublerec #2,stmp2
  doublerec #4,stmp3
  ldr stmp1, [spp], #4
  ldr stmp2, [spp], #4
  ldr stmp3, [spp], #4
  doublerec #6,stmp1

  // load c[32-35]
  ldr ctmp, [c], #32
  ldr stmp1, [spp], #4
  doublerec #8,stmp2
  doublerec #10,stmp3
  doublerec #12,stmp1
  ldr stmp1, [spp], #4
  ldr stmp2, [spp], #4
  ldr stmp3, [spp], #4
  doublerec #14,stmp1

  // load c[64-67]
  ldr ctmp, [c], #32
  ldr stmp1, [spp], #4
  doublerec #16,stmp2
  doublerec #18,stmp3
  doublerec #20,stmp1
  ldr stmp1, [spp], #4
  ldr stmp2, [spp], #4
  ldr stmp3, [spp], #4
  doublerec #22,stmp1

  // load c[96-99]
  ldr ctmp, [c]
  ldr stmp1, [spp], #4
  doublerec #24,stmp2
  doublerec #26,stmp3
  doublerec #28,stmp1
  ldr stmp1, [spp], #4
  doublerec #30,stmp1

  str result, [t]

  pop {r3-r12, pc}