.syntax unified
.cpu cortex-m4
.thumb

.macro compute_t shift t_position
	ldr stmp, [poly_s], #4

	lsr t_tmp, msg, \shift
	and t_tmp, t_tmp, #1

	mul t_tmp, t_tmp, mult_53280
	mla t_tmp, stmp, mult_32, t_tmp

  	mla t_tmp, t_tmp, magic, ad
  	lsr t_tmp, t_tmp, #26
  	and \t_position, t_tmp, #31 
.endm

.align 2
.global con5
.type con5, %function
con5:
	push {r3-r12, lr}

  	poly_c       .req r0
  	msg          .req r1 
  	poly_s       .req r2

  	magic        .req r3

  	mult_32         .req r4
  	mult_53280         .req r12
  	ad           .req r5

  	tmp1         .req r6
  	tmp2         .req r7
  	tmp3         .req r8
	stmp         .req r9
	t_tmp         .req r10

	ci_32        .req r11 

  	movw magic, #20159
  	//ldr mult, =#2150432
  	movw mult_32, #32
  	movw mult_53280, #53280
  	ldr ad, =#33544576

  	/// t[0-6] -> c[0,1,2,3]
  	movw ci_32, #0
  	// c[0]
  	compute_t #0 tmp1 // t0
  	compute_t #1 tmp2 // t1
  	lsl tmp3, tmp2, #5 // t1<<5
  	orr tmp1, tmp3
  	orr ci_32, tmp1
	// c[1]
  	lsr tmp1, tmp2, #3 // t1>>3
  	compute_t #2 tmp2
  	lsl tmp2, tmp2, #2 // t2<<2
  	orr tmp1, tmp2 // t1>>3 || t2<<2
  	compute_t #3 tmp2 // t3
  	lsl tmp3, tmp2, #7 // t3<<7
  	orr tmp1, tmp3 // t1>>3 || t2<<2 || t3<<7
  	lsl tmp1, tmp1, #8
  	orr ci_32, tmp1
	// c[2]
  	lsr tmp1, tmp2, #1 // t3>>1
  	compute_t #4 tmp2 // t4
  	lsl tmp3, tmp2, #4 // t4<<4
  	orr tmp1, tmp3
  	lsl tmp1, tmp1, #16
  	orr ci_32, tmp1
  	// c[3]
  	lsr tmp1, tmp2, #4 // t4>>4
  	compute_t #5 tmp2 // t5
  	lsl tmp2, tmp2, #1 // t5<<1
  	orr tmp1, tmp2 // t4>>4 || t5<<1
  	compute_t #6 tmp2 // t6
  	lsl tmp3, tmp2, #6 // t6<<6
  	orr tmp1, tmp3
  	lsl tmp1, tmp1, #24
  	orr ci_32, tmp1

  	str ci_32, [poly_c] 

  	// t[6-12] -> c[32,33,34,35]
  	movw ci_32, #0 
  	// c[32]
  	lsr tmp1, tmp2, #2 // t6>>2
  	compute_t #7 tmp2 // t7
  	lsl tmp2, tmp2, #3 // t7<<3
  	orr tmp1, tmp2
  	orr ci_32, tmp1
  	// c[33]
  	compute_t #8 tmp1 // t8
  	compute_t #9 tmp2 // t9
  	lsl tmp3, tmp2, #5 // t9<<5
  	orr tmp1, tmp3
  	lsl tmp1, tmp1, #8
  	orr ci_32, tmp1
	// c[34]
  	lsr tmp1, tmp2, #3 // t9>>3
  	compute_t #10 tmp2
  	lsl tmp2, tmp2, #2 // t10<<2
  	orr tmp1, tmp2 // t9>>3 || t10<<2
  	compute_t #11 tmp2 // t11
  	lsl tmp3, tmp2, #7 // t11<<7
  	orr tmp1, tmp3 // t9>>3 || t10<<2 || t11<<7
  	lsl tmp1, tmp1, #16
  	orr ci_32, tmp1
	// c[35]
  	lsr tmp1, tmp2, #1 // t11>>1
  	compute_t #12 tmp2 // t12
  	lsl tmp3, tmp2, #4 // t12<<4
  	orr tmp1, tmp3
  	lsl tmp1, tmp1, #24
  	orr ci_32, tmp1

  	str ci_32, [poly_c, 32]

  	// t[12-19] -> c[64-67]
  	movw ci_32, #0 
  	// c[64]
  	lsr tmp1, tmp2, #4 // t12>>4
  	compute_t #13 tmp2 // t13
  	lsl tmp2, tmp2, #1 // t13<<1
  	orr tmp1, tmp2 // t12>>4 || t13<<1
  	compute_t #14 tmp2 // t14
  	lsl tmp3, tmp2, #6 // t14<<6
  	orr tmp1, tmp3
  	orr ci_32, tmp1
  	// c[65]
  	lsr tmp1, tmp2, #2 // t14>>2
  	compute_t #15 tmp2 // t15
  	lsl tmp2, tmp2, #3 // t15<<3
  	orr tmp1, tmp2
  	lsl tmp1, tmp1, #8
  	orr ci_32, tmp1
  	// c[66]
  	compute_t #16 tmp1 // t16
  	compute_t #17 tmp2 // t17
  	lsl tmp3, tmp2, #5 // t17<<5
  	orr tmp1, tmp3
  	lsl tmp1, tmp1, #16
  	orr ci_32, tmp1
	// c[67]
  	lsr tmp1, tmp2, #3 // t17>>3
  	compute_t #18 tmp2
  	lsl tmp2, tmp2, #2 // t18<<2
  	orr tmp1, tmp2 // t17>>3 || t18<<2
  	compute_t #19 tmp2 // t19
  	lsl tmp3, tmp2, #7 // t19<<7
  	orr tmp1, tmp3 // t17>>3 || t18<<2 || t19<<7
  	lsl tmp1, tmp1, #24
  	orr ci_32, tmp1

  	str ci_32, [poly_c, 64] 

  	// t[19-25] -> c[96-99]
  	movw ci_32, #0 
  	// c[96]
  	lsr tmp1, tmp2, #1 // t19>>1
  	compute_t #20 tmp2 // t20
  	lsl tmp3, tmp2, #4 // t20<<4
  	orr tmp1, tmp3
  	orr ci_32, tmp1
  	// c[97]
  	lsr tmp1, tmp2, #4 // t20>>4
  	compute_t #21 tmp2 // t21
  	lsl tmp2, tmp2, #1 // t21<<1
  	orr tmp1, tmp2 // t20>>4 || t21<<1
  	compute_t #22 tmp2 // t22
  	lsl tmp3, tmp2, #6 // t22<<6
  	orr tmp1, tmp3
  	lsl tmp1, tmp1, #8
  	orr ci_32, tmp1
  	// c[98]
  	lsr tmp1, tmp2, #2 // t22>>2
  	compute_t #23 tmp2 // t23
  	lsl tmp2, tmp2, #3 // t23<<3
  	orr tmp1, tmp2
  	lsl tmp1, tmp1, #16
  	orr ci_32, tmp1
  	// c[99]
  	compute_t #24 tmp1 // t24
  	compute_t #25 tmp2 // t25
  	lsl tmp3, tmp2, #5 // t25<<5
  	orr tmp1, tmp3
  	lsl tmp1, tmp1, #24
  	orr ci_32, tmp1

  	str ci_32, [poly_c, 96] 

  	// t[25-31] -> c[128-131]
  	movw ci_32, #0 
  	// c[128]
  	lsr tmp1, tmp2, #3 // t25>>3
  	compute_t #26 tmp2
  	lsl tmp2, tmp2, #2 // t26<<2
  	orr tmp1, tmp2 // t25>>3 || t26<<2
  	compute_t #27 tmp2 // t27
  	lsl tmp3, tmp2, #7 // t27<<7
  	orr tmp1, tmp3 // t25>>3 || t26<<2 || t27<<7
  	orr ci_32, tmp1
	// c[129]
  	lsr tmp1, tmp2, #1 // t27>>1
  	compute_t #28 tmp2 // t28
  	lsl tmp3, tmp2, #4 // t28<<4
  	orr tmp1, tmp3
  	lsl tmp1, tmp1, #8
  	orr ci_32, tmp1
  	// c[130]
  	lsr tmp1, tmp2, #4 // t28>>4
  	compute_t #29 tmp2 // t29
  	lsl tmp2, tmp2, #1 // t29<<1
  	orr tmp1, tmp2 // t28>>4 || t29<<1
  	compute_t #30 tmp2 // t30
  	lsl tmp3, tmp2, #6 // t30<<6
  	orr tmp1, tmp3
  	lsl tmp1, tmp1, #16
  	orr ci_32, tmp1
  	// c[131]
  	lsr tmp1, tmp2, #2 // t30>>2
  	compute_t #31 tmp2 // t31
  	lsl tmp2, tmp2, #3 // t31<<3
  	orr tmp1, tmp2
  	lsl tmp1, tmp1, #24
  	orr ci_32, tmp1

  	str ci_32, [poly_c, 128]

  	pop {r3-r12, pc}
