void BSSN_ID(REAL xx0,REAL xx1,REAL xx2,REAL Cartxyz0,REAL Cartxyz1,REAL Cartxyz2,
	REAL *hDD00,REAL *hDD01,REAL *hDD02,REAL *hDD11,REAL *hDD12,REAL *hDD22,
	REAL *aDD00,REAL *aDD01,REAL *aDD02,REAL *aDD11,REAL *aDD12,REAL *aDD22,
	REAL *trK,
	REAL *lambdaU0,REAL *lambdaU1,REAL *lambdaU2,
	REAL *vetU0,REAL *vetU1,REAL *vetU2,
	REAL *betU0,REAL *betU1,REAL *betU2,
	REAL *alpha,REAL *cf) {
   {
         const double tmp0 = cos(xx1);
         const double tmp1 = pow(tmp0, 2);
         const double tmp2 = sin(xx1);
         const double tmp3 = pow(tmp2, 2);
         const double tmp4 = pow(sin(xx2), 2);
         const double tmp5 = tmp3*tmp4;
         const double tmp6 = pow(cos(xx2), 2);
         const double tmp7 = tmp3*tmp6;
         const double tmp8 = tmp1 + tmp5 + tmp7;
         const double tmp9 = tmp0*tmp2*xx0;
         const double tmp10 = tmp4*tmp9;
         const double tmp11 = tmp6*tmp9;
         const double tmp12 = tmp10 + tmp11 - tmp9;
         const double tmp13 = pow(xx0, 2);
         const double tmp14 = 1.0/tmp13;
         const double tmp15 = tmp13*tmp3;
         const double tmp16 = tmp1*tmp13;
         const double tmp17 = tmp15 + tmp16*tmp4 + tmp16*tmp6;
         const double tmp18 = tmp13*tmp5;
         const double tmp19 = tmp13*tmp7;
         const double tmp20 = tmp18 + tmp19;
         const double tmp21 = pow(psi_in, 4);
         const double tmp22 = tmp18*tmp21 + tmp19*tmp21;
         const double tmp23 = tmp1*tmp21;
         const double tmp24 = tmp13*tmp23;
         *hDD00 = tmp8 - 1;
         *hDD01 = tmp12/xx0;
         *hDD02 = 0;
         *hDD11 = tmp14*(-tmp13 + tmp17);
         *hDD12 = 0;
         *hDD22 = tmp14*(-tmp15 + tmp20)/tmp3;
         *aDD00 = 0;
         *aDD01 = 0;
         *aDD02 = 0;
         *aDD11 = 0;
         *aDD12 = 0;
         *aDD22 = 0;
         *trK = 0;
         *lambdaU0 = 0;
         *lambdaU1 = 0;
         *lambdaU2 = 0;
         *vetU0 = 0;
         *vetU1 = 0;
         *vetU2 = 0;
         *betU0 = 0;
         *betU1 = 0;
         *betU2 = 0;
         *alpha = alpha_in;
         *cf = pow((-tmp22*pow(tmp10*tmp21 + tmp11*tmp21 - tmp21*tmp9, 2) + tmp22*(tmp15*tmp21 + tmp24*tmp4 + tmp24*tmp6)*(tmp21*tmp5 + tmp21*tmp7 + tmp23))/(-pow(tmp12, 2)*tmp20 + tmp17*tmp20*tmp8), -1.0/6.0);
   }
}
