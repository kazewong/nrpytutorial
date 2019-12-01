void enforce_detgammabar_constraint(const int Nxx_plus_2NGHOSTS[3],REAL *xx[3], REAL *in_gfs) {

#pragma omp parallel for
for(int i2=0; i2<Nxx_plus_2NGHOSTS[2]; i2++) {
    const REAL xx2 = xx[2][i2];
    for(int i1=0; i1<Nxx_plus_2NGHOSTS[1]; i1++) {
        const REAL xx1 = xx[1][i1];
        for(int i0=0; i0<Nxx_plus_2NGHOSTS[0]; i0++) {
            const REAL xx0 = xx[0][i0];
            /* 
             * NRPy+ Finite Difference Code Generation, Step 1 of 1: Read from main memory and compute finite difference stencils:
             */
            const double hDD00 = in_gfs[IDX4(HDD00GF, i0,i1,i2)];
            const double hDD01 = in_gfs[IDX4(HDD01GF, i0,i1,i2)];
            const double hDD02 = in_gfs[IDX4(HDD02GF, i0,i1,i2)];
            const double hDD11 = in_gfs[IDX4(HDD11GF, i0,i1,i2)];
            const double hDD12 = in_gfs[IDX4(HDD12GF, i0,i1,i2)];
            const double hDD22 = in_gfs[IDX4(HDD22GF, i0,i1,i2)];
            /* 
             * NRPy+ Finite Difference Code Generation, Step 2 of 1: Evaluate SymPy expressions and write to main memory:
             */
            const double tmp0 = hDD00 + 1;
            const double tmp1 = sin(xx1);
            const double tmp2 = pow(tmp1, 2);
            const double tmp3 = tmp2*pow(xx0, 4);
            const double tmp4 = pow(xx0, 2);
            const double tmp5 = hDD11*tmp4 + tmp4;
            const double tmp6 = tmp2*tmp4;
            const double tmp7 = hDD22*tmp6 + tmp6;
            const double tmp8 = cbrt(1.0/(-pow(hDD01, 2)*tmp4*tmp7 + 2*hDD01*hDD02*hDD12*tmp3 - pow(hDD02, 2)*tmp5*tmp6 - pow(hDD12, 2)*tmp0*tmp3 + tmp0*tmp5*tmp7))*pow(fabs(tmp1), 2.0/3.0)*pow(fabs(xx0), 4.0/3.0);
            in_gfs[IDX4(HDD00GF, i0, i1, i2)] = tmp0*tmp8 - 1;
            in_gfs[IDX4(HDD01GF, i0, i1, i2)] = hDD01*tmp8;
            in_gfs[IDX4(HDD02GF, i0, i1, i2)] = hDD02*tmp8;
            in_gfs[IDX4(HDD11GF, i0, i1, i2)] = tmp8*(hDD11 + 1) - 1;
            in_gfs[IDX4(HDD12GF, i0, i1, i2)] = hDD12*tmp8;
            in_gfs[IDX4(HDD22GF, i0, i1, i2)] = tmp8*(hDD22 + 1) - 1;
            
            
        } // END LOOP: for(int i0=0; i0<Nxx_plus_2NGHOSTS[0]; i0++)
    } // END LOOP: for(int i1=0; i1<Nxx_plus_2NGHOSTS[1]; i1++)
} // END LOOP: for(int i2=0; i2<Nxx_plus_2NGHOSTS[2]; i2++)
}
