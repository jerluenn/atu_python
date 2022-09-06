/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_ ## ID
#endif

#include <math.h>

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int int
#endif

/* Add prefix to internal symbols */
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_s2 CASADI_PREFIX(s2)
#define casadi_s3 CASADI_PREFIX(s3)
#define casadi_s4 CASADI_PREFIX(s4)
#define casadi_s5 CASADI_PREFIX(s5)
#define casadi_sq CASADI_PREFIX(sq)

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

casadi_real casadi_sq(casadi_real x) { return x*x;}

static const casadi_int casadi_s0[21] = {17, 1, 0, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
static const casadi_int casadi_s1[5] = {1, 1, 0, 1, 0};
static const casadi_int casadi_s2[3] = {0, 0, 0};
static const casadi_int casadi_s3[95] = {17, 17, 0, 0, 0, 0, 12, 23, 34, 46, 49, 52, 55, 61, 67, 73, 74, 74, 75, 75, 0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 16, 0, 1, 3, 4, 5, 6, 10, 11, 12, 13, 16, 0, 1, 3, 4, 5, 6, 10, 11, 12, 13, 16, 0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 16, 11, 12, 16, 10, 12, 16, 10, 11, 16, 3, 4, 5, 6, 13, 16, 3, 4, 5, 6, 13, 16, 3, 4, 5, 6, 13, 16, 13, 13};
static const casadi_int casadi_s4[37] = {17, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
static const casadi_int casadi_s5[4] = {17, 1, 0, 0};

/* tetherunit_integrator_impl_dae_fun_jac_x_xdot_u:(i0[17],i1[17],i2,i3[],i4[])->(o0[17],o1[17x17,75nz],o2[17x17,17nz],o3[17x1,0nz]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a3, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a4, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a5, a50, a51, a52, a53, a54, a55, a56, a57, a58, a59, a6, a60, a61, a62, a63, a64, a65, a66, a67, a68, a69, a7, a70, a71, a72, a8, a9;
  a0=2.;
  a1=arg[0]? arg[0][4] : 0;
  a2=arg[0]? arg[0][6] : 0;
  a3=(a1*a2);
  a4=arg[0]? arg[0][3] : 0;
  a5=arg[0]? arg[0][5] : 0;
  a6=(a4*a5);
  a3=(a3+a6);
  a3=(a0*a3);
  a6=arg[1]? arg[1][0] : 0;
  a6=(a3-a6);
  if (res[0]!=0) res[0][0]=a6;
  a6=(a5*a2);
  a7=(a4*a1);
  a6=(a6-a7);
  a6=(a0*a6);
  a7=arg[1]? arg[1][1] : 0;
  a7=(a6-a7);
  if (res[0]!=0) res[0][1]=a7;
  a7=casadi_sq(a4);
  a8=casadi_sq(a2);
  a7=(a7+a8);
  a7=(a0*a7);
  a8=1.;
  a7=(a7-a8);
  a9=arg[1]? arg[1][2] : 0;
  a9=(a7-a9);
  if (res[0]!=0) res[0][2]=a9;
  a9=1.0000000000000001e-01;
  a10=casadi_sq(a4);
  a11=casadi_sq(a1);
  a10=(a10+a11);
  a11=casadi_sq(a5);
  a10=(a10+a11);
  a11=casadi_sq(a2);
  a10=(a10+a11);
  a10=(a8-a10);
  a10=(a9*a10);
  a11=(a10*a4);
  a12=5.0000000000000000e-01;
  a13=4.4570729321922933e+01;
  a14=casadi_sq(a4);
  a15=casadi_sq(a1);
  a14=(a14+a15);
  a14=(a0*a14);
  a14=(a14-a8);
  a15=(a13*a14);
  a16=arg[0]? arg[0][10] : 0;
  a17=(a15*a16);
  a18=(a1*a5);
  a19=(a4*a2);
  a18=(a18+a19);
  a18=(a0*a18);
  a19=(a13*a18);
  a20=arg[0]? arg[0][11] : 0;
  a21=(a19*a20);
  a17=(a17+a21);
  a21=(a1*a2);
  a22=(a4*a5);
  a21=(a21-a22);
  a21=(a0*a21);
  a22=(a13*a21);
  a23=arg[0]? arg[0][12] : 0;
  a24=(a22*a23);
  a17=(a17+a24);
  a24=(a17*a1);
  a25=(a1*a5);
  a26=(a4*a2);
  a25=(a25-a26);
  a25=(a0*a25);
  a26=(a13*a25);
  a27=(a26*a16);
  a28=casadi_sq(a4);
  a29=casadi_sq(a5);
  a28=(a28+a29);
  a28=(a0*a28);
  a28=(a28-a8);
  a29=(a13*a28);
  a30=(a29*a20);
  a27=(a27+a30);
  a30=(a5*a2);
  a31=(a4*a1);
  a30=(a30+a31);
  a30=(a0*a30);
  a31=(a13*a30);
  a32=(a31*a23);
  a27=(a27+a32);
  a32=(a27*a5);
  a24=(a24+a32);
  a32=5.3484875186307519e+01;
  a33=(a32*a3);
  a34=(a33*a16);
  a35=(a32*a6);
  a36=(a35*a20);
  a34=(a34+a36);
  a36=(a32*a7);
  a37=(a36*a23);
  a34=(a34+a37);
  a37=(a34*a2);
  a24=(a24+a37);
  a24=(a12*a24);
  a11=(a11-a24);
  a24=arg[1]? arg[1][3] : 0;
  a11=(a11-a24);
  if (res[0]!=0) res[0][3]=a11;
  a11=(a17*a4);
  a24=(a34*a5);
  a11=(a11+a24);
  a24=(a27*a2);
  a11=(a11-a24);
  a11=(a12*a11);
  a24=(a10*a1);
  a11=(a11+a24);
  a24=arg[1]? arg[1][4] : 0;
  a11=(a11-a24);
  if (res[0]!=0) res[0][4]=a11;
  a11=(a27*a4);
  a24=(a34*a1);
  a11=(a11-a24);
  a24=(a17*a2);
  a11=(a11+a24);
  a11=(a12*a11);
  a24=(a10*a5);
  a11=(a11+a24);
  a24=arg[1]? arg[1][5] : 0;
  a11=(a11-a24);
  if (res[0]!=0) res[0][5]=a11;
  a11=(a34*a4);
  a24=(a27*a1);
  a11=(a11+a24);
  a24=(a17*a5);
  a11=(a11-a24);
  a11=(a12*a11);
  a24=(a10*a2);
  a11=(a11+a24);
  a24=arg[1]? arg[1][6] : 0;
  a11=(a11-a24);
  if (res[0]!=0) res[0][6]=a11;
  a11=-3.4329750000000003e-01;
  a24=arg[1]? arg[1][7] : 0;
  a11=(a11-a24);
  if (res[0]!=0) res[0][7]=a11;
  a11=arg[1]? arg[1][8] : 0;
  a11=(-a11);
  if (res[0]!=0) res[0][8]=a11;
  a11=arg[1]? arg[1][9] : 0;
  a11=(-a11);
  if (res[0]!=0) res[0][9]=a11;
  a11=arg[0]? arg[0][9] : 0;
  a24=(a6*a11);
  a37=arg[0]? arg[0][8] : 0;
  a38=(a7*a37);
  a24=(a24-a38);
  a38=arg[1]? arg[1][10] : 0;
  a24=(a24+a38);
  a24=(-a24);
  if (res[0]!=0) res[0][10]=a24;
  a24=arg[0]? arg[0][7] : 0;
  a38=(a7*a24);
  a39=(a3*a11);
  a38=(a38-a39);
  a39=arg[1]? arg[1][11] : 0;
  a38=(a38+a39);
  a38=(-a38);
  if (res[0]!=0) res[0][11]=a38;
  a38=(a3*a37);
  a39=(a6*a24);
  a38=(a38-a39);
  a39=arg[1]? arg[1][12] : 0;
  a38=(a38+a39);
  a38=(-a38);
  if (res[0]!=0) res[0][12]=a38;
  a38=arg[0]? arg[0][15] : 0;
  a39=arg[0]? arg[0][13] : 0;
  a40=(a38*a39);
  a41=casadi_sq(a17);
  a42=casadi_sq(a27);
  a41=(a41+a42);
  a42=casadi_sq(a34);
  a41=(a41+a42);
  a41=sqrt(a41);
  a42=(a40*a41);
  a43=arg[1]? arg[1][13] : 0;
  a42=(a42+a43);
  a42=(-a42);
  if (res[0]!=0) res[0][13]=a42;
  a42=arg[1]? arg[1][14] : 0;
  a8=(a8-a42);
  if (res[0]!=0) res[0][14]=a8;
  a8=arg[1]? arg[1][15] : 0;
  a8=(-a8);
  if (res[0]!=0) res[0][15]=a8;
  a8=-4.4570729321922933e+01;
  a42=1.8696874518574300e-02;
  a43=(a42*a27);
  a44=(a43*a34);
  a45=2.2436249422289160e-02;
  a46=(a45*a34);
  a47=(a46*a27);
  a44=(a44-a47);
  a47=(a25*a24);
  a48=(a28*a37);
  a47=(a47+a48);
  a48=(a30*a11);
  a47=(a47+a48);
  a44=(a44-a47);
  a44=(a8*a44);
  a47=casadi_sq(a44);
  a48=(a45*a34);
  a49=(a48*a17);
  a50=(a42*a17);
  a51=(a50*a34);
  a49=(a49-a51);
  a51=(a14*a24);
  a52=(a18*a37);
  a51=(a51+a52);
  a52=(a21*a11);
  a51=(a51+a52);
  a49=(a49+a51);
  a49=(a8*a49);
  a51=casadi_sq(a49);
  a47=(a47+a51);
  a51=-5.3484875186307519e+01;
  a52=(a45*a17);
  a53=(a52*a27);
  a54=(a45*a27);
  a55=(a54*a17);
  a53=(a53-a55);
  a53=(a51*a53);
  a55=casadi_sq(a53);
  a47=(a47+a55);
  a47=sqrt(a47);
  a55=arg[1]? arg[1][16] : 0;
  a55=(a47-a55);
  if (res[0]!=0) res[0][16]=a55;
  a55=(a0*a5);
  if (res[1]!=0) res[1][0]=a55;
  a56=(a0*a1);
  a57=(-a56);
  if (res[1]!=0) res[1][1]=a57;
  a57=(a4+a4);
  a57=(a0*a57);
  if (res[1]!=0) res[1][2]=a57;
  a58=(a4+a4);
  a58=(a9*a58);
  a59=(a4*a58);
  a59=(a10-a59);
  a60=(a4+a4);
  a60=(a0*a60);
  a61=(a13*a60);
  a61=(a16*a61);
  a62=(a0*a2);
  a63=(a13*a62);
  a63=(a20*a63);
  a61=(a61+a63);
  a63=(a0*a5);
  a64=(a13*a63);
  a64=(a23*a64);
  a61=(a61-a64);
  a64=(a1*a61);
  a65=(a4+a4);
  a65=(a0*a65);
  a66=(a13*a65);
  a66=(a20*a66);
  a67=(a0*a2);
  a68=(a13*a67);
  a68=(a16*a68);
  a66=(a66-a68);
  a68=(a0*a1);
  a69=(a13*a68);
  a69=(a23*a69);
  a66=(a66+a69);
  a69=(a5*a66);
  a64=(a64+a69);
  a69=(a32*a55);
  a69=(a16*a69);
  a70=(a32*a56);
  a70=(a20*a70);
  a69=(a69-a70);
  a70=(a32*a57);
  a70=(a23*a70);
  a69=(a69+a70);
  a70=(a2*a69);
  a64=(a64+a70);
  a64=(a12*a64);
  a59=(a59-a64);
  if (res[1]!=0) res[1][3]=a59;
  a59=(a4*a61);
  a59=(a59+a17);
  a64=(a5*a69);
  a59=(a59+a64);
  a64=(a2*a66);
  a59=(a59-a64);
  a59=(a12*a59);
  a64=(a1*a58);
  a59=(a59-a64);
  if (res[1]!=0) res[1][4]=a59;
  a59=(a4*a66);
  a59=(a59+a27);
  a64=(a1*a69);
  a59=(a59-a64);
  a64=(a2*a61);
  a59=(a59+a64);
  a59=(a12*a59);
  a64=(a5*a58);
  a59=(a59-a64);
  if (res[1]!=0) res[1][5]=a59;
  a59=(a4*a69);
  a59=(a59+a34);
  a64=(a1*a66);
  a59=(a59+a64);
  a64=(a5*a61);
  a59=(a59-a64);
  a59=(a12*a59);
  a58=(a2*a58);
  a59=(a59-a58);
  if (res[1]!=0) res[1][6]=a59;
  a59=(a11*a56);
  a58=(a37*a57);
  a59=(a59+a58);
  if (res[1]!=0) res[1][7]=a59;
  a57=(a24*a57);
  a59=(a11*a55);
  a57=(a57-a59);
  a57=(-a57);
  if (res[1]!=0) res[1][8]=a57;
  a55=(a37*a55);
  a56=(a24*a56);
  a55=(a55+a56);
  a55=(-a55);
  if (res[1]!=0) res[1][9]=a55;
  a55=(a17+a17);
  a56=(a55*a61);
  a57=(a27+a27);
  a59=(a57*a66);
  a56=(a56+a59);
  a59=(a34+a34);
  a58=(a59*a69);
  a56=(a56+a58);
  a58=(a41+a41);
  a56=(a56/a58);
  a56=(a40*a56);
  a56=(-a56);
  if (res[1]!=0) res[1][10]=a56;
  a44=(a44+a44);
  a56=(a42*a66);
  a56=(a34*a56);
  a64=(a43*a69);
  a56=(a56+a64);
  a64=(a45*a69);
  a64=(a27*a64);
  a70=(a46*a66);
  a64=(a64+a70);
  a56=(a56-a64);
  a65=(a37*a65);
  a67=(a24*a67);
  a65=(a65-a67);
  a68=(a11*a68);
  a65=(a65+a68);
  a56=(a56-a65);
  a56=(a8*a56);
  a56=(a44*a56);
  a49=(a49+a49);
  a65=(a45*a69);
  a65=(a17*a65);
  a68=(a48*a61);
  a65=(a65+a68);
  a68=(a42*a61);
  a68=(a34*a68);
  a69=(a50*a69);
  a68=(a68+a69);
  a65=(a65-a68);
  a60=(a24*a60);
  a62=(a37*a62);
  a60=(a60+a62);
  a63=(a11*a63);
  a60=(a60-a63);
  a65=(a65+a60);
  a65=(a8*a65);
  a65=(a49*a65);
  a56=(a56+a65);
  a53=(a53+a53);
  a65=(a45*a61);
  a65=(a27*a65);
  a60=(a52*a66);
  a65=(a65+a60);
  a66=(a45*a66);
  a66=(a17*a66);
  a61=(a54*a61);
  a66=(a66+a61);
  a65=(a65-a66);
  a65=(a51*a65);
  a65=(a53*a65);
  a56=(a56+a65);
  a47=(a47+a47);
  a56=(a56/a47);
  if (res[1]!=0) res[1][11]=a56;
  a56=(a0*a2);
  if (res[1]!=0) res[1][12]=a56;
  a65=(a0*a4);
  a66=(-a65);
  if (res[1]!=0) res[1][13]=a66;
  a66=(a1+a1);
  a66=(a9*a66);
  a61=(a4*a66);
  a60=(a1+a1);
  a60=(a0*a60);
  a63=(a13*a60);
  a63=(a16*a63);
  a62=(a0*a5);
  a68=(a13*a62);
  a68=(a20*a68);
  a63=(a63+a68);
  a68=(a0*a2);
  a69=(a13*a68);
  a69=(a23*a69);
  a63=(a63+a69);
  a69=(a1*a63);
  a69=(a69+a17);
  a67=(a0*a5);
  a64=(a13*a67);
  a64=(a16*a64);
  a70=(a0*a4);
  a71=(a13*a70);
  a71=(a23*a71);
  a64=(a64+a71);
  a71=(a5*a64);
  a69=(a69+a71);
  a71=(a32*a56);
  a71=(a16*a71);
  a72=(a32*a65);
  a72=(a20*a72);
  a71=(a71-a72);
  a72=(a2*a71);
  a69=(a69+a72);
  a69=(a12*a69);
  a61=(a61+a69);
  a61=(-a61);
  if (res[1]!=0) res[1][14]=a61;
  a61=(a4*a63);
  a69=(a5*a71);
  a61=(a61+a69);
  a69=(a2*a64);
  a61=(a61-a69);
  a61=(a12*a61);
  a69=(a1*a66);
  a69=(a10-a69);
  a61=(a61+a69);
  if (res[1]!=0) res[1][15]=a61;
  a61=(a4*a64);
  a69=(a1*a71);
  a69=(a69+a34);
  a61=(a61-a69);
  a69=(a2*a63);
  a61=(a61+a69);
  a61=(a12*a61);
  a69=(a5*a66);
  a61=(a61-a69);
  if (res[1]!=0) res[1][16]=a61;
  a61=(a4*a71);
  a69=(a1*a64);
  a69=(a69+a27);
  a61=(a61+a69);
  a69=(a5*a63);
  a61=(a61-a69);
  a61=(a12*a61);
  a66=(a2*a66);
  a61=(a61-a66);
  if (res[1]!=0) res[1][17]=a61;
  a61=(a11*a65);
  if (res[1]!=0) res[1][18]=a61;
  a61=(a11*a56);
  if (res[1]!=0) res[1][19]=a61;
  a56=(a37*a56);
  a65=(a24*a65);
  a56=(a56+a65);
  a56=(-a56);
  if (res[1]!=0) res[1][20]=a56;
  a56=(a55*a63);
  a65=(a57*a64);
  a56=(a56+a65);
  a65=(a59*a71);
  a56=(a56+a65);
  a56=(a56/a58);
  a56=(a40*a56);
  a56=(-a56);
  if (res[1]!=0) res[1][21]=a56;
  a56=(a42*a64);
  a56=(a34*a56);
  a65=(a43*a71);
  a56=(a56+a65);
  a65=(a45*a71);
  a65=(a27*a65);
  a61=(a46*a64);
  a65=(a65+a61);
  a56=(a56-a65);
  a67=(a24*a67);
  a70=(a11*a70);
  a67=(a67+a70);
  a56=(a56-a67);
  a56=(a8*a56);
  a56=(a44*a56);
  a67=(a45*a71);
  a67=(a17*a67);
  a70=(a48*a63);
  a67=(a67+a70);
  a70=(a42*a63);
  a70=(a34*a70);
  a71=(a50*a71);
  a70=(a70+a71);
  a67=(a67-a70);
  a60=(a24*a60);
  a62=(a37*a62);
  a60=(a60+a62);
  a68=(a11*a68);
  a60=(a60+a68);
  a67=(a67+a60);
  a67=(a8*a67);
  a67=(a49*a67);
  a56=(a56+a67);
  a67=(a45*a63);
  a67=(a27*a67);
  a60=(a52*a64);
  a67=(a67+a60);
  a64=(a45*a64);
  a64=(a17*a64);
  a63=(a54*a63);
  a64=(a64+a63);
  a67=(a67-a64);
  a67=(a51*a67);
  a67=(a53*a67);
  a56=(a56+a67);
  a56=(a56/a47);
  if (res[1]!=0) res[1][22]=a56;
  a56=(a0*a4);
  if (res[1]!=0) res[1][23]=a56;
  a67=(a0*a2);
  if (res[1]!=0) res[1][24]=a67;
  a64=(a5+a5);
  a64=(a9*a64);
  a63=(a4*a64);
  a60=(a0*a1);
  a68=(a13*a60);
  a68=(a20*a68);
  a62=(a0*a4);
  a70=(a13*a62);
  a70=(a23*a70);
  a68=(a68-a70);
  a70=(a1*a68);
  a71=(a0*a1);
  a65=(a13*a71);
  a65=(a16*a65);
  a61=(a5+a5);
  a61=(a0*a61);
  a66=(a13*a61);
  a66=(a20*a66);
  a65=(a65+a66);
  a66=(a0*a2);
  a69=(a13*a66);
  a69=(a23*a69);
  a65=(a65+a69);
  a69=(a5*a65);
  a69=(a69+a27);
  a70=(a70+a69);
  a69=(a32*a56);
  a69=(a16*a69);
  a72=(a32*a67);
  a72=(a20*a72);
  a69=(a69+a72);
  a72=(a2*a69);
  a70=(a70+a72);
  a70=(a12*a70);
  a63=(a63+a70);
  a63=(-a63);
  if (res[1]!=0) res[1][25]=a63;
  a63=(a4*a68);
  a70=(a5*a69);
  a70=(a70+a34);
  a63=(a63+a70);
  a70=(a2*a65);
  a63=(a63-a70);
  a63=(a12*a63);
  a70=(a1*a64);
  a63=(a63-a70);
  if (res[1]!=0) res[1][26]=a63;
  a63=(a4*a65);
  a70=(a1*a69);
  a63=(a63-a70);
  a70=(a2*a68);
  a63=(a63+a70);
  a63=(a12*a63);
  a70=(a5*a64);
  a70=(a10-a70);
  a63=(a63+a70);
  if (res[1]!=0) res[1][27]=a63;
  a63=(a4*a69);
  a70=(a1*a65);
  a63=(a63+a70);
  a70=(a5*a68);
  a70=(a70+a17);
  a63=(a63-a70);
  a63=(a12*a63);
  a64=(a2*a64);
  a63=(a63-a64);
  if (res[1]!=0) res[1][28]=a63;
  a63=(a11*a67);
  a63=(-a63);
  if (res[1]!=0) res[1][29]=a63;
  a63=(a11*a56);
  if (res[1]!=0) res[1][30]=a63;
  a56=(a37*a56);
  a67=(a24*a67);
  a56=(a56-a67);
  a56=(-a56);
  if (res[1]!=0) res[1][31]=a56;
  a56=(a55*a68);
  a67=(a57*a65);
  a56=(a56+a67);
  a67=(a59*a69);
  a56=(a56+a67);
  a56=(a56/a58);
  a56=(a40*a56);
  a56=(-a56);
  if (res[1]!=0) res[1][32]=a56;
  a56=(a42*a65);
  a56=(a34*a56);
  a67=(a43*a69);
  a56=(a56+a67);
  a67=(a45*a69);
  a67=(a27*a67);
  a63=(a46*a65);
  a67=(a67+a63);
  a56=(a56-a67);
  a71=(a24*a71);
  a61=(a37*a61);
  a71=(a71+a61);
  a66=(a11*a66);
  a71=(a71+a66);
  a56=(a56-a71);
  a56=(a8*a56);
  a56=(a44*a56);
  a71=(a45*a69);
  a71=(a17*a71);
  a66=(a48*a68);
  a71=(a71+a66);
  a66=(a42*a68);
  a66=(a34*a66);
  a69=(a50*a69);
  a66=(a66+a69);
  a71=(a71-a66);
  a60=(a37*a60);
  a62=(a11*a62);
  a60=(a60-a62);
  a71=(a71+a60);
  a71=(a8*a71);
  a71=(a49*a71);
  a56=(a56+a71);
  a71=(a45*a68);
  a71=(a27*a71);
  a60=(a52*a65);
  a71=(a71+a60);
  a65=(a45*a65);
  a65=(a17*a65);
  a68=(a54*a68);
  a65=(a65+a68);
  a71=(a71-a65);
  a71=(a51*a71);
  a71=(a53*a71);
  a56=(a56+a71);
  a56=(a56/a47);
  if (res[1]!=0) res[1][33]=a56;
  a56=(a0*a1);
  if (res[1]!=0) res[1][34]=a56;
  a71=(a0*a5);
  if (res[1]!=0) res[1][35]=a71;
  a65=(a2+a2);
  a65=(a0*a65);
  if (res[1]!=0) res[1][36]=a65;
  a68=(a2+a2);
  a9=(a9*a68);
  a68=(a4*a9);
  a60=(a0*a4);
  a62=(a13*a60);
  a62=(a20*a62);
  a66=(a0*a1);
  a69=(a13*a66);
  a69=(a23*a69);
  a62=(a62+a69);
  a69=(a1*a62);
  a61=(a0*a5);
  a67=(a13*a61);
  a67=(a23*a67);
  a0=(a0*a4);
  a13=(a13*a0);
  a13=(a16*a13);
  a67=(a67-a13);
  a13=(a5*a67);
  a69=(a69+a13);
  a13=(a32*a56);
  a16=(a16*a13);
  a13=(a32*a71);
  a20=(a20*a13);
  a16=(a16+a20);
  a32=(a32*a65);
  a23=(a23*a32);
  a16=(a16+a23);
  a23=(a2*a16);
  a23=(a23+a34);
  a69=(a69+a23);
  a69=(a12*a69);
  a68=(a68+a69);
  a68=(-a68);
  if (res[1]!=0) res[1][37]=a68;
  a68=(a4*a62);
  a69=(a5*a16);
  a68=(a68+a69);
  a69=(a2*a67);
  a69=(a69+a27);
  a68=(a68-a69);
  a68=(a12*a68);
  a69=(a1*a9);
  a68=(a68-a69);
  if (res[1]!=0) res[1][38]=a68;
  a68=(a4*a67);
  a69=(a1*a16);
  a68=(a68-a69);
  a69=(a2*a62);
  a69=(a69+a17);
  a68=(a68+a69);
  a68=(a12*a68);
  a69=(a5*a9);
  a68=(a68-a69);
  if (res[1]!=0) res[1][39]=a68;
  a68=(a4*a16);
  a69=(a1*a67);
  a68=(a68+a69);
  a69=(a5*a62);
  a68=(a68-a69);
  a68=(a12*a68);
  a9=(a2*a9);
  a10=(a10-a9);
  a68=(a68+a10);
  if (res[1]!=0) res[1][40]=a68;
  a68=(a11*a71);
  a10=(a37*a65);
  a68=(a68-a10);
  a68=(-a68);
  if (res[1]!=0) res[1][41]=a68;
  a65=(a24*a65);
  a68=(a11*a56);
  a65=(a65-a68);
  a65=(-a65);
  if (res[1]!=0) res[1][42]=a65;
  a56=(a37*a56);
  a71=(a24*a71);
  a56=(a56-a71);
  a56=(-a56);
  if (res[1]!=0) res[1][43]=a56;
  a56=(a55*a62);
  a71=(a57*a67);
  a56=(a56+a71);
  a71=(a59*a16);
  a56=(a56+a71);
  a56=(a56/a58);
  a56=(a40*a56);
  a56=(-a56);
  if (res[1]!=0) res[1][44]=a56;
  a56=(a42*a67);
  a56=(a34*a56);
  a71=(a43*a16);
  a56=(a56+a71);
  a71=(a45*a16);
  a71=(a27*a71);
  a65=(a46*a67);
  a71=(a71+a65);
  a56=(a56-a71);
  a61=(a11*a61);
  a24=(a24*a0);
  a61=(a61-a24);
  a56=(a56-a61);
  a56=(a8*a56);
  a56=(a44*a56);
  a61=(a45*a16);
  a61=(a17*a61);
  a24=(a48*a62);
  a61=(a61+a24);
  a24=(a42*a62);
  a24=(a34*a24);
  a16=(a50*a16);
  a24=(a24+a16);
  a61=(a61-a24);
  a37=(a37*a60);
  a11=(a11*a66);
  a37=(a37+a11);
  a61=(a61+a37);
  a61=(a8*a61);
  a61=(a49*a61);
  a56=(a56+a61);
  a61=(a45*a62);
  a61=(a27*a61);
  a37=(a52*a67);
  a61=(a61+a37);
  a67=(a45*a67);
  a67=(a17*a67);
  a62=(a54*a62);
  a67=(a67+a62);
  a61=(a61-a67);
  a61=(a51*a61);
  a61=(a53*a61);
  a56=(a56+a61);
  a56=(a56/a47);
  if (res[1]!=0) res[1][45]=a56;
  a56=(-a7);
  if (res[1]!=0) res[1][46]=a56;
  if (res[1]!=0) res[1][47]=a6;
  a14=(a8*a14);
  a14=(a49*a14);
  a25=(a8*a25);
  a25=(a44*a25);
  a14=(a14-a25);
  a14=(a14/a47);
  if (res[1]!=0) res[1][48]=a14;
  if (res[1]!=0) res[1][49]=a7;
  a7=(-a3);
  if (res[1]!=0) res[1][50]=a7;
  a18=(a8*a18);
  a18=(a49*a18);
  a28=(a8*a28);
  a28=(a44*a28);
  a18=(a18-a28);
  a18=(a18/a47);
  if (res[1]!=0) res[1][51]=a18;
  a6=(-a6);
  if (res[1]!=0) res[1][52]=a6;
  if (res[1]!=0) res[1][53]=a3;
  a21=(a8*a21);
  a21=(a49*a21);
  a30=(a8*a30);
  a30=(a44*a30);
  a21=(a21-a30);
  a21=(a21/a47);
  if (res[1]!=0) res[1][54]=a21;
  a21=(a1*a15);
  a30=(a5*a26);
  a21=(a21+a30);
  a30=(a2*a33);
  a21=(a21+a30);
  a21=(a12*a21);
  a21=(-a21);
  if (res[1]!=0) res[1][55]=a21;
  a21=(a4*a15);
  a30=(a5*a33);
  a21=(a21+a30);
  a30=(a2*a26);
  a21=(a21-a30);
  a21=(a12*a21);
  if (res[1]!=0) res[1][56]=a21;
  a21=(a4*a26);
  a30=(a1*a33);
  a21=(a21-a30);
  a30=(a2*a15);
  a21=(a21+a30);
  a21=(a12*a21);
  if (res[1]!=0) res[1][57]=a21;
  a21=(a4*a33);
  a30=(a1*a26);
  a21=(a21+a30);
  a30=(a5*a15);
  a21=(a21-a30);
  a21=(a12*a21);
  if (res[1]!=0) res[1][58]=a21;
  a21=(a55*a15);
  a30=(a57*a26);
  a21=(a21+a30);
  a30=(a59*a33);
  a21=(a21+a30);
  a21=(a21/a58);
  a21=(a40*a21);
  a21=(-a21);
  if (res[1]!=0) res[1][59]=a21;
  a21=(a42*a26);
  a21=(a34*a21);
  a30=(a43*a33);
  a21=(a21+a30);
  a30=(a45*a33);
  a30=(a27*a30);
  a3=(a46*a26);
  a30=(a30+a3);
  a21=(a21-a30);
  a21=(a8*a21);
  a21=(a44*a21);
  a30=(a45*a33);
  a30=(a17*a30);
  a3=(a48*a15);
  a30=(a30+a3);
  a3=(a42*a15);
  a3=(a34*a3);
  a33=(a50*a33);
  a3=(a3+a33);
  a30=(a30-a3);
  a30=(a8*a30);
  a30=(a49*a30);
  a21=(a21+a30);
  a30=(a45*a15);
  a30=(a27*a30);
  a3=(a52*a26);
  a30=(a30+a3);
  a26=(a45*a26);
  a26=(a17*a26);
  a15=(a54*a15);
  a26=(a26+a15);
  a30=(a30-a26);
  a30=(a51*a30);
  a30=(a53*a30);
  a21=(a21+a30);
  a21=(a21/a47);
  if (res[1]!=0) res[1][60]=a21;
  a21=(a1*a19);
  a30=(a5*a29);
  a21=(a21+a30);
  a30=(a2*a35);
  a21=(a21+a30);
  a21=(a12*a21);
  a21=(-a21);
  if (res[1]!=0) res[1][61]=a21;
  a21=(a4*a19);
  a30=(a5*a35);
  a21=(a21+a30);
  a30=(a2*a29);
  a21=(a21-a30);
  a21=(a12*a21);
  if (res[1]!=0) res[1][62]=a21;
  a21=(a4*a29);
  a30=(a1*a35);
  a21=(a21-a30);
  a30=(a2*a19);
  a21=(a21+a30);
  a21=(a12*a21);
  if (res[1]!=0) res[1][63]=a21;
  a21=(a4*a35);
  a30=(a1*a29);
  a21=(a21+a30);
  a30=(a5*a19);
  a21=(a21-a30);
  a21=(a12*a21);
  if (res[1]!=0) res[1][64]=a21;
  a21=(a55*a19);
  a30=(a57*a29);
  a21=(a21+a30);
  a30=(a59*a35);
  a21=(a21+a30);
  a21=(a21/a58);
  a21=(a40*a21);
  a21=(-a21);
  if (res[1]!=0) res[1][65]=a21;
  a21=(a42*a29);
  a21=(a34*a21);
  a30=(a43*a35);
  a21=(a21+a30);
  a30=(a45*a35);
  a30=(a27*a30);
  a26=(a46*a29);
  a30=(a30+a26);
  a21=(a21-a30);
  a21=(a8*a21);
  a21=(a44*a21);
  a30=(a45*a35);
  a30=(a17*a30);
  a26=(a48*a19);
  a30=(a30+a26);
  a26=(a42*a19);
  a26=(a34*a26);
  a35=(a50*a35);
  a26=(a26+a35);
  a30=(a30-a26);
  a30=(a8*a30);
  a30=(a49*a30);
  a21=(a21+a30);
  a30=(a45*a19);
  a30=(a27*a30);
  a26=(a52*a29);
  a30=(a30+a26);
  a29=(a45*a29);
  a29=(a17*a29);
  a19=(a54*a19);
  a29=(a29+a19);
  a30=(a30-a29);
  a30=(a51*a30);
  a30=(a53*a30);
  a21=(a21+a30);
  a21=(a21/a47);
  if (res[1]!=0) res[1][66]=a21;
  a21=(a1*a22);
  a30=(a5*a31);
  a21=(a21+a30);
  a30=(a2*a36);
  a21=(a21+a30);
  a21=(a12*a21);
  a21=(-a21);
  if (res[1]!=0) res[1][67]=a21;
  a21=(a4*a22);
  a30=(a5*a36);
  a21=(a21+a30);
  a30=(a2*a31);
  a21=(a21-a30);
  a21=(a12*a21);
  if (res[1]!=0) res[1][68]=a21;
  a21=(a4*a31);
  a30=(a1*a36);
  a21=(a21-a30);
  a2=(a2*a22);
  a21=(a21+a2);
  a21=(a12*a21);
  if (res[1]!=0) res[1][69]=a21;
  a4=(a4*a36);
  a1=(a1*a31);
  a4=(a4+a1);
  a5=(a5*a22);
  a4=(a4-a5);
  a12=(a12*a4);
  if (res[1]!=0) res[1][70]=a12;
  a55=(a55*a22);
  a57=(a57*a31);
  a55=(a55+a57);
  a59=(a59*a36);
  a55=(a55+a59);
  a55=(a55/a58);
  a40=(a40*a55);
  a40=(-a40);
  if (res[1]!=0) res[1][71]=a40;
  a40=(a42*a31);
  a40=(a34*a40);
  a43=(a43*a36);
  a40=(a40+a43);
  a43=(a45*a36);
  a43=(a27*a43);
  a46=(a46*a31);
  a43=(a43+a46);
  a40=(a40-a43);
  a40=(a8*a40);
  a44=(a44*a40);
  a40=(a45*a36);
  a40=(a17*a40);
  a48=(a48*a22);
  a40=(a40+a48);
  a42=(a42*a22);
  a34=(a34*a42);
  a50=(a50*a36);
  a34=(a34+a50);
  a40=(a40-a34);
  a8=(a8*a40);
  a49=(a49*a8);
  a44=(a44+a49);
  a49=(a45*a22);
  a27=(a27*a49);
  a52=(a52*a31);
  a27=(a27+a52);
  a45=(a45*a31);
  a17=(a17*a45);
  a54=(a54*a22);
  a17=(a17+a54);
  a27=(a27-a17);
  a51=(a51*a27);
  a53=(a53*a51);
  a44=(a44+a53);
  a44=(a44/a47);
  if (res[1]!=0) res[1][72]=a44;
  a38=(a41*a38);
  a38=(-a38);
  if (res[1]!=0) res[1][73]=a38;
  a41=(a41*a39);
  a41=(-a41);
  if (res[1]!=0) res[1][74]=a41;
  a41=-1.;
  if (res[2]!=0) res[2][0]=a41;
  if (res[2]!=0) res[2][1]=a41;
  if (res[2]!=0) res[2][2]=a41;
  if (res[2]!=0) res[2][3]=a41;
  if (res[2]!=0) res[2][4]=a41;
  if (res[2]!=0) res[2][5]=a41;
  if (res[2]!=0) res[2][6]=a41;
  if (res[2]!=0) res[2][7]=a41;
  if (res[2]!=0) res[2][8]=a41;
  if (res[2]!=0) res[2][9]=a41;
  if (res[2]!=0) res[2][10]=a41;
  if (res[2]!=0) res[2][11]=a41;
  if (res[2]!=0) res[2][12]=a41;
  if (res[2]!=0) res[2][13]=a41;
  if (res[2]!=0) res[2][14]=a41;
  if (res[2]!=0) res[2][15]=a41;
  if (res[2]!=0) res[2][16]=a41;
  return 0;
}

CASADI_SYMBOL_EXPORT int tetherunit_integrator_impl_dae_fun_jac_x_xdot_u(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_release(int mem) {
}

CASADI_SYMBOL_EXPORT void tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_incref(void) {
}

CASADI_SYMBOL_EXPORT void tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_n_in(void) { return 5;}

CASADI_SYMBOL_EXPORT casadi_int tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_n_out(void) { return 4;}

CASADI_SYMBOL_EXPORT casadi_real tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_name_in(casadi_int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    case 3: return "i3";
    case 4: return "i4";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_name_out(casadi_int i){
  switch (i) {
    case 0: return "o0";
    case 1: return "o1";
    case 2: return "o2";
    case 3: return "o3";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    case 2: return casadi_s1;
    case 3: return casadi_s2;
    case 4: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s3;
    case 2: return casadi_s4;
    case 3: return casadi_s5;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int tetherunit_integrator_impl_dae_fun_jac_x_xdot_u_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 5;
  if (sz_res) *sz_res = 4;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
