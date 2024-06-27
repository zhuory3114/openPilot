#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1941098188152223888) {
   out_1941098188152223888[0] = delta_x[0] + nom_x[0];
   out_1941098188152223888[1] = delta_x[1] + nom_x[1];
   out_1941098188152223888[2] = delta_x[2] + nom_x[2];
   out_1941098188152223888[3] = delta_x[3] + nom_x[3];
   out_1941098188152223888[4] = delta_x[4] + nom_x[4];
   out_1941098188152223888[5] = delta_x[5] + nom_x[5];
   out_1941098188152223888[6] = delta_x[6] + nom_x[6];
   out_1941098188152223888[7] = delta_x[7] + nom_x[7];
   out_1941098188152223888[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8037839476743691077) {
   out_8037839476743691077[0] = -nom_x[0] + true_x[0];
   out_8037839476743691077[1] = -nom_x[1] + true_x[1];
   out_8037839476743691077[2] = -nom_x[2] + true_x[2];
   out_8037839476743691077[3] = -nom_x[3] + true_x[3];
   out_8037839476743691077[4] = -nom_x[4] + true_x[4];
   out_8037839476743691077[5] = -nom_x[5] + true_x[5];
   out_8037839476743691077[6] = -nom_x[6] + true_x[6];
   out_8037839476743691077[7] = -nom_x[7] + true_x[7];
   out_8037839476743691077[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_923930073261382940) {
   out_923930073261382940[0] = 1.0;
   out_923930073261382940[1] = 0;
   out_923930073261382940[2] = 0;
   out_923930073261382940[3] = 0;
   out_923930073261382940[4] = 0;
   out_923930073261382940[5] = 0;
   out_923930073261382940[6] = 0;
   out_923930073261382940[7] = 0;
   out_923930073261382940[8] = 0;
   out_923930073261382940[9] = 0;
   out_923930073261382940[10] = 1.0;
   out_923930073261382940[11] = 0;
   out_923930073261382940[12] = 0;
   out_923930073261382940[13] = 0;
   out_923930073261382940[14] = 0;
   out_923930073261382940[15] = 0;
   out_923930073261382940[16] = 0;
   out_923930073261382940[17] = 0;
   out_923930073261382940[18] = 0;
   out_923930073261382940[19] = 0;
   out_923930073261382940[20] = 1.0;
   out_923930073261382940[21] = 0;
   out_923930073261382940[22] = 0;
   out_923930073261382940[23] = 0;
   out_923930073261382940[24] = 0;
   out_923930073261382940[25] = 0;
   out_923930073261382940[26] = 0;
   out_923930073261382940[27] = 0;
   out_923930073261382940[28] = 0;
   out_923930073261382940[29] = 0;
   out_923930073261382940[30] = 1.0;
   out_923930073261382940[31] = 0;
   out_923930073261382940[32] = 0;
   out_923930073261382940[33] = 0;
   out_923930073261382940[34] = 0;
   out_923930073261382940[35] = 0;
   out_923930073261382940[36] = 0;
   out_923930073261382940[37] = 0;
   out_923930073261382940[38] = 0;
   out_923930073261382940[39] = 0;
   out_923930073261382940[40] = 1.0;
   out_923930073261382940[41] = 0;
   out_923930073261382940[42] = 0;
   out_923930073261382940[43] = 0;
   out_923930073261382940[44] = 0;
   out_923930073261382940[45] = 0;
   out_923930073261382940[46] = 0;
   out_923930073261382940[47] = 0;
   out_923930073261382940[48] = 0;
   out_923930073261382940[49] = 0;
   out_923930073261382940[50] = 1.0;
   out_923930073261382940[51] = 0;
   out_923930073261382940[52] = 0;
   out_923930073261382940[53] = 0;
   out_923930073261382940[54] = 0;
   out_923930073261382940[55] = 0;
   out_923930073261382940[56] = 0;
   out_923930073261382940[57] = 0;
   out_923930073261382940[58] = 0;
   out_923930073261382940[59] = 0;
   out_923930073261382940[60] = 1.0;
   out_923930073261382940[61] = 0;
   out_923930073261382940[62] = 0;
   out_923930073261382940[63] = 0;
   out_923930073261382940[64] = 0;
   out_923930073261382940[65] = 0;
   out_923930073261382940[66] = 0;
   out_923930073261382940[67] = 0;
   out_923930073261382940[68] = 0;
   out_923930073261382940[69] = 0;
   out_923930073261382940[70] = 1.0;
   out_923930073261382940[71] = 0;
   out_923930073261382940[72] = 0;
   out_923930073261382940[73] = 0;
   out_923930073261382940[74] = 0;
   out_923930073261382940[75] = 0;
   out_923930073261382940[76] = 0;
   out_923930073261382940[77] = 0;
   out_923930073261382940[78] = 0;
   out_923930073261382940[79] = 0;
   out_923930073261382940[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_5423253775197792676) {
   out_5423253775197792676[0] = state[0];
   out_5423253775197792676[1] = state[1];
   out_5423253775197792676[2] = state[2];
   out_5423253775197792676[3] = state[3];
   out_5423253775197792676[4] = state[4];
   out_5423253775197792676[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5423253775197792676[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5423253775197792676[7] = state[7];
   out_5423253775197792676[8] = state[8];
}
void F_fun(double *state, double dt, double *out_4142146171863206567) {
   out_4142146171863206567[0] = 1;
   out_4142146171863206567[1] = 0;
   out_4142146171863206567[2] = 0;
   out_4142146171863206567[3] = 0;
   out_4142146171863206567[4] = 0;
   out_4142146171863206567[5] = 0;
   out_4142146171863206567[6] = 0;
   out_4142146171863206567[7] = 0;
   out_4142146171863206567[8] = 0;
   out_4142146171863206567[9] = 0;
   out_4142146171863206567[10] = 1;
   out_4142146171863206567[11] = 0;
   out_4142146171863206567[12] = 0;
   out_4142146171863206567[13] = 0;
   out_4142146171863206567[14] = 0;
   out_4142146171863206567[15] = 0;
   out_4142146171863206567[16] = 0;
   out_4142146171863206567[17] = 0;
   out_4142146171863206567[18] = 0;
   out_4142146171863206567[19] = 0;
   out_4142146171863206567[20] = 1;
   out_4142146171863206567[21] = 0;
   out_4142146171863206567[22] = 0;
   out_4142146171863206567[23] = 0;
   out_4142146171863206567[24] = 0;
   out_4142146171863206567[25] = 0;
   out_4142146171863206567[26] = 0;
   out_4142146171863206567[27] = 0;
   out_4142146171863206567[28] = 0;
   out_4142146171863206567[29] = 0;
   out_4142146171863206567[30] = 1;
   out_4142146171863206567[31] = 0;
   out_4142146171863206567[32] = 0;
   out_4142146171863206567[33] = 0;
   out_4142146171863206567[34] = 0;
   out_4142146171863206567[35] = 0;
   out_4142146171863206567[36] = 0;
   out_4142146171863206567[37] = 0;
   out_4142146171863206567[38] = 0;
   out_4142146171863206567[39] = 0;
   out_4142146171863206567[40] = 1;
   out_4142146171863206567[41] = 0;
   out_4142146171863206567[42] = 0;
   out_4142146171863206567[43] = 0;
   out_4142146171863206567[44] = 0;
   out_4142146171863206567[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_4142146171863206567[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_4142146171863206567[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4142146171863206567[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4142146171863206567[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_4142146171863206567[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_4142146171863206567[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_4142146171863206567[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_4142146171863206567[53] = -9.8000000000000007*dt;
   out_4142146171863206567[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_4142146171863206567[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_4142146171863206567[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4142146171863206567[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4142146171863206567[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_4142146171863206567[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_4142146171863206567[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_4142146171863206567[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4142146171863206567[62] = 0;
   out_4142146171863206567[63] = 0;
   out_4142146171863206567[64] = 0;
   out_4142146171863206567[65] = 0;
   out_4142146171863206567[66] = 0;
   out_4142146171863206567[67] = 0;
   out_4142146171863206567[68] = 0;
   out_4142146171863206567[69] = 0;
   out_4142146171863206567[70] = 1;
   out_4142146171863206567[71] = 0;
   out_4142146171863206567[72] = 0;
   out_4142146171863206567[73] = 0;
   out_4142146171863206567[74] = 0;
   out_4142146171863206567[75] = 0;
   out_4142146171863206567[76] = 0;
   out_4142146171863206567[77] = 0;
   out_4142146171863206567[78] = 0;
   out_4142146171863206567[79] = 0;
   out_4142146171863206567[80] = 1;
}
void h_25(double *state, double *unused, double *out_1889040035443611617) {
   out_1889040035443611617[0] = state[6];
}
void H_25(double *state, double *unused, double *out_187517270873638696) {
   out_187517270873638696[0] = 0;
   out_187517270873638696[1] = 0;
   out_187517270873638696[2] = 0;
   out_187517270873638696[3] = 0;
   out_187517270873638696[4] = 0;
   out_187517270873638696[5] = 0;
   out_187517270873638696[6] = 1;
   out_187517270873638696[7] = 0;
   out_187517270873638696[8] = 0;
}
void h_24(double *state, double *unused, double *out_2543769906956926137) {
   out_2543769906956926137[0] = state[4];
   out_2543769906956926137[1] = state[5];
}
void H_24(double *state, double *unused, double *out_282940211169700028) {
   out_282940211169700028[0] = 0;
   out_282940211169700028[1] = 0;
   out_282940211169700028[2] = 0;
   out_282940211169700028[3] = 0;
   out_282940211169700028[4] = 1;
   out_282940211169700028[5] = 0;
   out_282940211169700028[6] = 0;
   out_282940211169700028[7] = 0;
   out_282940211169700028[8] = 0;
   out_282940211169700028[9] = 0;
   out_282940211169700028[10] = 0;
   out_282940211169700028[11] = 0;
   out_282940211169700028[12] = 0;
   out_282940211169700028[13] = 0;
   out_282940211169700028[14] = 1;
   out_282940211169700028[15] = 0;
   out_282940211169700028[16] = 0;
   out_282940211169700028[17] = 0;
}
void h_30(double *state, double *unused, double *out_3637017401107411496) {
   out_3637017401107411496[0] = state[4];
}
void H_30(double *state, double *unused, double *out_58178323730398626) {
   out_58178323730398626[0] = 0;
   out_58178323730398626[1] = 0;
   out_58178323730398626[2] = 0;
   out_58178323730398626[3] = 0;
   out_58178323730398626[4] = 1;
   out_58178323730398626[5] = 0;
   out_58178323730398626[6] = 0;
   out_58178323730398626[7] = 0;
   out_58178323730398626[8] = 0;
}
void h_26(double *state, double *unused, double *out_9166996801521183285) {
   out_9166996801521183285[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3553986048000417528) {
   out_3553986048000417528[0] = 0;
   out_3553986048000417528[1] = 0;
   out_3553986048000417528[2] = 0;
   out_3553986048000417528[3] = 0;
   out_3553986048000417528[4] = 0;
   out_3553986048000417528[5] = 0;
   out_3553986048000417528[6] = 0;
   out_3553986048000417528[7] = 1;
   out_3553986048000417528[8] = 0;
}
void h_27(double *state, double *unused, double *out_5367103681649840092) {
   out_5367103681649840092[0] = state[3];
}
void H_27(double *state, double *unused, double *out_2116584988070026285) {
   out_2116584988070026285[0] = 0;
   out_2116584988070026285[1] = 0;
   out_2116584988070026285[2] = 0;
   out_2116584988070026285[3] = 1;
   out_2116584988070026285[4] = 0;
   out_2116584988070026285[5] = 0;
   out_2116584988070026285[6] = 0;
   out_2116584988070026285[7] = 0;
   out_2116584988070026285[8] = 0;
}
void h_29(double *state, double *unused, double *out_8826802833295998367) {
   out_8826802833295998367[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3829947714939577318) {
   out_3829947714939577318[0] = 0;
   out_3829947714939577318[1] = 1;
   out_3829947714939577318[2] = 0;
   out_3829947714939577318[3] = 0;
   out_3829947714939577318[4] = 0;
   out_3829947714939577318[5] = 0;
   out_3829947714939577318[6] = 0;
   out_3829947714939577318[7] = 0;
   out_3829947714939577318[8] = 0;
}
void h_28(double *state, double *unused, double *out_365430413988994447) {
   out_365430413988994447[0] = state[0];
}
void H_28(double *state, double *unused, double *out_1866317443374251067) {
   out_1866317443374251067[0] = 1;
   out_1866317443374251067[1] = 0;
   out_1866317443374251067[2] = 0;
   out_1866317443374251067[3] = 0;
   out_1866317443374251067[4] = 0;
   out_1866317443374251067[5] = 0;
   out_1866317443374251067[6] = 0;
   out_1866317443374251067[7] = 0;
   out_1866317443374251067[8] = 0;
}
void h_31(double *state, double *unused, double *out_2605511692490013296) {
   out_2605511692490013296[0] = state[8];
}
void H_31(double *state, double *unused, double *out_218163232750599124) {
   out_218163232750599124[0] = 0;
   out_218163232750599124[1] = 0;
   out_218163232750599124[2] = 0;
   out_218163232750599124[3] = 0;
   out_218163232750599124[4] = 0;
   out_218163232750599124[5] = 0;
   out_218163232750599124[6] = 0;
   out_218163232750599124[7] = 0;
   out_218163232750599124[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_1941098188152223888) {
  err_fun(nom_x, delta_x, out_1941098188152223888);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8037839476743691077) {
  inv_err_fun(nom_x, true_x, out_8037839476743691077);
}
void car_H_mod_fun(double *state, double *out_923930073261382940) {
  H_mod_fun(state, out_923930073261382940);
}
void car_f_fun(double *state, double dt, double *out_5423253775197792676) {
  f_fun(state,  dt, out_5423253775197792676);
}
void car_F_fun(double *state, double dt, double *out_4142146171863206567) {
  F_fun(state,  dt, out_4142146171863206567);
}
void car_h_25(double *state, double *unused, double *out_1889040035443611617) {
  h_25(state, unused, out_1889040035443611617);
}
void car_H_25(double *state, double *unused, double *out_187517270873638696) {
  H_25(state, unused, out_187517270873638696);
}
void car_h_24(double *state, double *unused, double *out_2543769906956926137) {
  h_24(state, unused, out_2543769906956926137);
}
void car_H_24(double *state, double *unused, double *out_282940211169700028) {
  H_24(state, unused, out_282940211169700028);
}
void car_h_30(double *state, double *unused, double *out_3637017401107411496) {
  h_30(state, unused, out_3637017401107411496);
}
void car_H_30(double *state, double *unused, double *out_58178323730398626) {
  H_30(state, unused, out_58178323730398626);
}
void car_h_26(double *state, double *unused, double *out_9166996801521183285) {
  h_26(state, unused, out_9166996801521183285);
}
void car_H_26(double *state, double *unused, double *out_3553986048000417528) {
  H_26(state, unused, out_3553986048000417528);
}
void car_h_27(double *state, double *unused, double *out_5367103681649840092) {
  h_27(state, unused, out_5367103681649840092);
}
void car_H_27(double *state, double *unused, double *out_2116584988070026285) {
  H_27(state, unused, out_2116584988070026285);
}
void car_h_29(double *state, double *unused, double *out_8826802833295998367) {
  h_29(state, unused, out_8826802833295998367);
}
void car_H_29(double *state, double *unused, double *out_3829947714939577318) {
  H_29(state, unused, out_3829947714939577318);
}
void car_h_28(double *state, double *unused, double *out_365430413988994447) {
  h_28(state, unused, out_365430413988994447);
}
void car_H_28(double *state, double *unused, double *out_1866317443374251067) {
  H_28(state, unused, out_1866317443374251067);
}
void car_h_31(double *state, double *unused, double *out_2605511692490013296) {
  h_31(state, unused, out_2605511692490013296);
}
void car_H_31(double *state, double *unused, double *out_218163232750599124) {
  H_31(state, unused, out_218163232750599124);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
