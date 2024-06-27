#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_1941098188152223888);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8037839476743691077);
void car_H_mod_fun(double *state, double *out_923930073261382940);
void car_f_fun(double *state, double dt, double *out_5423253775197792676);
void car_F_fun(double *state, double dt, double *out_4142146171863206567);
void car_h_25(double *state, double *unused, double *out_1889040035443611617);
void car_H_25(double *state, double *unused, double *out_187517270873638696);
void car_h_24(double *state, double *unused, double *out_2543769906956926137);
void car_H_24(double *state, double *unused, double *out_282940211169700028);
void car_h_30(double *state, double *unused, double *out_3637017401107411496);
void car_H_30(double *state, double *unused, double *out_58178323730398626);
void car_h_26(double *state, double *unused, double *out_9166996801521183285);
void car_H_26(double *state, double *unused, double *out_3553986048000417528);
void car_h_27(double *state, double *unused, double *out_5367103681649840092);
void car_H_27(double *state, double *unused, double *out_2116584988070026285);
void car_h_29(double *state, double *unused, double *out_8826802833295998367);
void car_H_29(double *state, double *unused, double *out_3829947714939577318);
void car_h_28(double *state, double *unused, double *out_365430413988994447);
void car_H_28(double *state, double *unused, double *out_1866317443374251067);
void car_h_31(double *state, double *unused, double *out_2605511692490013296);
void car_H_31(double *state, double *unused, double *out_218163232750599124);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}