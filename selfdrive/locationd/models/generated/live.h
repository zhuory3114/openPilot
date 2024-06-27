#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_4157293951668813942);
void live_err_fun(double *nom_x, double *delta_x, double *out_4215917468253811954);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_2607967756817656059);
void live_H_mod_fun(double *state, double *out_763305028412363756);
void live_f_fun(double *state, double dt, double *out_1666349005085250907);
void live_F_fun(double *state, double dt, double *out_4142580856902382037);
void live_h_4(double *state, double *unused, double *out_3045987888290929190);
void live_H_4(double *state, double *unused, double *out_2097197187971950069);
void live_h_9(double *state, double *unused, double *out_2344445991804938187);
void live_H_9(double *state, double *unused, double *out_5190021747292497401);
void live_h_10(double *state, double *unused, double *out_3346531318811293775);
void live_H_10(double *state, double *unused, double *out_9201430287860347608);
void live_h_12(double *state, double *unused, double *out_2940339153109469531);
void live_H_12(double *state, double *unused, double *out_2922259220060011726);
void live_h_35(double *state, double *unused, double *out_2031125974194172514);
void live_H_35(double *state, double *unused, double *out_5667822252385025435);
void live_h_32(double *state, double *unused, double *out_8181692891130292531);
void live_H_32(double *state, double *unused, double *out_1400629654341397776);
void live_h_13(double *state, double *unused, double *out_7923644008632447371);
void live_H_13(double *state, double *unused, double *out_2582932054664047527);
void live_h_14(double *state, double *unused, double *out_2344445991804938187);
void live_H_14(double *state, double *unused, double *out_5190021747292497401);
void live_h_33(double *state, double *unused, double *out_1334516723763015635);
void live_H_33(double *state, double *unused, double *out_8818379257023883039);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}