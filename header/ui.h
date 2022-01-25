#ifndef UI_H_INCLUDED
#define UI_H_INCLUDED

void ui_init_only(void);

double ui_get_a(void);
double ui_get_b(void);
double ui_get_h_forward(void);
double ui_get_h_backguard(void);
int ui_get_nodes(void);
int ui_get_n(void);
double ui_get_wMin(void);
double ui_get_wMax(void);
double ui_get_hW(void);
int ui_get_nW(void);
double ui_get_rS(void);
double ui_get_l(void);
void ui_set_l(double l);
double ui_get_nL(void);
double ui_get_aa(void);

#endif // UI_H_INCLUDED
