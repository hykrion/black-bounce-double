#ifndef TORTOISE_H_INCLUDED
#define TORTOISE_H_INCLUDED

struct tortoise_xyParams
{
  double rS;
  double aa;
};

void tortoise_init(void);
void tortoise_destroy(void);

void tortoise_calculate_yx_xy(void);
int
tortoise_xy_integration(double a,
                        double b,
                        int nodes,
                        double ci,
                        void *param);
int
tortoise_yx_integration(double a,
                        double b,
                        int nodes,
                        double ci,
                        void *param);
void tortoise_reverse_xy(void);
void tortoise_reverse_yx(void);
double* tortoise_get_xy(void);
double* tortoise_get_yx(void);
double tortoise_get_xy_i(int i);
double tortoise_get_yx_i(int i);
void tortoise_xy_analytical(void);
double* tortoise_get_xy_analytical(void);
void tortoise_yx_analytical(void);

#endif // TORTOISE_H_INCLUDED
