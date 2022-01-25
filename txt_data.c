#include "header/txt_data.h"

#include "header/ui.h"
#include "header/phi.h"
#include "header/tortoise.h"

#include <stdio.h>

static int m_nodes;

/* -----------------------------------------------------------------------------
  PUBLIC
----------------------------------------------------------------------------- */
void
txt_data_init_only(int nodes)
{
  m_nodes = nodes;
}

void
txt_data_wave(void)
{
  FILE *fp;
  fp = fopen("wave.txt", "w");

  double a = ui_get_a();
  double h = ui_get_h_forward();
  double x = a;
  double* phiRea = phi_get_rea();
  double* phiImg = phi_get_img();
  double* phiV = phi_get_v();

  int i;
  for (i = 0; i < m_nodes; i++)
  {
    fprintf(fp, "%.32f\t%.32f\t%.32f\t%.32f\n", x, phiRea[i], phiImg[i], phiV[i]);
    x += h;
  }
  fclose(fp);
}

/**
  @brief  Comparar el valor calculado con el analítico
*/
void
txt_data_turtle()
{
  FILE *fp;
  fp = fopen("turtle.txt", "w");

  double a = ui_get_a();
  double h = ui_get_h_forward();
  double x = a;
  /* Valor calculado */
  double *xy = tortoise_get_xy();
  /* Valor analítico. Lo invierto en gnuplot para que sea xy */
  double *yx = tortoise_get_yx();

  int i;
  for (i = 0; i < m_nodes; i++)
  {
    fprintf(fp, "%.32f\t%.32f\t%.32f\n", x, xy[i], yx[i]);
    x += h;
  }
  fclose(fp);
}

void
txt_data_rt()
{
  FILE *fp;
  fp = fopen("coefficients.txt", "w");

  int nW = ui_get_nW();
  double wMin = ui_get_wMin();
  double hW = ui_get_hW();
  double *R = phi_get_R();
  double *T = phi_get_T();
  double w;

  int i;
  for (i = 0; i < nW; i++)
  {
    w = wMin + hW*i;
    fprintf(fp, "%.32f\t%.32f\t%.32f\t%.32f\n", w, R[i], T[i], R[i] + T[i]);
  }
  fclose(fp);
}

void
txt_data_delta()
{

  FILE *fp;
  fp = fopen("sigma-l.txt", "w");

  double wMin = ui_get_wMin();
  int nW = ui_get_nW();
  double hW = ui_get_hW();
  double *deltaL = phi_get_deltaL();
  double w;

  int i;
  for (i = 0; i < nW; i++)
  {
    w = wMin + hW*i;
    fprintf(fp, "%.32f\t%.32f\t\n", w, deltaL[i]);
  }
  fclose(fp);
}
