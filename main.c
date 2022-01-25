#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "header/ui.h"
#include "header/phi.h"
#include "header/tortoise.h"
#include "header/txt_data.h"

/**
  @brief  Calcular phi resolviendo la EDO y'(x) para con un xFijo, obtener
          y(xFijo) = yFijo, y usando yFijo resolver la edo x'(y) con
          x(yFijo) = xFijo como CI de la EDO.
*/
void
calculate_phi_yx_xy(void)
{
  double l = ui_get_l();
  tortoise_calculate_yx_xy();
  /* Realizar el mismo cálculo para una serie de momentos angulares */
  int nL = ui_get_nL();
  int i;
  for(i = 0; i < nL; i++)
  {
    phi_wave_yx_xy(l);
    phi_calculate_delta_l(l);
    l += 1.0;
    /*
    Como no es el momento angular que hay en la GUI... tenemos que guardar
    la modificación
    */
    ui_set_l(l);
  }
}
/* -----------------------------------------------------------------------------
----------------------------------------------------------------------------- */
int
main()
{
  /* Leer los parámetros creados por el GUI */
  ui_init_only();
  int nodes = ui_get_nodes();

  /*
  Calcular phi resolviendo y'(x) para una vez resuelto obtener una CI para
  resolver x(y')
  */
  phi_init();
  tortoise_init();
  calculate_phi_yx_xy();

  /* Obtener los valores del potencial para dibujarlos */
  phi_calculate_v();

  txt_data_init_only(nodes);
  /* phi Re, Im y potencial */
  txt_data_wave();

  /* Coordenadas tortuga */
  txt_data_turtle();

  /* Valores de los coeficientes R, T y R + T; dependiendo de w */
  txt_data_rt();

  /* Valores delta_l */
  txt_data_delta();

  tortoise_destroy();
  phi_destroy();

  return 0;
}
