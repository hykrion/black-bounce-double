#include "header/ui.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Número de parámetros + 1 */
#define UI_PARAM_SIZE 10
/* Longitud máxima + 1 del valor de los parámetros */
#define UI_PARAM_LEN 16

/*
NOTE  El número y orden de los parámetros ha de coincidir con los declarados
      en 'cmd_visualizar' del fichero 'gui.tcl'.
*/
#define UI_PARAM_A 0
#define UI_PARAM_B 1
#define UI_PARAM_L 2
#define UI_PARAM_RS 3
#define UI_PARAM_WMIN 4
#define UI_PARAM_WMAX 5
#define UI_PARAM_NW 6
#define UI_PARAM_NL 7
#define UI_PARAM_AA 8
#define UI_PARAM_NODES 9

static double m_uiParams[UI_PARAM_SIZE];

/**
  @brief  Leer los diferentes parámetros que modificaresmos desde la GUI
          Por simplicidad estos parámetros están en modo texto en un fichero,
          parameters.txt, que tendrá un parámetro por cada línea.
  */
void
ui_init_only(void)
{
  char str[UI_PARAM_SIZE][UI_PARAM_LEN];

  int i;
  for(i = 0; i < UI_PARAM_SIZE; i++)
    memset(str[i], '\0', UI_PARAM_LEN);

  /* Leer los parámetros de fichero */
  FILE *fp;
  fp = fopen("parameters.txt", "r");

  i = 0;
  while(!feof(fp))
  {
    fgets(str[i], UI_PARAM_LEN - 1, fp);
    i++;
  }

  fclose(fp);

  i = 0;
  while(i < UI_PARAM_SIZE)
  {
    m_uiParams[i] = atof(str[i]);
    i++;
  }
}

double
ui_get_a(void)
{
  return m_uiParams[UI_PARAM_A];
}

double
ui_get_b(void)
{
  return m_uiParams[UI_PARAM_B];
}

double
ui_get_h_forward(void)
{
  double a = ui_get_a();
  double b = ui_get_b();
  double n = ui_get_n();

  return (b - a) / n;
}

double
ui_get_h_backguard(void)
{
  double a = ui_get_a();
  double b = ui_get_b();
  double n = ui_get_n();

  return (a - b) / n;
}

int
ui_get_nodes(void)
{
  return m_uiParams[UI_PARAM_NODES];
}

int
ui_get_n(void)
{
  return m_uiParams[UI_PARAM_NODES] - 1;
}

double
ui_get_wMin(void)
{
  return m_uiParams[UI_PARAM_WMIN];
}

double
ui_get_wMax(void)
{
  return m_uiParams[UI_PARAM_WMAX];
}

double
ui_get_hW(void)
{
  double wm = ui_get_wMin();
  double wM = ui_get_wMax();
  double nW = ui_get_nW();

  return (wM - wm) / nW;
}

int
ui_get_nW(void)
{
  return m_uiParams[UI_PARAM_NW];
}

double
ui_get_rS(void)
{
  return m_uiParams[UI_PARAM_RS];
}

double
ui_get_l(void)
{
  return m_uiParams[UI_PARAM_L];
}

double
ui_get_nL(void)
{
  return m_uiParams[UI_PARAM_NL];
}

void
ui_set_l(double l)
{
  m_uiParams[UI_PARAM_L] = l;
}

double
ui_get_aa(void)
{
  return m_uiParams[UI_PARAM_AA];
}
