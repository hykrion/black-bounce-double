#include "header/phi.h"

#include "header/ui.h"
#include "header/tortoise.h"

#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_deriv.h>

#define SYS_DIM 2
#define PHI_SYS_IND_0 15
#define PHI_SYS_IND_1 20

struct rsSqrtParams
{
  double A;
  double B;
  double rx;
};

struct rxParams
{
  double aa;
};

static double *m_phiRea;
static double *m_phiImg;
static double *m_T;
static double *m_R;
static double *m_v;
static double *m_deltaL;
static int m_nodes;

/**
  @brief  f(x) para calcular d(sqrt(A*B))/dx
*/
static double
funcSqrtABx(double r,
            void *params)
{
  (void)(params);
  double rS = ui_get_rS();
  double aa = ui_get_aa();
  double r_x = sqrt(pow(r, 2) + pow(aa, 2));

  /* sqrt(A*B) = A... porque A = B */
  return 1 - rS/r_x;
}

/**
  @brief  Función para calcular r''(x). En este caso aproximamos r'(x) por
          diferencias centrales.
*/
static double
funcRxx(double r,
        void *params)
{
  struct rxParams *par = (struct rxParams *)params;
  double aa = (par->aa);
  double h = 1e-3;

  double plus = sqrt(pow(r + h, 2) + pow(aa, 2));
  double minus = sqrt(pow(r - h, 2) + pow(aa, 2));

  return (plus - minus)/(2*h);
}

/**
  @brief  Función para calcular dr[x]/dx
*/
static double
funcRx(double r,
   void *params)
{
  struct rxParams *par = (struct rxParams *)params;
  double aa = (par->aa);

  return sqrt(r*r + aa*aa);
}

/*
  @brief  Ecuación de la barrera/pozo potencial. Es necesario haber calculado
          antes xy[]
*/
static double
phi_v(double r)
{
  /* Para las coordenadas tortuga */
  double a = ui_get_a();
  double hh = ui_get_h_forward();

  double aa = ui_get_aa();
  double l = ui_get_l();
  double rS = ui_get_rS();
  /* Para las derivadas */
  gsl_function Frx;
  gsl_function Frxx;
  gsl_function FSqrtABx;
  double h = 1e-8;
  double rx;
  double rxx;
  double sqrtABx;
  double absErr;

  /* Nos pasan un double... pero necesitamos un int */
  int i = (r - a) / hh;
  double xy = tortoise_get_xy_i(i);

  double r_x = sqrt(pow(xy, 2) + pow(aa, 2));
  double A = 1.0 - rS/r_x;
  struct rxParams rxPar = {aa};
  /* Calcular r'(x) */
  Frx.function = &funcRx;
  Frx.params = &rxPar;
  gsl_deriv_central(&Frx, xy, h, &rx, &absErr);

  /*
  Voy a calcular d(r'(x)*sqrt(A*B))/dx como
  r''(x)*sqrt(A*B) + r'(x)*d(sqrt(A*B))/dx
  */
  /* Calcular r''(x) */
  Frxx.function = &funcRxx;
  Frxx.params = &rxPar;
  gsl_deriv_central(&Frxx, xy, h, &rxx, &absErr);

  /*
  Calcular d(sqrt(A*B))/dx
  */
  FSqrtABx.function = &funcSqrtABx;
  FSqrtABx.params = NULL;
  gsl_deriv_central(&FSqrtABx, xy, h, &sqrtABx, &absErr);

  double ryy = A*(rxx*A + rx*sqrtABx);

  return (l*(1 + l)*A)/pow(r_x, 2) + ryy/r_x;
}

/**
  @brief  Jacobiana de phi
*/
static int
jacoPhi(double x,
        const double y[],
        double *dfdy,
        double dfdt[],
        void *params)
{
  (void)(y);
  (void)(params);
  double a = ui_get_a();
  double h = ui_get_h_forward();

  double aa = ui_get_aa();
  double rS = ui_get_rS();
  double l = ui_get_l();

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 1, 1);
  gsl_matrix *m = &dfdy_mat.matrix;
  /* Nos pasan un double... pero necesitamos un int */
  int i = (x - a) / h;
  double xy = tortoise_get_xy_i(i);

  double r = sqrt(pow(xy, 2) + pow(aa, 2));
  double A = 1.0 - rS/r;
  double x2aa2 = pow(xy, 2) + pow(aa, 2);
  double one = (l*rS*xy*(l + 1))/pow(x2aa2, 2.5);
  double two = (2*l*xy*(l + 1)*A)/pow(x2aa2, 2);
  double three = r*rS*(rS*pow(xy, 2)/pow(x2aa2, 2) + (A*(1.0 - pow(xy, 2)/x2aa2))/sqrt(x2aa2))/pow(x2aa2, 2);
  double four = xy*A*(rS*pow(xy, 2)/pow(x2aa2, 2) + A*(1.0 - (1.0 - pow(xy, 2)/x2aa2))/sqrt(x2aa2))/pow(x2aa2, 1.5);
  double five = A*(-4*rS*pow(xy, 3)/pow(x2aa2, 3) + rS*xy*(1.0 - pow(xy, 2)/pow(x2aa2, 2))/pow(x2aa2, 2) + 2*rS*xy/pow(x2aa2, 2) - xy*A*(1.0 - pow(xy, 2)/x2aa2)/pow(x2aa2, 1.5) + A*(2*pow(xy, 3)/pow(x2aa2, 2) - 2*xy/x2aa2)/sqrt(x2aa2))/sqrt(x2aa2);
  gsl_matrix_set (m, 0, 0, one - two + three - four + five);
  dfdt[0] = 0.0;

  return GSL_SUCCESS;
}

/**
  @brief  Sistema de ecuaciones

  @param  r: variable independiente
  @param  y: variables del sistema que permiten tener una ecuación de 1er orden
  @param  dydx: derivadas del sistema: 1a, 2a, 3a...
  @param  params: parámetros para las ecuaciones
*/
static int
funcPhi (double r,
         const double y[],
         double dydx[],
         void* params)
{
  struct TPhiParams *par = (struct TPhiParams *)params;
  double w = (par->w);

  double u = y[0];
  double v = y[1];

  dydx[0] = v;
  dydx[1] = (phi_v(r) - w*w)*u;

  return GSL_SUCCESS;
}

/**
  @brief  Realizar la integración de phi

  @param  IN, a: límite inferior del intervalo
  @param  IN, b: límite superior
  @param  IN, nodes: número de nodos
  @param  IN, ci: condición inicial
  @param  IN, intParams: parámetros necesarios para la integración
  @param  IN, OUT phi: arreglo con los valores de phi calculados

  @return Estado de la integración
*/
static int
phi_integration(double a,
                double b,
                int nodes,
                void *intParams,
                double *phi)
{
  int status = GSL_SUCCESS;

  struct TPhiParams *intPar = (struct TPhiParams*)intParams;
  double w = (intPar->w);
  double l = (intPar->l);
  double u = (intPar->u);
  double v = (intPar->v);

  double h = (b - a)/(nodes - 1);

  double x0 = a;
  double x1 = x0 + h;
  double epsAbs = 0;
  double epsRel = 1e-6;

  struct TPhiParams par = {w, l, 0, 0};

  gsl_odeiv2_system sys = {funcPhi, jacoPhi, SYS_DIM, &par};
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd; /* rkf45 rk8pd , rk4imp bsimp msadams msbdf*/
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new (&sys, T, h, epsAbs, epsRel);

  double y[2] = {u, v};
  int i;
  for (i = 0; i < nodes; i++)
  {
    status = gsl_odeiv2_driver_apply(d, &x0, x1, y);
    x0 = x1;
    x1 = x0 + h;

    if (status != GSL_SUCCESS)
    {
      printf ("error, return value = %d\n", status);
      printf("x = %f, y = %f\n", x0, y[0]);
      break;
    }
    else
      phi[i] = y[0];
  }
  gsl_odeiv2_driver_free(d);

  return status;
}

/* -----------------------------------------------------------------------------
  PUBLIC
----------------------------------------------------------------------------- */
void phi_init(void)
{
  m_nodes = ui_get_nodes();
  m_phiRea = calloc(m_nodes, sizeof(double));
  m_phiImg = calloc(m_nodes, sizeof(double));
  m_T = calloc(m_nodes, sizeof(double));
  m_R = calloc(m_nodes, sizeof(double));
  m_v = calloc(m_nodes, sizeof(double));
  m_deltaL = calloc(m_nodes, sizeof(double));
}

void phi_destroy(void)
{
  free(m_phiRea);
  free(m_phiImg);
  free(m_T);
  free(m_R);
  free(m_v);
  free(m_deltaL);
}

/**
  @brief  Calcula los coeficientes R y T entre un intervalo [a, b] para un rango
          de frecuencias angulares [wMin, wMax] para un momento angular concreto
          ,l, resolviendo directamente la EDO x'(y)

  @param  IN, l:  momento angular
  @param  IN, ci: condición inicial
*/
void
phi_wave_xy(double l,
            double ci)
{
  double a = ui_get_a();
  double b = ui_get_b();
  double wMin = ui_get_wMin();
  double hW = ui_get_hW();
  int nW = ui_get_nW();
  double rS = ui_get_rS();
  double aa = ui_get_aa();
  double w;
  struct tortoise_xyParams torParam = {rS, aa};

  /*
  Puntos alejados de la barrera donde resolveremos el sistema de ecuaciones
  para obtener los valores de R y T
  */
  tortoise_xy_integration(b, a, m_nodes, ci, &torParam);
  tortoise_reverse_xy();

  int i;
  for(i = 0; i < nW; i++)
  {
    /* Calcular 'q' de la ecuación y buscar índices de la barrera */
    w = wMin + hW*i;

    /* Integración hacia atrás */
    /* Parte real */
    double u = cos(w*b);
    double v = -w*sin(w*b);
    struct TPhiParams intParams = {w, l, u, v};
    phi_integration_rea(b, a, m_nodes, &intParams);

    /* Parte imaginaria */
    u = sin(w*b);
    v = w*cos(w*b);
    intParams.u = u;
    intParams.w = w;
    phi_integration_img(b, a, m_nodes, &intParams);

    /* Ordenar los vectores como si hubiéramos integrado hacia delante */
    phi_rea_fwd();
    phi_img_fwd();

    phi_calculate_RT(w, i);
  }
}

/**
  @brief  Calcula los coeficientes R y T entre un intervalo [a, b] para un rango
          de frecuencias angulares [wMin, wMax] para un momento angular concreto
          ,l, resolviendo la EDO y'(x), que nos servirá una vez resuelta para
          obtener la CI de la EDO x'(y)

  @param  IN, l:  momento angular
*/
void
phi_wave_yx_xy(double l)
{
  double a = ui_get_a();
  double b = ui_get_b();
  double wMin = ui_get_wMin();
  double hW = ui_get_hW();
  int nW = ui_get_nW();
  double w;

  /*
  Puntos alejados de la barrera donde resolveremos el sistema de ecuaciones
  para obtener los valores de R y T
  */
  int i;
  for(i = 0; i < nW; i++)
  {
    /* Calcular 'q' de la ecuación y buscar índices de la barrera */
    w = wMin + hW*i;

    /* Integración hacia atrás */

    /* Parte real */
    double u = cos(w*b);
    double v = -w*sin(w*b);
    struct TPhiParams intParams = {w, l, u, v};
    phi_integration_rea(b, a, m_nodes, &intParams);

    /* Parte imaginaria */
    u = sin(w*b);
    v = w*cos(w*b);
    intParams.u = u;
    intParams.v = v;
    phi_integration_img(b, a, m_nodes, &intParams);

    /* Ordenar los vectores como si hubiéramos integrado hacia delante */
    phi_rea_fwd();
    phi_img_fwd();

    phi_calculate_RT(w, i);
  }
}

/**
  @brief  Realizar la integración de phi

  @param  IN, a: límite inferior del intervalo
  @param  IN, b: límite superior
  @param  IN, nodes: número de nodos
  @param  IN, ci: condición inicial
  @param  IN, intParams: parámetros necesarios para la integración
  @param  IN, OUT phi: arreglo con los valores de phi calculados

  @return Estado de la integración
*/
int
phi_integration_rea(double a,
                    double b,
                    int nodes,
                    void *intParams)
{
  return phi_integration(a, b, nodes, intParams, m_phiRea);
}

/**
  @brief  Realizar la integración de phi

  @param  IN, a: límite inferior del intervalo
  @param  IN, b: límite superior
  @param  IN, nodes: número de nodos
  @param  IN, ci: condición inicial
  @param  IN, intParams: parámetros necesarios para la integración
  @param  IN, OUT phi: arreglo con los valores de phi calculados

  @return Estado de la integración
*/
int
phi_integration_img(double a,
                    double b,
                    int nodes,
                    void *intParams)
{
  return phi_integration(a, b, nodes, intParams, m_phiImg);
}

/**
  @brief  Ordenar phiRea como si se hubiera integrado hacia delante
*/
void
phi_rea_fwd(void)
{
  /* TODO Probar phi_reverse */
  double *tmp = calloc(m_nodes, sizeof(double));

  int i;
  for(i = 0; i < m_nodes; i++)
    tmp[i] = m_phiRea[m_nodes - 1 - i];
  for(i = 0; i < m_nodes; i++)
    m_phiRea[i] = tmp[i];

  free(tmp);
}

/**
  @brief  Ordenar phiRea como si se hubiera integrado hacia delante
*/
void
phi_img_fwd(void)
{
  /* TODO Probar phi_reverse */
  double* tmp = calloc(m_nodes, sizeof(double));

  int i;
  for(i = 0; i < m_nodes; i++)
    tmp[i] = m_phiImg[m_nodes - 1 - i];
  for(i = 0; i < m_nodes; i++)
    m_phiImg[i] = tmp[i];

  free(tmp);
}

/**
  @brief  Calcular los coeficientes R y T de la función de onda
*/
void
phi_calculate_RT(double w,
                 int i)
{
  double a = ui_get_a();
  double h = ui_get_h_forward();
  double B0 = m_phiRea[PHI_SYS_IND_0];
  double B1 = m_phiImg[PHI_SYS_IND_0];
  double B2 = m_phiRea[PHI_SYS_IND_1];
  double B3 = m_phiImg[PHI_SYS_IND_1];

  /* Dejamos las constantes en valor positivo */
  double k1 = -1*(a + PHI_SYS_IND_0*h);
  double k2 = -1*(a + PHI_SYS_IND_1*h);
  gsl_complex z1 = gsl_complex_polar(1, k1*w);
  gsl_complex z2 = gsl_complex_polar(1, k2*w);

  double e0r = GSL_REAL(z1);
  double e0i = GSL_IMAG(z1);
  double e1r = GSL_REAL(z2);
  double e1i = GSL_IMAG(z2);

  double Ar;
  double Ai;
  double Br;
  double Bi;
  gsl_complex AA;
  gsl_complex BB;

  double A[] = { e0r, e0i, e0r, -e0i,
                -e0i, e0r, e0i, e0r,
                 e1r, e1i, e1r, -e1i,
                -e1i, e1r, e1i, e1r };

  double B[] = { B0,
                 B1,
                 B2,
                 B3 };
  gsl_matrix_view m = gsl_matrix_view_array (A, 4, 4);
  gsl_vector_view b = gsl_vector_view_array (B, 4);
  gsl_vector *x = gsl_vector_alloc (4);
  int signum;
  gsl_permutation *permutation = gsl_permutation_alloc (4);
  gsl_linalg_LU_decomp (&m.matrix, permutation, &signum);
  gsl_linalg_LU_solve (&m.matrix, permutation, &b.vector, x);

  Ar = gsl_vector_get(x, 0);
  Ai = gsl_vector_get(x, 1);
  Br = gsl_vector_get(x, 2);
  Bi = gsl_vector_get(x, 3);
  AA = gsl_complex_rect(Ar, Ai);
  BB = gsl_complex_rect(Br, Bi);

  m_T[i] = 1.0 / gsl_complex_abs2(AA);
  m_R[i] = gsl_complex_abs2(BB) / gsl_complex_abs2(AA);

  gsl_permutation_free(permutation);
  gsl_vector_free(x);
}

/**
  @brief  Valores del potencial. Así podemos graficarlo
*/
void
phi_calculate_v(void)
{
  double a = ui_get_a();
  double h = ui_get_h_forward();

  double x;
  int i;
  for (i = 0; i < m_nodes; i++)
  {
    x = a + i*h;
    m_v[i] = phi_v(x);
  }
}

double*
phi_get_rea(void)
{
  return m_phiRea;
}

double*
phi_get_img(void)
{
  return m_phiImg;
}

double*
phi_get_v(void)
{
  return m_v;
}

double*
phi_get_R(void)
{
  return m_R;
}

double*
phi_get_T(void)
{
  return m_T;
}

/**
  @brief  Calcular delta_l para los diferentes momentos angulares, l.
*/
void
phi_calculate_delta_l(double l)
{
  double wMin = ui_get_wMin();
  double nW = ui_get_nW();
  double hW = ui_get_hW();
  double *R = phi_get_R();
  double w;

  int i;
  for(i = 0; i < nW; i++)
  {
    w = wMin + hW*i;
    m_deltaL[i] += (M_PI*(2*l + 1)*(1 - R[i]))/(w*w);
  }
}

double*
phi_get_deltaL(void)
{
  return m_deltaL;
}
