#ifndef NUMERIC_INTEGRATION_H_INCLUDED
#define NUMERIC_INTEGRATION_H_INCLUDED


struct schwarzschildParams
{
  double w;
};

int
backguardIntegrationRKSchwarzschild(double x,
                                  double u,
                                  double v,
                                  double uu[],
                                  double h,
                                  double w);

#endif
