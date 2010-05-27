#ifndef _CPU_POT_H
#define	_CPU_POT_H

void cpu_pot(float dens, double& ex, double& ec, double& y2a);
void cpu_potg(float dens, const double3& grad, const double3& hess1, const double3& hess2, double& ex, double& ec, double& y2a);

#endif	/* _CPU_POT_H */

