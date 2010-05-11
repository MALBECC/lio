#ifndef _CPU_POT_H
#define	_CPU_POT_H

void cpu_pot(float dens, float& ex, float& ec, float& y2a);
void cpu_potg(float dens, const double3& grad, const double3& hess1, const double3& hess2, float& ex, float& ec, float& y2a);

#endif	/* _CPU_POT_H */

