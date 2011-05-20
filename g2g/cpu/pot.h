#ifndef _CPU_POT_H
#define	_CPU_POT_H

void cpu_pot(real dens, real& ex, real& ec, real& y2a);
void cpu_potg(real dens, const real3& grad, const real3& hess1, const real3& hess2, real& ex, real& ec, real& y2a);

#endif	/* _CPU_POT_H */

