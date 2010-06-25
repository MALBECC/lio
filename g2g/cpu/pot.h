#ifndef _CPU_POT_H
#define	_CPU_POT_H

void cpu_pot(float dens, float& ex, float& ec, float& y2a);
void cpu_potg(float dens, const float3& grad, const float3& hess1, const float3& hess2, float& ex, float& ec, float& y2a);

#endif	/* _CPU_POT_H */

