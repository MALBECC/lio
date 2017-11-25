#ifndef _EXCHNUM_H
#define _EXCHNUM_H

template <bool compute_energy, bool do_forces>
void cpu_compute_density_forces(float* energy, float* point_weights,
                                uint points, float* rdm, float* rmm_output,
                                float* function_values, float4* gradient_values,
                                float4* forces, uint* nuc, uint nucleii_count,
                                uint m, Timer& t, Timer& trmm);

#endif /* _EXCHNUM_H */
