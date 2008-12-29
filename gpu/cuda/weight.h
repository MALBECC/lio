#define MAX_ATOMS_PER_CUBE 20

// TODO: precomputar/precargar las cosas por atomo una sola vez para todos los puntos
// TODO: pasar esto a nucleii_count, etc

__global__ void gpu_compute_weights(uint points, uint* nuc_atom, float3* point_positions, uint nucleii_count, float* weights, uint* atom_of_points)
{
  uint3 pos = index(blockDim, blockIdx, threadIdx);
  uint point = pos.x;

  bool valid_thread = (point < points);
	if (!valid_thread) return;

  float3 point_position = point_positions[point];
  uint atom_of_point = atom_of_points[point];

  float P_total = 0.0f;
  float P_atom = 0.0f;

  #if GPU_WEIGHT_CUTOFF
  for (uint nuc_i = 0; nuc_i < nucleii_count; nuc_i++) {
    uint atom_i = nuc_atom[nuc_i];
  #else
  for (uint atom_i = 0; atom_i < gpu_atoms; atom_i++) {
  #endif
    float3 atom_i_pos = gpu_atom_positions[atom_i];
    float rm_atom_i = gpu_rm[atom_i];

    float P_curr = 1.0f;

  #if GPU_WEIGHT_CUTOFF
    for (uint nuc_j = 0; nuc_j < nucleii_count; nuc_j++) {
      if (nuc_i == nuc_j) continue;
      uint atom_j = nuc_atom[nuc_j];
  #else
    for (uint atom_j = 0; atom_j < gpu_atoms; atom_j++) {
  #endif
      if (atom_i == atom_j) continue;
      float3 atom_j_pos = gpu_atom_positions[atom_j];
      float rm_atom_j = gpu_rm[atom_j];

      float u = (::distance(point_position,atom_i_pos) - ::distance(point_position, atom_j_pos)) / ::distance(atom_i_pos, atom_j_pos);

			float x;
			x = rm_atom_i / rm_atom_j;
			x = (x - 1.0f) / (x + 1.0f);
			u += (x / (x * x - 1.0f)) * (1.0f - u * u);

      u = 1.5f * u - 0.5f * (u * u * u);
			u = 1.5f * u - 0.5f * (u * u * u);
			u = 1.5f * u - 0.5f * (u * u * u);
			u = 0.5f * (1.0f - u);

      P_curr *= u;
    }

    if (atom_i == atom_of_point) P_atom = P_curr;
    P_total += P_curr;
  }

  weights[point] = P_atom / P_total;
}
