// TODO: precomputar/precargar las cosas por atomo una sola vez para todos los puntos, pasar esto a nucleii_count, etc, coalescing

__global__ void gpu_compute_weights(uint points, float4* point_positions, float* weights)
{
  uint point = index(blockDim, blockIdx, threadIdx).x;

  bool valid_thread = (point < points);
  if (!valid_thread) return;

  float4 point_position4 = point_positions[point];
  float3 point_position = make_float3(point_position4.x, point_position4.y, point_position4.z);
  uint atom_of_point = (uint)floor(point_position4.w);

  float P_total = 0.0f;
  float P_atom = 0.0f;

  for (uint atom_i = 0; atom_i < gpu_atoms; atom_i++) {
    float3 atom_i_pos = gpu_atom_positions[atom_i];
    float rm_atom_i = gpu_rm[atom_i];

    float P_curr = 1.0f;

    for (uint atom_j = 0; atom_j < gpu_atoms; atom_j++) {
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
