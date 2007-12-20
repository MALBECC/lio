/*
 * function for evaluating exchange correlation density and
 * potential, depends of Iexch, index that says which potential
 * will be used
 * - 1 X alpha
 * - 2 Gunnarson-Lundqvist
 * - 3 Vosko, Wilk and Nusair
 * Exchange part is in all cases the same, for the time being LD
 * Self interaction corrections are used in correlation part
 */
	
__device__ void pot_kernel(float dens, float& ex, float& ec, float& y2a, uint Iexch, uint big_index)
{
	// data X alpha

	if (dens == 0) { ex = 0.0f; ec = 0.0f; y2a = 0.0f; return; }
	
	float y = powf(dens, 0.333333333333333333f);
	float e0 = pot_alpha * y;	
	float v0 = -0.984745021842697f * y; // 4/3 * pot_alpha * y
	
	switch(Iexch) {
		case 1:
		{
			ex = e0;
			ec = 0;
			y2a = v0;
		}
		break;
		case 2:
		{
			ex = e0;
			float rs = pot_gl / y;
			float x1 = rs / 11.4f;
			float vc;
			
			if (x1 > 1.0f) {
				ec = -0.0333f * (0.5f * x1 - 0.33333333333333f);
				vc = 0.0111f * x1 * 0.5f;
			}
			else {
				float t1 = (1.0f + x1 * x1 * x1);
				float t2 = logf(1.0f + 1.0f / x1);
				float t3 = x1 * x1;
        ec = -0.0333f * (t1 * t2 - t3 + 0.5f * x1 - 0.33333333333333f);
        vc = 0.0111f * x1 * (3.0f * t3 * t2 - t1 / (x1 * (x1 + 1.0f)) - 2.0f * x1 + 0.5f);
			}
			y2a = v0 + ec + vc;
		}
		break;
		case 3:
		{
			ex = e0;
			float rs = pot_gl / y;
			float x1 = sqrtf(rs);
			float Xx = rs + pot_vosko_b1 * x1 + pot_vosko_c1;
			float Xxo = pot_vosko_x0 * pot_vosko_x0 + pot_vosko_b1 * pot_vosko_x0 + pot_vosko_c1;
			float t1 = 2 * x1 + pot_vosko_b1;
			float t2 = logf(Xx);
			float t3 = atanf(pot_vosko_q/t1);
			float t4 = pot_vosko_b1 * pot_vosko_x0 / Xxo;
      float t5 = (pot_vosko_b1 * x1 + 2.0f * pot_vosko_c1) / x1;
			float t6 = pot_vosko_x0 / Xxo;			
			
      ec = pot_vosko_a1 * (2 * logf(x1) - t2 + 2 * pot_vosko_b1 / pot_vosko_q * t3 - t4 *
								 (2 * logf(x1 - pot_vosko_x0) - t2 + 2 * (pot_vosko_b1 + 2 * pot_vosko_x0) / pot_vosko_q * t3));
			
      float vc = ec - pot_vosko_a16 * x1 * (t5 / Xx - 4.0f * pot_vosko_b1 / (t1 * t1 + pot_vosko_q * pot_vosko_q) * (1.0f - t6 * (pot_vosko_b1 - 2.0f * pot_vosko_x0)) -
																						t4 * (2.0f / (x1 - pot_vosko_x0) - t1 / Xx));
			
			y2a = v0 + vc;
		}
		break;		
	}
}
