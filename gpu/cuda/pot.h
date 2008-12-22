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

/* pot_kernel constants */
#define POT_ALPHA 		-0.738558766382022447
#define POT_GL 				0.620350490899400087
#define POT_VOSKO_A1 	0.03109205
#define POT_VOSKO_B1 	3.72744
#define POT_VOSKO_C1 	12.9352
#define POT_VOSKO_X0 	-0.10498
#define POT_VOSKO_Q 	6.15199066246304849
#define POT_VOSKO_A16 0.005182008333
#define POT_VOSKO_A2 	0.015546025
#define POT_VOSKO_B2 	7.06042
#define POT_VOSKO_C2 	18.0578
#define POT_VOSKO_X02	-0.32500
#define POT_VOSKO_Q2 	4.7309269
#define POT_VOSKO_A26 0.0025910042
	
template<bool compute_exc, bool compute_y2a> __device__ void gpu_pot(float dens, float& ex, float& ec, float& y2a)
{
	// data X alpha

	if (dens == 0) {
		if (compute_exc) { ex = 0.0f; ec = 0.0f; }
		if (compute_y2a) y2a = 0.0f;
		return;
	}
	
	float y = powf(dens, 0.333333333333333333f);
	float e0 = POT_ALPHA * y;	
	float v0 = -0.984745021842697f * y; // 4/3 * pot_alpha * y

	if (compute_exc) ex = e0;
	
	switch(gpu_Iexch) {
		case 1:
		{		
			if (compute_exc) ec = 0;
			if (compute_y2a) y2a = v0;
		}
		break;
		case 2:
		{
			float rs = POT_GL / y;
			float x1 = rs / 11.4f;
			float vc;
			
			if (x1 > 1.0f) {
				ec = -0.0333f * (0.5f * x1 - 0.33333333333333f);
				if (compute_y2a) vc = 0.0111f * x1 * 0.5f;
			}
			else {
				float t1 = (1.0f + x1 * x1 * x1);
				float t2 = logf(1.0f + 1.0f / x1);
				float t3 = x1 * x1;
        ec = -0.0333f * (t1 * t2 - t3 + 0.5f * x1 - 0.33333333333333f);
        if (compute_y2a) vc = 0.0111f * x1 * (3.0f * t3 * t2 - t1 / (x1 * (x1 + 1.0f)) - 2.0f * x1 + 0.5f);
			}
			if (compute_y2a) y2a = v0 + ec + vc;
		}
		break;
		case 3:
		{
			float rs = POT_GL / y;
			float x1 = sqrtf(rs);
			float Xx = rs + POT_VOSKO_B1 * x1 + POT_VOSKO_C1;
			float Xxo = POT_VOSKO_X0 * POT_VOSKO_X0 + POT_VOSKO_B1 * POT_VOSKO_X0 + POT_VOSKO_C1;
			float t1 = 2 * x1 + POT_VOSKO_B1;
			float t2 = logf(Xx);
			float t3 = atanf(POT_VOSKO_Q/t1);
			float t4 = POT_VOSKO_B1 * POT_VOSKO_X0 / Xxo;
      float t5 = (POT_VOSKO_B1 * x1 + 2.0f * POT_VOSKO_C1) / x1;
			float t6 = POT_VOSKO_X0 / Xxo;			
			
      ec = POT_VOSKO_A1 * (2 * logf(x1) - t2 + 2 * POT_VOSKO_B1 / POT_VOSKO_Q * t3 - t4 *
								 (2 * logf(x1 - POT_VOSKO_X0) - t2 + 2 * (POT_VOSKO_B1 + 2 * POT_VOSKO_X0) / POT_VOSKO_Q * t3));
			
			float vc;
      if (compute_y2a) {
				vc = ec - POT_VOSKO_A16 * x1 * (t5 / Xx - 4.0f * POT_VOSKO_B1 / (t1 * t1 + POT_VOSKO_Q * POT_VOSKO_Q) * (1.0f - t6 * (POT_VOSKO_B1 - 2.0f * POT_VOSKO_X0)) -
																						t4 * (2.0f / (x1 - POT_VOSKO_X0) - t1 / Xx));
				y2a = v0 + vc;
			}
		}
		break;		
	}
}
