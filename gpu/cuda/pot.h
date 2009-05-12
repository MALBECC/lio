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
#define POT_ALPHA 		-0.738558766382022447f // -(3/PI)^(1/3)
#define POT_GL 				0.620350490899400087f

#define POT_VOSKO_A1 	0.03109205f
#define POT_VOSKO_B1 	3.72744f
#define POT_VOSKO_C1 	12.9352f
#define POT_VOSKO_X0 	-0.10498f

#define POT_VOSKO_Q 	6.15199066246304849f
#define POT_VOSKO_A16 0.005182008333f
#define POT_VOSKO_A2 	0.015546025f
#define POT_VOSKO_B2 	7.06042f
#define POT_VOSKO_C2 	18.0578f
#define POT_VOSKO_X02 -0.32500f
#define POT_VOSKO_Q2 	4.7309269f
#define POT_VOSKO_A26 0.0025910042f

#define POT_XX0 12.5549141492f // POT_VOSKO_X0 * POT_VOSKO_X0 + POT_VOSKO_B1 * POT_VOSKO_X0 + POT_VOSKO_C1
#define POT_T6 -0.00836166609762834f // POT_VOSKO_X0 / POT_XX0
#define POT_T4 -0.0311676086789438f // POT_VOSKO_B1 * POT_VOSKO_X0 / POT_XX0
#define POT_VOSKO_2C1 25.8704f // 2 * POT_VOSKO_C1

#define POT_VOSKO_2B1Q 1.21178337371132f // 2 * POT_VOSKO_B1 / POT_VOSKO_Q
#define POT_VOSKO_B2X0Q 1.14352579286644f // 2 * (POT_VOSKO_B1 + 2 * POT_VOSKO_X0) / POT_VOSKO_Q
#define POT_VOSKO_QSQ 37.8469891110325f // POT_VOSKO_Q * POT_VOSKO_Q
#define POT_VOSKO_4B1B1X0 15.4006373696499f // 4.0 * POT_VOSKO_B1 * (1.0f - t6 * (POT_VOSKO_B1 - 2.0f * POT_VOSKO_X0))
	
template<bool compute_exc, bool compute_y2a> __device__ void gpu_pot(float dens, float& exc_corr, float& y2a)
{
	// data X alpha

	if (dens == 0) {
		if (compute_exc) { exc_corr = 0.0f; }
		if (compute_y2a) y2a = 0.0f;
		return;
	}
	
	float y = powf(dens, 0.333333333333333333f);  // rho^(1/3)
	float v0 = -0.984745021842697f * y; // -4/3 * (3/PI)^(1/3) * rho^(1/3)

  float ec;
	if (compute_exc) exc_corr = POT_ALPHA * y; // -(3/PI)^(1/3) * rho^(1/3)
	
	switch(gpu_Iexch) {
		case 1:
		{		
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
      if (compute_exc) exc_corr += ec;
		}
		break;
		case 3:
		{
			float rs = POT_GL / y;
			float x1 = sqrtf(rs);
			float Xx = rs + POT_VOSKO_B1 * x1 + POT_VOSKO_C1;
			float t1 = 2 * x1 + POT_VOSKO_B1;
			float t2 = logf(Xx);
			
      ec = POT_VOSKO_A1 * (2 * logf(x1) - t2 + POT_VOSKO_2B1Q * atanf(POT_VOSKO_Q/t1)
        - POT_T4 * (2 * logf(x1 - POT_VOSKO_X0) - t2 + POT_VOSKO_B2X0Q * atanf(POT_VOSKO_Q/t1)));
			
      if (compute_y2a) {
        float vc;
				vc = ec - POT_VOSKO_A16 * x1 * (((POT_VOSKO_B1 * x1 + POT_VOSKO_2C1) / (x1 * Xx)) -
          POT_VOSKO_4B1B1X0 / (t1 * t1 + POT_VOSKO_QSQ) - POT_T4 * (2.0f / (x1 - POT_VOSKO_X0) - t1 / Xx));
				y2a = v0 + vc;
			}
      if (compute_exc) exc_corr += ec;
		}
		break;		
	}
}
