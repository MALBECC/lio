/*
 * function for evaluating exchange correlation density and
 * potential , depends of Iexch, index that says which potential
 * will be used
 * - 1 X alpha
 * - 2 Gunnarson-Lundqvist
 * - 3 Vosko, Wilk and Nusair
 * Exchange part is in all cases the same, for the time being LD
 * Self interaction corrections are used in correlation part
 */

void cal_pot(...) {
	// data X alpha
	const double alpha = -0.738558766382022447;
	const double gl = 0.620350490899400087;
	const double vosko_a1 = 0.03109205, vokso_b1 = 3.72744, vosko_c1 = 12.9352, vosko_x0 = -0.10498, vosko_q = 6.15199066246304849,
		vosko_a16 = 0.005182008333, vosko_a2 = 0.015546025, vosko_b2 = 7.06042, vosko_c2 = 18.0578, vosko_x02 = -0.32500, vosko_q2 = 4.7309269, vosko_a26 = 0.0025910042;
	

	if (dens == 0) {
		ex = 0;
		ec = 0;
		v = 0;
		return;		
	}
	
	double y = pow(dens, 0.333333333333333333);
	double e0 = alpha * y;
	double v0 = 1.33333333333333 * e0;
	
	switch(Iexch) {
		case 1:
		{
			ex = e0;
			ec = 0;
			v = v0;
		}
		break;
		case 2:
		{
			ex = e0;
			rs = gl / y;
			x1 = rs / 11.4; 
			
			if (x1 > 1) // DUDA: aca decia 1.D20
			{
				ec = -0.0333 * (0.5 * x1 - 0.33333333333333);
				vc = 0.0111 * x1 * 0.5;
			}
			else {
				t1 = (1 + x1 * x1 * x1);
				t2 = log(1 + 1 / x1);
				t3 = x1 * x1;
        ec = -0.0333 * (t1 * t2 - t3 + 0.5 * x1 - 0.33333333333333);
        vc = 0.0111 * x1 * (3 * t3 * t2 - t1 / (x1 * (x1 + 1)) - 2 * x1 + 0.5);				
			}

			v = v0 + ec + vc;
		}
		break;
		case 3:
		{
			ex = e0;
			rs = gl / y;
			x1 = sqrt(rs);
			Xx = rs + b1 * x1 + c1;
			Xxo = x0 * x0 + b1 * x0 + c1;
			t1 = 2 * x1 + b1;
			t2 = log(Xx);
			t3 = atan(Q/t1);
			t4 = b1 * x0 / Xxo;
			
      ec = A1 * (2 * log(x1) - t2 + 2 * b1 / Q * t3 - t4 * (2 * log(x1 - x0) - t2 + 2 * (b1 + 2 * x0) / Q * t3));
			t5 = (b1 * x1 + 2 * c1) / x1;
			t6 = x0 / Xxo;
			vc = ec - A16 * x1 * (t5 / Xx - 4 * b1 / (t1 * t1 + Q * Q) * (1 - t6 * (b1 - 2 * x0)) - t4 * (2 / (x1 - x0) - t1 / Xx));
			v = v0 + vc;
		}
		break;		
	}
}
