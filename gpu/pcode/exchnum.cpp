
// Energia de Intercambio y Correlacion
// Input:  posicion de los atomos ('pos_atomos'), cantidad de capas segun tipo de atomos ('capas'), tipo de los atomos ('tipos'),
//         el tipo de grilla a usar ('tipo_grilla'), la posicion relativa al centro de los puntos de
//         la grilla ('pos_puntos'), arreglo Rm, un valor por cada tipo de atomo ('Rm'), la cantidad de funciones s, p y d ('cant_funcs'), M ('m'), Nucleo del orbital (numero
//         de atomo en el centro de la funcion) ('nuc')
// Output: energia de intercambio y correlacion

double exch_corr(array<double3> pos_atomos, array<int> capas, array<int> tipos, TipoGrilla tipo_grilla, array<double3> pos_puntos,
								 array<double> Rm, int3 cant_funcs, double m, array<double> nuc, array<int> contracciones, array<double> factor_a, array<double> factor_c)
{
	double exc_total = 0;
  double ecorr_total = 0;
	double ss = 0; // ?
	int atomos = pos_atomos.size;
	double3 distancias_al_resto(atomos, atomos);
	double3 distancias_punto_atomos(atomos);

	//  do 43 l=1,3
	//  Ll(l)= l*(l-1)/2
	
	// distancias todos entre todos
	// DUDA: esto tiene sentido recalcularlo aca? --> en principio se va afuera
	for each atomo_i in 1,atomos do {
		for each atomo_j in 1,atomos do {
			distancias_al_resto(atomo_i, atomo_j) = dist(pos_atomos(atomo_i), pos_atomos(atomo_j))
		}
	}
	
	int puntos_grilla = (tipo_grilla == TipoGrilla1 ? 116 : 194);

	for each atomo_i in 1,atomos do {
		tipo_atomo_i = tipos(atomo_i);
		capas_tomo_i = capas(tipo_atomo_i);
		
		for each capa in 1,capas_atomo_i {
			double tmp0 = (PI / (capas_atomo_i + 1));
			double tmp1 = tmp0 * capa;
			double x = cos(tmp1)
			double w = tmp0 * abs(sin(tmp1))
			double r1 = Rm(tipo_atomo_i) * (1 + x) / (1 - x)
			double wrad = w * (r1 * r1) * Rm(tipo_atomo_i) * 2 / ((1 - x) * (1 - x));
										
			for each punto in 1,puntos_grilla {
				double posicion_punto = pos_atomos(atomo_i) + r1 * (tipo_grilla == TipoGrilla1 ? pos_puntos(punto) : pos_puntos_alt(punto));
				double peso_integracion = wrad * (tipo_grilla == TipoGrilla1 ? wang(punto) : wang_alt(punto));

				double3 distancias_punto_atomos(atomos);
				for each atomo_j in 1,atomos {
					distancias_punto_atomos(atomo_j) = dist_cuadrada(posicion_punto, pos_atomos(atomo_j))
				}

				double exc_actual, ecorr_actual;
				calcular_densidad(&exc_actual, &ecorr_actual, cant_funcs, m, distancias_punto_atomos, nuc, contracciones, factor_a, factor_c, posicion_punto, posiciones_atomos, ...);
				calcular_pot(...);

				// integracion numerica
				double P_total = 0;
				
				for each atomo_j in 1,atomos {
					P(atomo_j) = 1

					for each atomo_k in 1,atomos {
						if (atomo_k == atomo_j) continue;
						double u = sqrt(distancias_punto_atomos(atomo_j)) - sqrt(distancias_punto_atomos(atomo_k));
						double x = Rm(tipos(atomo_j)) / Rm(tipos(atomo_k));
						double x1 = (x - 1) / (x + 1);
						double Aij = x1 / (x1 * x1 - 1);
						u += aij * (1 + u * u);

						double p1 = 1.5 * u - 0.5 * u ** 3;
						double p2 = 1.5 * p1 - 0.5 * p1 ** 3;
						double p3 = 1.5 * p2 - 0.5 * p2 ** 3;
						double p4 = 1.5 * p3 - 0.5 * p3 ** 3;
						double p5 = 1.5 * p4 - 0.5 * p4 ** 3;
						double s = 0.5 * (1 - p3);
						P(atomo_j) *= s;
					}
					P_total += P(atomo_j);
				}

				tmp0 = (P(atomo_i) / P_total) * dens * peso_integracion;
				exc_total += tmp0 * exc_actual;
				ecorr_total += tmp0 * ecorr_actual;
				ss += tmp0;
				// Npt++;	// DUDA: esta no esta usada en fortran
			} // puntos grilla
		} // capas
	}	// atomos

	return exc_total + ecorr_total;
}
