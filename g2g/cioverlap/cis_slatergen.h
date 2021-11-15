void do_cis_slater(bool uhf, int ncore0, int (&nocc0)[2], int (&nvirt0)[2],
                   int inactive_occ, int inactive_virt, bool include_occ = 1,
                   bool number_occ = 1) {
  // cout << "Is this UHF run? "<<uhf<<endl;
  // cout <<"core occa occb virta virtb freeze disc "<<ncore0<<" "<<nocc0[0]<<"
  // "<<nocc0[1]<<" "<<nvirt0[0]<<" "<<nvirt0[1]<<" "<<inactive_occ<<"
  // "<<inactive_virt<<endl;

  int ncore = ncore0 + inactive_occ;
  int nocc[2];
  nocc[0] = nocc0[0] - inactive_occ;
  nocc[1] = nocc0[1] - inactive_occ;
  int nvirt[2];
  nvirt[0] = nvirt0[0] - inactive_virt;
  nvirt[1] = nvirt0[1] - inactive_virt;

  // cout <<"nocc_eff = "<<nocc[0] <<" "<<nocc[1]<<endl;
  // cout <<"virt_eff = "<<nvirt[0] <<" "<<nvirt[1]<<endl;

  slaterdet tmp((include_occ ? 2 * ncore : 0) + nocc[0] +
                nocc[1]);  // general uhf

  int slaterfile = open("slaterfile", O_WRONLY | O_CREAT | O_LARGEFILE, 0777);
  if (slaterfile < 0) {
    perror("cannot open");
    laerror("IO error");
  }

  if (include_occ)
    for (int iocc = 1; iocc <= ncore;
         ++iocc)  // generate fixed core occupation alternating alpha and beta
    {
      tmp[2 * iocc - 2] = iocc;
      tmp[2 * iocc - 1] = -iocc;
    }

  int shiftcore = include_occ ? ncore : 0;
  int numcore = number_occ ? ncore : 0;

  // fermi vacuum determinant

  NRVec<int> positions[2];
  if (uhf) {
    positions[0].resize(nocc[0]);
    positions[1].resize(nocc[1]);
  }

  int noccmin = nocc[0];
  if (nocc[1] < noccmin) noccmin = nocc[1];
  int noccmax = nocc[0];
  if (nocc[1] > noccmax) noccmax = nocc[1];
  // generate part with alternating alpha and beta
  for (int iocc = 1; iocc <= noccmin; ++iocc) {
    int pos = 2 * shiftcore + 2 * iocc - 2;
    tmp[2 * shiftcore + 2 * iocc - 2] = iocc + numcore;
    tmp[pos + 1] = -iocc - numcore;
    if (uhf) {
      positions[0][iocc - 1] = pos;
      positions[1][iocc - 1] = pos + 1;
    }
  }
  // and add extra alpha or beta orbitals
  if (nocc[0] != nocc[1]) {
    int spin = nocc[0] > nocc[1] ? 1 : -1;
    for (int iocc = noccmin + 1; iocc <= noccmax; ++iocc) {
      int pos = 2 * shiftcore + 2 * noccmin + iocc - noccmin - 1;
      tmp[pos] = spin * (iocc + numcore);
      if (uhf) positions[spin == 1 ? 0 : 1][iocc - 1] = pos;
    }
  }

  // store the vacuum
  tmp.put(slaterfile, false);

  // now generate all possible excitations

  // order of indices will be occ,virt
  // notice that we make replacement in place - parity of permutation to bring a
  // given spinorb in front and after excitation to bring it back compensates

  if (uhf) {
    for (int spin = 0; spin <= 1;
         ++spin)  // first all alpha, then all beta excitations
    {
      for (int iocc = 1; iocc <= nocc[spin]; ++iocc)
        for (int ivirt = nocc[spin] + 1; ivirt <= nocc[spin] + nvirt[spin];
             ++ivirt) {
          // use stored positions to make replacements
          int pos = positions[spin][iocc - 1];
          int original = tmp[pos];
          tmp[pos] = (spin ? -1 : 1) * (ivirt + numcore);
          tmp.put(slaterfile, false);
          tmp[pos] = original;
        }
    }
  } else  // not uhf ... keep original order of determinants for full backward
          // compatibility
  {
    for (int iocc = 1; iocc <= nocc[0]; ++iocc)
      for (int ivirt = nocc[0] + 1; ivirt <= nocc[0] + nvirt[0]; ++ivirt) {
        tmp[2 * shiftcore + 2 * iocc - 2] = ivirt + numcore;
        tmp.put(slaterfile, false);  // alpha excitation
        tmp[2 * shiftcore + 2 * iocc - 2] = iocc + numcore;

        tmp[2 * shiftcore + 2 * iocc - 1] = -ivirt - numcore;
        tmp.put(slaterfile, false);  // beta excitation
        tmp[2 * shiftcore + 2 * iocc - 1] = -iocc - numcore;
      }
  }

  // cout <<"Number of determinants generated =
  // "<<1+nocc[0]*nvirt[0]+nocc[1]*nvirt[1]<<endl;

  close(slaterfile);
}
