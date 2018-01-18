!========================================================
!= Sample program using the f90 fdf module : CJP 1/3/99 =
!========================================================
Program sample

  Use fdf

  Implicit None

  Integer  maxa
  Parameter ( maxa = 100 )

  Character         fname*20, symbol(maxa)*2
  Integer           i, ia, isa(maxa), na, na_default, iblk
  Real              wmix

  Double Precision  factor, xa(3,maxa), cutoff, phonon_energy

  Logical doit, debug

  Call fdf_init('sample.fdf','sample.out')

  If (fdf_defined('new-style')) Write(6,*) 'New-style stuff'

  na_default = 0
  na = fdf_integer('NumberOfAtoms', na_default )
  Write(6,*) 'examples: na =', na

  fname = fdf_string('NameOfFile','calimero')
  Write(6,*) fname

  cutoff = fdf_physical('MeshCutoff',8.d0,'Ry')
  Write(6,*) cutoff

  phonon_energy = fdf_physical('phonon-energy',0.01d0,'eV')
  Write(6,*) phonon_energy

  i = fdf_integer('SomeInt',34)
  Write(6,*) i

  wmix = fdf_single('WmixValue',0.55)
  Write(6,*) wmix

  factor = fdf_double('FactorValue',1.d-10)
  Write(6,*) factor

  debug = fdf_boolean('Debug',.True.)
  Write(6,*) debug

  doit = fdf_boolean('DoIt',.False.)
  Write(6,*) doit

  If (fdf_block('AtomicCoordinatesAndAtomicSpecies',iblk)) Then
     Do ia = 1,na
        Read(iblk,*) (xa(i,ia),i=1,3), isa(ia)
     Enddo
  Endif

  If (fdf_block('AtomicSymbolsAndAtomicCoordinates',iblk)) Then
     Do ia = 1,na
        Read(iblk,*) symbol(ia), (xa(i,ia),i=1,3)
     Enddo
  Endif

  Do ia = 1,na
     Write(6,*) (xa(i,ia),i=1,3)
  Enddo

  If (fdf_block('AtomicInfo',iblk)) Then
     Do ia = 1,na
        Read(iblk,*) (xa(i,ia),i=1,3)
     Enddo
  Endif

  Do ia = 1,na
     Write(6,*) (xa(i,ia),i=1,3)
  Enddo

End Program sample
