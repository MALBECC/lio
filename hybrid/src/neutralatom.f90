        module neutralatom 

        integer, parameter       :: nemax = 100  ! max. number of chem. species
        double precision, save   :: delta(nemax) ! Dist between points in table
        private nemax,delta
        end module
