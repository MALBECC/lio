interface read_matrix
   module procedure read_matrix_cd
   module procedure read_matrix_cs
   module procedure read_matrix_d
   module procedure read_matrix_s
end interface read_matrix

interface read_sqmatrix
   module procedure read_sqmatrix_cd
   module procedure read_sqmatrix_cs
   module procedure read_sqmatrix_d
   module procedure read_sqmatrix_s
end interface read_sqmatrix

interface write_matrix
   module procedure write_matrix_cd
   module procedure write_matrix_cs
   module procedure write_matrix_d
   module procedure write_matrix_s
end interface write_matrix

interface write_sqmatrix
   module procedure write_sqmatrix_cd
   module procedure write_sqmatrix_cs
   module procedure write_sqmatrix_d
   module procedure write_sqmatrix_s
end interface write_sqmatrix
