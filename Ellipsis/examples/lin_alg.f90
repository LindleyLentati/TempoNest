module linearAlgebraUtils
  implicit none
  
contains
  
  !------------------------------------------------------------!
  !               Inverse of a postive definite matrix                   !
  !------------------------------------------------------------!
  subroutine SymPosDefMatInv(mat,ndim)
    integer :: ndim
    double precision, dimension(ndim,ndim) :: mat
    integer :: info
    character ::uplo
    integer :: lda
    integer i,j
    
    lda=ndim
    info=0
    uplo='L'
    
    call dpotrf(uplo,ndim,mat,lda,info)
    
    if(.not.(info .eq. 0)) then
       print *, "@SymPosDefMatInv: Error! Cannot find the inverse"
       stop
    endif
    
    call dpotri(uplo,ndim,mat,lda,info)
    
    if(.not.(info .eq. 0)) then
       print *, "@SymPosDefMatInv: Error! Cannot find the inverse"
       stop
    endif
    
    do i=1,ndim
       do j=1,i-1
          mat(j,i)=mat(i,j)
       enddo
    enddo
  end subroutine SymPosDefMatInv
  
  !------------------------------------------------------------!
  !    eigen decomposition of a postive definite matrix         !
  !------------------------------------------------------------!
  subroutine SymPosDefEigenDecomp(mat,ev,ndim)
    integer :: ndim
    double precision, dimension(ndim,ndim) :: mat,Z
    double precision, dimension(ndim) :: ev
    double precision, dimension(26*ndim) :: WORK
    double precision VL,VU,ABSTOL
    integer :: info
    integer, dimension(2*ndim)  :: ISUPPZ
    integer, dimension(26*ndim) :: IWORK
    character ::uplo,rang,jobs,absfact
    integer lda,ldz,n,m,IL,IU,LWORK,LIWORK
    integer i,j
    double precision dlamch
    
    info=0
    jobs='V'
    rang='A'
    uplo='L'
    n=ndim
    lda=ndim
    VL=0
    VU=0
    IL=0
    IU=0
    absfact='S'
    ABSTOL=dlamch(absfact)
    m=n
    ldz=ndim
    LWORK=26*ndim
    LIWORK=10*ndim
    
    call dsyevr(jobs, rang, uplo, n, mat,lda, VL, VU, IL, IU,ABSTOL,m, ev, &
         Z,LDZ, ISUPPZ,WORK, LWORK,IWORK, LIWORK,info)
    
    if(.not.(info .eq. 0)) then
       print *, "@SymPosDefMatInv: Error! Cannot find the eigen decomposition"
       stop
    endif
    
    do i=1,ndim
       do j=1,ndim
          mat(i,j)=Z(i,j)
       enddo
    enddo
    
  end subroutine SymPosDefEigenDecomp
  
end module linearAlgebraUtils
