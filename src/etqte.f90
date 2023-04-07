!==============================================================================
!      MM\ MMMMMMMM\ MMMMMMMM\ MM\  MMMMMM\  MM\ MMMMMMMM\ MMMMMMMM\ MM\      =
!     MM  |MM  _____|\__MM  __|MM |MM  __MM\ MM |\__MM  __|MM  _____|\MM\     =
!    MM  / MM |         MM |   MM |MM /  MM |MM |   MM |   MM |       \MM\    =
!   MM  /  MMMMM\       MM |   MM |MM |  MM |MM |   MM |   MMMMM\      \MM\   =
!   \MM<   MM  __|      MM |   MM |MM |  MM |MM |   MM |   MM  __|     MM  |  =
!    \MM\  MM |         MM |   MM |MM MM\MM |MM |   MM |   MM |       MM  /   =
!     \MM\ MMMMMMMM\    MM |   MM |\MMMMMM / MM |   MM |   MMMMMMMM\ MM  /    =
!      \__|\________|   \__|   \__| \___MMM\ \__|   \__|   \________|\__/     =
!                                       \___|                                 =
!==============================================================================
!==============================================================================
!========== Energy    Transfer       Quantum        Time   Evolution ==========
!================================ Version 1.0.0 ===============================
!================================== April 2023 ================================
!==============================================================================
!===================== J.-R. Li, Y. Zhai, Z. Qu, and H. Li ====================
!==============================================================================
! This program performs the time evolution based on quantum mechanics.
! The detailed theory can be found in the manual.
! The advantage of this program is that it is free from complex computations,
! and utilize Chebyshev propagator to get the population evolution.
!
! This program is released under MIT licence.
!
! Copyright 2023 the authors of ETQTE
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the “Software”), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

program main
  implicit none
  external :: dsyev ! symmetric matrix diagonalization subroutine from lapack
  character(len=99), external::numberedfile
  double precision, allocatable :: h(:,:), e(:)
  double precision, allocatable :: rho(:,:)
  character(len=99) :: inputfile
  integer :: sys, step, fstate

  ! Input information
  integer :: nstate, nsys, istate, nstep
  character(len=99) :: prefix, suffix
  double precision :: tstep
  common /etqtesetting/ nstate, prefix, suffix, nsys, nstep, tstep, istate

  interface
    subroutine readh(filename, n, heff)
      ! read in matrix file
      implicit none
      character(len=99), intent(in):: filename
      integer, intent(in):: n
      double precision, intent(out) :: heff(:,:)
    end subroutine
  end interface

  write(*,*) "==============================================================================="
  write(*,*) "=      MM\ MMMMMMMM\ MMMMMMMM\ MM\  MMMMMM\  MM\ MMMMMMMM\ MMMMMMMM\ MM\      ="
  write(*,*) "=     MM  |MM  _____|\__MM  __|MM |MM  __MM\ MM |\__MM  __|MM  _____|\MM\     ="
  write(*,*) "=    MM  / MM |         MM |   MM |MM /  MM |MM |   MM |   MM |       \MM\    ="
  write(*,*) "=   MM  /  MMMMM\       MM |   MM |MM |  MM |MM |   MM |   MMMMM\      \MM\   ="
  write(*,*) "=   \MM<   MM  __|      MM |   MM |MM |  MM |MM |   MM |   MM  __|     MM  |  ="
  write(*,*) "=    \MM\  MM |         MM |   MM |MM MM\MM |MM |   MM |   MM |       MM  /   ="
  write(*,*) "=     \MM\ MMMMMMMM\    MM |   MM |\MMMMMM / MM |   MM |   MMMMMMMM\ MM  /    ="
  write(*,*) "=      \__|\________|   \__|   \__| \___MMM\ \__|   \__|   \________|\__/     ="
  write(*,*) "=                                       \___|                                 ="
  write(*,*) "==============================================================================="
  write(*,*) "==============================================================================="
  write(*,*) "=========== Energy    Transfer       Quantum        Time   Evolution =========="
  write(*,*) "================================= Version 1.0.0 ==============================="
  write(*,*) "=================================== April 2023 ================================"
  write(*,*) "==============================================================================="
  write(*,*) "====================== J.-R. Li, Y. Zhai, Z. Qu, and H. Li ===================="
  write(*,*) "==============================================================================="
  call getarg(1, inputfile)
  write(*,*) "The input file is ", inputfile
  write(*,*) "==============================================================================="
  call readinput(inputfile)
  write(*,*) "General settings"
  write(*,*) "==============================================================================="
  write(*,"(a30,5x )", advance="no") "number_of_states" 
  write(*,*)  nstate
  write(*,"(a30,5x )", advance="no") "prefix_of_matrix_files" 
  write(*,*)  prefix
  write(*,"(a30,5x )", advance="no") "suffix_of_matrix_files" 
  write(*,*)  suffix
  write(*,"(a30,5x )", advance="no") "number_of_systems" 
  write(*,*)  nsys
  write(*,"(a30,5x )", advance="no") "number_of_steps" 
  write(*,*)  nstep
  write(*,"(a30,5x )", advance="no") "step_length" 
  write(*,"(f8.5)")  tstep
  write(*,"(a30,5x )", advance="no") "initial_state" 
  write(*,*)  istate
  write(*,*) "==============================================================================="


  tstep=1.d0
  nstep=1000
  istate=1 ! initial state
  allocate(h(nstate,nstate))
  allocate(e(nstate))
  allocate(rho(nstate, 0:nstep))

  write(*,*) "*** Start time evolution..."
  do sys = 1, nsys
    ! write(*,*) "***system ",  sys, " is running."
    call readh(numberedfile(prefix,1,suffix), nstate, h)
    call timeevolution(h, nstate, istate, nstep, tstep, rho)
    if(.true.) then
      open(unit = 1001, file = trim(numberedfile(prefix, sys, suffix))//".etout")
      write(1001, "(a18)", advance = "no") "#             time"
      do fstate = 1, nstate
        write(1001, "(a10,i3)", advance = "no") "rho_", fstate
      enddo
      write(1001, *) 
      do step = 0, nstep
        write(1001,'(f18.5,10(x,f12.6))') step*tstep , rho(:,step)
      enddo
      close(1001)
    endif
  enddo

  write(*,*) "*** End time evolution... "
  write(*,*) "*** Time evolution of each sampled system has been written to .etout files."

end program main

subroutine readinput(fname)
  implicit none
  character(len=99), intent(in) :: fname
  character(len=99) :: tit, val
  character(len=200) :: line
  integer :: io, i

  ! Input information
  integer :: nstate, nsys, istate, nstep
  character(len=99) :: prefix, suffix
  double precision :: tstep
  common /etqtesetting/ nstate, prefix, suffix, nsys, nstep, tstep, istate

  open(unit = 99, file = trim(adjustl(fname)))
  rewind(99)
  do while(.true.)
    read(99,"(a200)",iostat = io) line
    if(io .eq. 0) then
      ! We have some problem with fortran... it cannot read in string containing slash...
      ! This is a possible solution, who knows...
      do i = 1,200
        if(line(i:i) .eq. "/") then
          ! I do not think there is any nut want his/her file name has a char "#"
          ! so it should be safe to use it as a placeholder
          line(i:i) = "#" 
        endif
      enddo
      read(line, *) tit, val
      do i = 1,99
        if(val(i:i) .eq. "#") then
          val(i:i) = "/"
        endif
        ! I think we will not have problem with tit
      enddo
      if (tit .eq. "number_of_states") then
        read(val,*) nstate
      elseif( tit .eq. "prefix_of_matrix_files") then
        prefix = val
      elseif( tit .eq. "suffix_of_matrix_files") then
        suffix = val
      elseif( tit .eq. "number_of_systems") then
        read(val,*) nsys
      elseif( tit .eq. "number_of_steps") then
        read(val,*) nstep
      elseif( tit .eq. "step_length") then
        read(val,*) tstep
      elseif( tit .eq. "initial_state") then
        read(val,*) istate
      endif
    elseif(io .gt. 0) then
      write(*,*) "*** Something went wrong! Check "// fname // "."
      stop
    else
      exit
    endif
  enddo
  close(99)

  return
end subroutine readinput



! this function read in the elements of effective hamiltonian
! here i follow jia-rui's convention.
subroutine readh(filename, n, heff)
  ! read in matrix file
  implicit none
  character(len=99), intent(in):: filename
  integer, intent(in):: n
  double precision, intent(out) :: heff(:,:)
  integer :: i,j

  open(999,file=filename)
  rewind(999)
  do i=1,n
    do j=1,n
      read(999,'(16x,f15.13)') heff(i,j)
    enddo
  enddo
  close(999)     
  return 

end subroutine readh

! this function returns the file name of matrices.
! here i follow jia-rui's convention.
! we assume the input matrices are in atomic units
character(len=99) function numberedfile(prefix, i, suffix)
  implicit none
  character(len=50), intent(in) :: prefix
  integer, intent(in) :: i
  character(len=50), intent(in) :: suffix
  character*20 :: ai
  write(ai, *) i
  numberedfile  = trim(adjustl(prefix)) // trim(adjustl(ai)) // trim(adjustl(suffix))
  return 
end function numberedfile

! the quantum time evolution subroutine
subroutine timeevolution(h, nstate, istate, nstep, tstep, rho)
  implicit none
  double precision, intent(inout) :: h(nstate,nstate) ! in lapack, h will be changed
  integer, intent(in) :: nstate, istate
  integer, intent(in) :: nstep
  double precision, intent(in):: tstep
  double precision, intent(out) :: rho(nstate,0:nstep)
  double precision  :: e(nstate)
  double precision  :: rho0(nstate), param(nstate,nstate*(nstate-1)/2), deltae(nstate*(nstate-1)/2)
  double precision  :: chebyshev_k(nstate*(nstate-1)/2)
  double precision  :: chebyshev_k_1(nstate*(nstate-1)/2)
  double precision  :: omega(nstate*(nstate-1)/2)
  integer :: i, j, fstate, iii, step
  integer :: lwork=-1, info=0, oldnstate=-1
  double precision, allocatable :: work(:)
  save oldnstate, lwork
  character*1 :: jobz = "v", uplo = "u"
  e=0.0d0

  ! some preparation 
  if (oldnstate .ne. nstate) then
    if (allocated(work)) deallocate(work)
    allocate(work(1))
    lwork = -1
    call dsyev(jobz, uplo, nstate, h, nstate, e, work, lwork, info)
    lwork = int(work(1))
    ! write(*,*) "optimized work space length is", lwork
    deallocate(work)
    allocate(work(lwork))
  endif

  call dsyev(jobz, uplo, nstate, h, nstate, e, work, lwork, info)
  ! n.b., h is now the eigenvector of the h passed in!
  ! e is the eigenvalue.

  rho0 = 0.d0
  do fstate = 1, nstate
    do i = 1,nstate
      rho0(fstate) = rho0(fstate) + h(istate, i)**2 * h(fstate,i)**2
    enddo
    iii = 0 
    do i =1, nstate
      do j = i+1, nstate
        iii = iii +1
        param(fstate, iii) = h(istate, i)*h(istate, j)*h(fstate, i)*h(fstate, j)
      enddo
    enddo
  enddo
  param = 2.d0 * param ! we move the 2 here
  iii = 0
  do i =1, nstate
    do j = i+1, nstate
      iii = iii +1
      deltae(iii) = e(j) - e(i)
    enddo
  enddo

  ! n.b. 
  ! here i use the chebyshev propagator, because the cosine function is extremely slow.
  ! ref: chen & guo, comput. phys. comm. 119, 19-31 (1999)
  omega = dcos(tstep*deltae) 

  ! the following is for k = step = 0...
  chebyshev_k = 1.d0
  chebyshev_k_1 = omega ! t_{-1}(omega)
  do step=0,nstep

    rho(:, step) = matmul(param, chebyshev_k) + rho0 
    ! the coefficient 2 have been included in param
    chebyshev_k_1 = 2.d0*omega*chebyshev_k -chebyshev_k_1
    chebyshev_k = chebyshev_k + chebyshev_k_1
    chebyshev_k_1 = chebyshev_k - chebyshev_k_1
    chebyshev_k = chebyshev_k - chebyshev_k_1

  enddo

  return
end subroutine

