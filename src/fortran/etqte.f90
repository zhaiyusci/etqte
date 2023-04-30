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

module etqte
  implicit none
  ! Input information
  integer :: nstate, nsys, istate, nstep
  character(len=99) :: prefix, suffix
  double precision :: tstep

  ! Workspace
  double precision, allocatable :: h(:,:,:), e(:,:), psi(:,:,:)
  double precision, allocatable :: rho(:,:,:), rho0(:,:)
  double precision, allocatable :: averrho(:,:), averrho0(:)

contains

  subroutine writetimeseries(sys, instval, infval)
    implicit none

    integer, intent(in) :: sys
    double precision, intent(in) :: instval(nstate, 0:nstep), infval(nstate)
    integer :: fstate, step, state

    ! write(timeformat,"(a, i5, a)")'(f18.5,', nstate, '(x,f12.6))'

    open(unit = 1001, file = trim(numberedfile(sys))//".etout")
    write(1001, "(a18)", advance = "no") "#             time"
    do fstate = 1, nstate
      write(1001, "(a10,i3)", advance = "no") "phi_", fstate
    enddo
    write(1001, *) 
    do step = 0, nstep
      write(1001,"(f18.5)", advance = "no") step*tstep 
      do state = 1, nstate
        write(1001,"(f13.5)", advance = "no") instval(state,step)
      enddo
      write(1001,*)
    enddo
    ! We use a long enough time to present an infinite long time
    write(1001,"(f18.5)", advance = "no") nstep*tstep*2
    do state = 1, nstate
      write(1001,"(f13.5)", advance = "no") infval(state)
    enddo
    write(1001,*)
    close(1001)
    return 
  end subroutine writetimeseries

  ! Read the input file, a sample of which is provided.
  subroutine readsettings(fname)
    implicit none
    character(len=99), intent(in) :: fname
    character(len=99) :: tit, val
    character(len=200) :: line
    integer :: io, i

    open(unit = 99, file = trim(adjustl(fname)))
    rewind(99)
    do while(.true.)
      read(99,"(a200)",iostat = io) line
      if(io .eq. 0) then
        ! We have some problem with fortran... it cannot read in string containing slash...
        ! This is a possible solution, who knows...
        do i = 1,200
          if(line(i:i) .eq. "/") then
            ! I do not think there is any nut wanting the file name has a char "#",
            ! which is the comment char in most of script languages, includeing bash, zsh, python, etc.
            ! So it should be safe to use it as a placeholder.
            line(i:i) = "#" 
          endif
        enddo
        read(line, *) tit, val
        do i = 1,99
          if(val(i:i) .eq. "#") then
            val(i:i) = "/"
          endif
          ! I think we will not have problem with tit, we set it.
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
        elseif( tit(1:1) .eq. "#") then
          read(val,*) ! do nothing
        else
          write(*,*) "*** Unknown title ", tit 
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
  end subroutine readsettings

  ! this function read in the elements of effective hamiltonian
  ! here i follow jia-rui's convention.
  subroutine readhamiltonian(sys)
    ! read in matrix file
    implicit none
    integer, intent(in) :: sys
    integer :: i,j

    open(999,file=numberedfile(sys))
    rewind(999)
    do i=1,nstate
      do j=1,nstate
        read(999,'(16x,f15.13)') h(i,j,sys)
      enddo
    enddo
    close(999)     
    return 

  end subroutine readhamiltonian

  ! this function returns the file name of matrices.
  ! here i follow jia-rui's convention.
  ! we assume the input matrices are in atomic units
  character(len=99) function numberedfile(i)
    implicit none
    integer, intent(in) :: i
    character*20 :: ai
    write(ai, *) i
    numberedfile  = trim(adjustl(prefix)) // trim(adjustl(ai)) // trim(adjustl(suffix))
    return 
  end function numberedfile

  ! the quantum time evolution subroutine
  subroutine timeevolution(sys)
    implicit none
    integer , intent(in) :: sys
    double precision  :: param(nstate,nstate*(nstate-1)/2), deltae(nstate*(nstate-1)/2)
    double precision  :: chebyshev_k(nstate*(nstate-1)/2)
    double precision  :: chebyshev_k_1(nstate*(nstate-1)/2)
    double precision  :: omega(nstate*(nstate-1)/2)
    integer :: i, j, fstate, iii, step
    integer :: lwork=-1, info=0, oldnstate=-1
    double precision, allocatable :: work(:)
    save oldnstate, lwork
    character*1 :: jobz = "v", uplo = "u"
    e(:,sys)=0.0d0

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

    ! Have a copy of h so lapack does not distroy them
    psi(:, :, sys)=h(:,:,sys)

    rho0(:, sys) = 0.d0

    call dsyev(jobz, uplo, nstate, psi(:,:,sys), nstate, e(:, sys), work, lwork, info)
    ! N.B.
    ! h is now the eigenvectors of the h passed in!
    ! e is the eigenvalue.

    do fstate = 1, nstate
      do i = 1,nstate
        rho0(fstate,sys) = rho0(fstate,sys) + psi(istate, i, sys)**2 * psi(fstate,i, sys)**2
      enddo
      iii = 0 
      do i =1, nstate
        do j = i+1, nstate
          iii = iii +1
          param(fstate, iii) = psi(istate, i, sys)*psi(istate, j,sys)*psi(fstate, i,sys)*psi(fstate, j, sys)
        enddo
      enddo
    enddo
    param = 2.d0 * param ! we move the 2 here
    iii = 0
    do i =1, nstate
      do j = i+1, nstate
        iii = iii +1
        deltae(iii) = e(j, sys) - e(i, sys)
      enddo
    enddo

    ! N.B. 
    ! Here I use the Chebyshev propagator, because the cosine function is extremely slow.
    ! Ref: Chen & Guo, Comput. Phys. Comm. 119, 19-31 (1999)
    omega = dcos(tstep*deltae) 

    ! the following is for k = step = 0...
    chebyshev_k = 1.d0
    chebyshev_k_1 = omega ! t_{-1}(omega)
    do step=0,nstep

      rho(:, step, sys) = matmul(param, chebyshev_k) + rho0(:,sys)
      ! the coefficient 2 have been included in param
      chebyshev_k_1 = 2.d0*omega*chebyshev_k -chebyshev_k_1
      chebyshev_k = chebyshev_k + chebyshev_k_1
      chebyshev_k_1 = chebyshev_k - chebyshev_k_1
      chebyshev_k = chebyshev_k - chebyshev_k_1

    enddo


    return
  end subroutine

  subroutine printbanner(io)
    implicit none
    integer, intent(in) :: io
    write(io,*) "==============================================================================="
    write(io,*) "=      MM\ MMMMMMMM\ MMMMMMMM\ MM\  MMMMMM\  MM\ MMMMMMMM\ MMMMMMMM\ MM\      ="
    write(io,*) "=     MM  |MM  _____|\__MM  __|MM |MM  __MM\ MM |\__MM  __|MM  _____|\MM\     ="
    write(io,*) "=    MM  / MM |         MM |   MM |MM /  MM |MM |   MM |   MM |       \MM\    ="
    write(io,*) "=   MM  /  MMMMM\       MM |   MM |MM |  MM |MM |   MM |   MMMMM\      \MM\   ="
    write(io,*) "=   \MM<   MM  __|      MM |   MM |MM |  MM |MM |   MM |   MM  __|     MM  |  ="
    write(io,*) "=    \MM\  MM |         MM |   MM |MM MM\MM |MM |   MM |   MM |       MM  /   ="
    write(io,*) "=     \MM\ MMMMMMMM\    MM |   MM |\MMMMMM / MM |   MM |   MMMMMMMM\ MM  /    ="
    write(io,*) "=      \__|\________|   \__|   \__| \___MMM\ \__|   \__|   \________|\__/     ="
    write(io,*) "=                                       \___|                                 ="
    write(io,*) "==============================================================================="
    write(io,*) "==============================================================================="
    write(io,*) "=========== Energy    Transfer       Quantum        Time   Evolution =========="
    write(io,*) "================================= Version 1.0.0 ==============================="
    write(io,*) "=================================== April 2023 ================================"
    write(io,*) "==============================================================================="
    write(io,*) "====================== J.-R. Li, Y. Zhai, Z. Qu, and H. Li ===================="
    write(io,*) "==============================================================================="
  end subroutine printbanner

  subroutine printsettings(io)
    implicit none
    integer, intent(in) :: io
    write(io,*) "General settings"
    write(io,*) "==============================================================================="
    write(io,"(a30,5x )", advance="no") "number_of_states" 
    write(io,*)  nstate
    write(io,"(a30,5x )", advance="no") "prefix_of_matrix_files" 
    write(io,*)  prefix
    write(io,"(a30,5x )", advance="no") "suffix_of_matrix_files" 
    write(io,*)  suffix
    write(io,"(a30,5x )", advance="no") "number_of_systems" 
    write(io,*)  nsys
    write(io,"(a30,5x )", advance="no") "number_of_steps" 
    write(io,*)  nstep
    write(io,"(a30,5x )", advance="no") "step_length" 
    write(io,"(f12.5)")  tstep
    write(io,"(a30,5x )", advance="no") "initial_state" 
    write(io,*)  istate
    write(io,*) "==============================================================================="
  end subroutine printsettings

  subroutine allocatearrays()
    implicit none
    allocate(h(nstate,nstate,nsys))
    allocate(psi(nstate,nstate,nsys))
    allocate(e(nstate,nsys))
    allocate(rho(nstate, 0:nstep, nsys))
    allocate(rho0(nstate, nsys))
    allocate(averrho(nstate, 0:nstep))
    allocate(averrho0(nstate))
  end subroutine allocatearrays


end module etqte


program main
  use etqte
  implicit none
  external :: dsyev ! symmetric matrix diagonalization subroutine from lapack
  character(len=99) :: inputfile
  integer :: sys, i ,j
  DOUBLE PRECISION :: t1, t2

  call printbanner(6)


  call getarg(1, inputfile)
  write(*,*) "*** The input file is ", inputfile
  write(*,*) "==============================================================================="
  call readsettings(inputfile)
  call printsettings(6)

  call allocatearrays()

  write(*,*) "*** Read Hamiltonian matrices in..."
  do sys = 1, nsys
    call readhamiltonian(sys)
  enddo
  write(*,*) "*** Done."

  ! Apparently the matrices are in eV, Errrrr...
  ! Turn hamiltonian in eV to hartree. Idiot.
  h = h / 27.21138602d0 

  write(*,*) "*** Start time evolution..."
  call CPU_TIME(t1)
  do sys = 1, nsys
    call timeevolution(sys)
  enddo
  call CPU_TIME(t2)
  write(*,*) "*** End time evolution... "
  write(*,"(3f12.6)") t1, t2, t2-t1
  do sys = 1, nsys
    call writetimeseries(sys, rho(:,:,sys), rho0(:,sys))
  enddo
  write(*,*) "*** Time evolution of each sampled system has been written to .etout files."

  if(.false.) then
    do sys = 1, nsys
      write(*,*) "Hamiltonian of system ",sys
      do i = 1,nstate
        ! Anyway, h is a real symmetric matrix
        do j = 1,nstate
          write(*,"(f12.6)", advance = "no") h(j, i, sys) 
        enddo
        write(*,*)
      enddo
      write(*,*) 

      write(*,*) "Eigenvectors of the system ", sys, " (phi representation)"
      do i = 1,nstate
        do j = 1,nstate
          write(*,"(f12.6)", advance = "no") psi(i, j, sys) 
        enddo
        write(*,*)
      enddo
      write(*,*) 

      if(.true.) then
      endif
    enddo
  endif


  ! Ensemble average...
  write(*,*) "*** Start ensemble avearge... "
  averrho = sum(rho, 3)/nsys
  averrho0 = sum(rho0, 2)/nsys
  call writetimeseries(0, averrho(:,:), averrho0(:))
  write(*,*) "*** Time evolution of the whole ensemble has been written to the zeroth .etout files."
  write(*,*) "*** Print out the adiabatic energies... "

  do sys = 1, nsys
    do i = 1, nstate
      write(1002, "(f13.5)", advance = "no") e(i, sys)
    enddo
    write(1002, "(f13.5)") 
  enddo


end program main

