!MS$DECLARE
module      comp_parameter_mod
    save
    integer(8),parameter   ::  time_max=nint(8e14_8,8)!80000000 !0!nint(1d6)
    integer,parameter   ::  time_mul=1!nint(1d6)
    integer,parameter   ::  halfsize=66
    integer,parameter   ::  nroot   =5
    integer             ::  wo_cnt  =0
endmodule   comp_parameter_mod

module      deformations_mod
    use comp_parameter_mod
    save
    real(8),dimension(-halfsize:halfsize,-halfsize:halfsize) :: &
    displacement_x,displacement_y,strain_xx,strain_yy!,strain_y
endmodule   deformations_mod

module      diffusion_mod
    use comp_parameter_mod
    save
    real(8),dimension(-halfsize:halfsize,-halfsize:halfsize) :: d1_s,d2_s!,dx_s,dy_s
endmodule   diffusion_mod

module      concentration_mod
    use comp_parameter_mod
    save
    real(8),dimension(-halfsize:halfsize,-halfsize:halfsize) :: c_iia,cdx,cdy,c_add !concentration of carbon
endmodule   concentration_mod

module      phys_constants_mod
    save
    real(8),parameter   ::  lattice_parameter   = 2.86d0
    real(8),parameter   ::  volume_per_atom     = lattice_parameter**3*5d-1
    real(8),parameter   ::  diffusion_const     = 5d-8*1d20
    real(8),parameter   ::  barrier             = 0.86d0 !ev
    real(8),parameter   ::  boltzmann_const     = 8.6173303d-5 !
    real(8),parameter   ::  carbon_concentration= 5d-4 !initial concentration

!    real(8),parameter   ::  time_step_duration  = 16d-13*2d-2!/14d1
!    real(8),parameter   ::  temperature         = 900d0+273d0

!    real(8),parameter   ::  time_step_duration  = 16d-13*1d-4
!    real(8),parameter   ::  temperature         = 823d0

    real(8),parameter   ::  time_step_duration  = 37d-13*3d-2* 12d-2!* 1d-2
    real(8),parameter   ::  temperature         = 1500d0+273d0

!    real(8),parameter   ::  time_step_duration  = 21d-13*1d-3
!    real(8),parameter   ::  temperature         = 800d0+273d0
    real(8),parameter   ::  k1                  = 5.6693838153357401d0 !halfackland
    real(8),parameter   ::  k2                  =-4.0988435809159078d0 !halfackland
    real(8),parameter   ::  k3                  = 2.8860853075335791d0 !halfackland
!    real(8),parameter   ::  k1                  = 7.3065386856839769d0 !halfackland
!    real(8),parameter   ::  k2                  =-3.4342361908421379d0 !halfackland
!    real(8),parameter   ::  k3                  = 3.3308963651920958d0 !halfackland
    real(8),parameter   ::  h_x                 =lattice_parameter!*5d-1
    real(8),parameter   ::  h_y                 =h_x!lattice_parameter!*5d-1
endmodule   phys_constants_mod

PROGRAM    COTTRELL_100_HA_SIMPLE !___________________
    use omp_lib
    use comp_parameter_mod

    integer i,j
    integer(8) i_time
!    real ancil1
  !STOP" nothing is done "
    call    po_info
    call    obtain____matrices_e   !; zagruzitq napriazxenija
    call    wo_matrices_e
    !STOP "obtain matrices deformation"
    call    calculate_matrices_d   !; rasscxitatq matricy D
        !CALL    WO_MATRICES
    !STOP "obtain matrices diffusion"
    call    initiate__matrices_c   !; rasscxitatq matricy D
    !STOP "obtain matrices diffusion"
        CALL    CALCULATE_CONCENTRATION_D_PRODUCT !WRAP IN THE PARTIALS
    print'(A,$)',"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
    do i_time=1_8,1_8+time_max
        !EXIT
        call    calculate_concentration_d_product !wrap in the partials
        call    calculate_concentration__addition(i_time) !wrap in the partials
        call    apply_concentration_addition
        if (mod(i_time-1,1500000).eq.0) print*," ",i_time," of ",time_max
        if (mod(i_time-1,0050000).eq.0) print'(A,$)',"."
        if (mod(i_time  ,0050000).eq.1) call wo_matrices(i_time)
    enddo


ENDPROGRAM COTTRELL_100_HA_SIMPLE !___________________

subroutine      obtain____matrices_e
    use comp_parameter_mod
    use deformations_mod
    use phys_constants_mod
    integer i_x,i_y
    real(8),parameter :: reci_a=1d0/lattice_parameter
    real(8) strain_xx_symmetric(-halfsize:halfsize,-halfsize:halfsize)
    !REAL(8) E11,E22
!    open(301, file = 'Matrix-stage11-z0-Ux-100.txt')
!    open(302, file = 'Matrix-stage11-z0-Uy-100.txt')
    open(301, file = 'Matrix-stage11-un-Ux-100.txt')
    open(302, file = 'Matrix-stage11-un-Uy-100.txt')
!    open(301, file = 'Matrix-stage01-z0-Ux-100.txt')
!    open(302, file = 'Matrix-stage01-z0-Uy-100.txt')
!    open(301, file = 'Matrix-stage02-z0-Ux-100.txt')
!    open(302, file = 'Matrix-stage02-z0-Uy-100.txt')
    do  i_y=halfsize,-halfsize,-1
        !do  i_x=-halfsize,halfsize
            read(301,401) displacement_x(-halfsize:halfsize,i_y)
            read(302,401) displacement_y(-halfsize:halfsize,i_y)
            !i_y is slow changing here, as it is in wo_matrices
        !enddo
    enddo
    !displacement_y(0,2:halfsize)=displacement_y(0,2:halfsize)+lattice_parameter*5d-1
      !displacement_y(0,1:halfsize)=displacement_y(0,1:halfsize)+lattice_parameter*5d-1 !initial
    !displacement_y(0,0:halfsize)=displacement_y(0,0:halfsize)+lattice_parameter*5d-1
    do  i_y=-halfsize+1,halfsize-1
        do  i_x=-halfsize+1,halfsize-1
            strain_xx(i_x,i_y)=(displacement_x(i_x+1,i_y  )-displacement_x(i_x-1,i_y  ))*&
            reci_a*5d-1
            strain_yy(i_x,i_y)=(displacement_y(i_x  ,i_y+1)-displacement_y(i_x  ,i_y-1))*&
            reci_a*5d-1
        enddo
    enddo
    !STOP"continue from there"
!    do  i_y= 1,halfsize-1
    do  i_y= -halfsize+1,-1
        i_x= 1        !?
        strain_xx(i_x,i_y)=(displacement_x(i_x  ,i_y  )-displacement_x(i_x+1,i_y  ))*&
        reci_a
        strain_xx(i_x,i_y)=strain_xx(i_x+1,i_y)*2d0-strain_xx(i_x+2,i_y)        !?
!        strain_yy(i_x,i_y)=(displacement_y(i_x  ,i_y  )-displacement_y(i_x+1,i_y  ))*&
!        reci_a
!        strain_yy(i_x,i_y)=strain_yy(i_x+1,i_y)*2d0-strain_yy(i_x+2,i_y)
    enddo
    !do  i_y= 1,halfsize-1
    do  i_y= -halfsize+1,-1
        i_x=-1
        strain_xx(i_x,i_y)=(displacement_x(i_x  ,i_y  )-displacement_x(i_x-1,i_y  ))*&
        reci_a
        strain_xx(i_x,i_y)=strain_xx(i_x,i_y)*2d0-strain_xx(i_x-1,i_y)
!        strain_yy(i_x,i_y)=(displacement_y(i_x  ,i_y  )-displacement_y(i_x-1,i_y  ))*&
!        reci_a
!        strain_yy(i_x,i_y)=strain_yy(i_x,i_y)*2d0-strain_yy(i_x-1,i_y)
    enddo

    !do  i_y= 1,halfsize-1
    do  i_y= -halfsize+1,-1
        i_x= 0
!        PRINT*, i_y, STRAIN_XX(I_X-2,I_Y), STRAIN_XX(I_X-1,I_Y), STRAIN_XX(I_X ,I_Y), STRAIN_XX(I_X+1,I_Y), STRAIN_XX(I_X+2,I_Y)
!        PRINT*, i_y, displacement_x(I_X-2,I_Y), displacement_x(I_X-1,I_Y), displacement_x(I_X ,I_Y), &
!          displacement_x(I_X+1,I_Y), displacement_x(I_X+2,I_Y)
        !exit
        strain_xx(i_x,i_y)=(strain_xx(i_x-1,i_y)+strain_xx(i_x+1,i_y))*5d-1
        !PRINT*, i_y, STRAIN_XX(I_X-2,I_Y), STRAIN_XX(I_X-1,I_Y), STRAIN_XX(I_X ,I_Y), STRAIN_XX(I_X+1,I_Y), STRAIN_XX(I_X+2,I_Y)
        !strain_yy(i_x,i_y)=(strain_yy(i_x-1,i_y)+strain_yy(i_x+1,i_y))*5d-1
    enddo

    do  i_y=-halfsize+1,halfsize-1
        i_x= halfsize
        strain_xx(i_x,i_y)=strain_xx(i_x-1,i_y)*2d0-strain_xx(i_x-2,i_y)
        strain_yy(i_x,i_y)=strain_yy(i_x-1,i_y)*2d0-strain_yy(i_x-2,i_y)
    enddo
    do  i_y=-halfsize+1,halfsize-1
        i_x=-halfsize
        strain_xx(i_x,i_y)=strain_xx(i_x+1,i_y)*2d0-strain_xx(i_x+2,i_y)
        strain_yy(i_x,i_y)=strain_yy(i_x+1,i_y)*2d0-strain_yy(i_x+2,i_y)
    enddo

    do  i_x=-halfsize ,halfsize
        i_y= halfsize
        strain_xx(i_x,i_y)=strain_xx(i_x,i_y-1)*2d0-strain_xx(i_x,i_y-2)
        strain_yy(i_x,i_y)=strain_yy(i_x,i_y-1)*2d0-strain_yy(i_x,i_y-2)
    enddo
    do  i_x=-halfsize ,halfsize
        i_y=-halfsize
        strain_xx(i_x,i_y)=strain_xx(i_x,i_y+1)*2d0-strain_xx(i_x,i_y+2)
        strain_yy(i_x,i_y)=strain_yy(i_x,i_y+1)*2d0-strain_yy(i_x,i_y+2)
    enddo
    strain_xx_symmetric(-halfsize:halfsize,-halfsize:halfsize)=&
    strain_xx(-halfsize:halfsize,-halfsize:halfsize)+&
    strain_xx(halfsize:-halfsize:-1,-halfsize:halfsize)
    strain_xx=5d-1*strain_xx_symmetric
    close(301)
    close(302)
401 format(200(1x,E11.4))
    !    DO  I_Y=-HALFSIZE,HALFSIZE
    !        IF(I_Y.EQ.-HALFSIZE)PRINT*,"FAKE STRAIN FIELDS ARE USED"
    !        DO  I_X=-HALFSIZE,HALFSIZE
    !            CALL STRAIN_DAMPER(I_X*LATTICE_PARAMETER,I_Y*LATTICE_PARAMETER,E11,E22)
    !            !IF(E11.lt.-1d30)E11=-1d30
    !            !IF(E22.lt.-1d30)E22=-1d30
    !            STRAIN_XX(I_X,I_Y)=E11
    !            STRAIN_YY(I_X,I_Y)=E22
    !         !   print'((1x,E11.4)$)',E22
    !        ENDDO
    !        !print*,""
    !    ENDDO
    !    !stop
endsubroutine   obtain____matrices_e

subroutine      calculate_matrices_d
    use phys_constants_mod
    use diffusion_mod
    use deformations_mod
    integer i_x,i_y
    real(8) d1,d2
    print*, "diffusion coefficients marices are obtained"
    print*, "multiplied by timestep, divided by omega and h_x and h_y, implying h_x=h_y"
    do  i_y=-halfsize,halfsize
        do  i_x=-halfsize,halfsize
            d1_s(i_x,i_y)=d1(strain_xx(i_x,i_y),strain_yy(i_x,i_y),0d0)
            d2_s(i_x,i_y)=d2(strain_xx(i_x,i_y),strain_yy(i_x,i_y),0d0)
        enddo
    enddo
    d1_s=d1_s*time_step_duration!/(volume_per_atom*h_x*h_x)
    d2_s=d2_s*time_step_duration!/(volume_per_atom*h_x*h_x)
    !volume mer atom division deprecated as AVN said
endsubroutine   calculate_matrices_d

subroutine      calculate_concentration_d_product
    use diffusion_mod
    use concentration_mod
    integer i,j
    !!$OMP PARALLEL
    !!$OMP WORKSHARE
    !!$OMP DO
    do j=-halfsize,halfsize
        do i=-halfsize,halfsize
            !IF(I*I .GT. 10) CYCLE ; IF(J*J .GT. 10) CYCLE
            !print*,c_iia(i,j),d1_s(i,j),d2_s(i,j),i,j
            cdx(i,j)=c_iia(i,j)*d1_s(i,j)
            cdy(i,j)=c_iia(i,j)*d2_s(i,j)
        enddo
    enddo
    !print*,"*"
    !STOP "asddsa"
    !!$OMP END DO
    !!$OMP END WORKSHARE
    !!$OMP END PARALLEL
endsubroutine   calculate_concentration_d_product

subroutine      calculate_concentration__addition(i_time_local)
    use concentration_mod
    integer i,j
    integer(8),intent(in)::i_time_local
    !!$OMP PARALLEL
    !!$OMP WORKSHARE
    !!$OMP DO
    do j=-halfsize+2,halfsize-2
        do i=-halfsize+2,halfsize-2
!            if(j*j .gt. 325) then
!              c_add(i,j)=cdy(i,j-1)-cdy(i,j)*4d0+cdy(i,j+1)  +cdx(i-1,j)        +cdx(i+1,j)
!              cycle
!            endif
!            if(i*i .gt. 325) then
!              c_add(i,j)=cdy(i,j-1)-cdy(i,j)*4d0+cdy(i,j+1)  +cdx(i-1,j)        +cdx(i+1,j)
!              cycle
!            endif

!          if( (abs(c_add(i,j-1))+abs(c_add(i,j+1)).lt.1e-30).and. &
!              (i.gt.0) ) exit
            !IF(I*I .GT. 10) CYCLE ; IF(J*J .GT. 10) CYCLE
            !PRINT'(SP,1x,E12.3e2,$)',C_ADD(I,J)
!            PRINT*,CDY(I,J-1)
!            PRINT*,CDY(I,J)
!            PRINT*,CDY(I,J+1)
!            PRINT*,CDX(I-1,J)
!            PRINT*,CDX(I,J)
!            PRINT*,CDX(I+1,J)
            c_add(i,j)=cdy(i,j-1)-cdy(i,j)*2d0+cdy(i,j+1) +cdx(i-1,j)-cdx(i,j)*2d0+cdx(i+1,j)
!            c_add(i,j)=(&
!            cdy(i,j-1)+cdy(i,j-2)-cdy(i,j)*8d0+cdy(i,j+1)+cdy(i,j+2)+&
!            cdx(i-1,j)+cdx(i-2,j)             +cdx(i+1,j)+cdx(i+2,j))*5d-1
!more points, full formula
!            c_add(i,j)=(&
!            cdy(i,j-1) -cdy(i,j)*2d0 +cdy(i,j+1)+&
!            cdy(i,j-2)*25d-2 -cdy(i,j)*2d0*25d-2 +cdy(i,j+2)*25d-2+&
!            cdy(i-1,j) -cdy(i,j)*2d0 +cdy(i+1,j)+&
!            cdy(i-2,j)*25d-2 -cdy(i,j)*2d0*25d-2 +cdy(i+2,j)*25d-2+&
!            )


!more points, shortened formula
!              c_add(i,j)=(&
!              cdy(i,j-1) -cdy(i,j)*5d0 +cdy(i,j+1)  +&
!              cdy(i-1,j)               +cdy(i+1,j)+ (&
!              cdy(i,j-2)   +cdy(i,j+2)              +&
!              cdy(i-2,j)   +cdy(i+2,j)  )*25d-2      &
!              )*5d-1
!wiki
!              c_add(i,j)=(&
!              -cdy(i,j-2)+16d0*cdy(i,j-1) -60d0*cdy(i,j) +16d0*cdy(i,j+1)-cdy(i,j+2) +&
!              -cdy(i-2,j)+16d0*cdy(i-1,j)                +16d0*cdy(i+1,j)-cdy(i+2,j)  &
!              )/12d0 !10065s without shortings
!            if(i_time_local .lt. 8000) then
!            else
!              c_add(i,j)=cdy(i,j-1)-cdy(i,j)*4d0+cdy(i,j+1)  +cdx(i-1,j)        +cdx(i+1,j)
!            endif

        enddo
    enddo
    !!$OMP END DO
    !!$OMP END WORKSHARE
    !!$OMP END PARALLEL
endsubroutine   calculate_concentration__addition

subroutine      apply_concentration_addition
    use concentration_mod
    c_iia(-halfsize+1:halfsize-1,-halfsize+1:halfsize-1)=&
    c_iia(-halfsize+1:halfsize-1,-halfsize+1:halfsize-1)+&
    c_add(-halfsize+1:halfsize-1,-halfsize+1:halfsize-1)
endsubroutine   apply_concentration_addition

subroutine      initiate__matrices_c
    use concentration_mod
    use phys_constants_mod
    c_iia=carbon_concentration
    c_add=0d0
    cdx=0d0
    cdy=0d0
endsubroutine   initiate__matrices_c

subroutine      wo_matrices(current_time)
    use concentration_mod
    use phys_constants_mod
    USE DIFFUSION_MOD
    USE DEFORMATIONS_MOD
    integer i_x,i_y
    integer(8),intent(in) :: current_time
    character (LEN = 17) file_name

    wo_cnt=wo_cnt+1
   !write( filename_trajectory,'(A,A,A,I5.5,A)'),"test_",unique_random_string,"_grid_frame",wo_cnt,".txt"
    write( file_name,'(A,I5.5,A)'),"c_grid_",wo_cnt,"_.txt"
    open (302, file = file_name)
    write(302,'(A,I10.10,A,A,I4.4,A)')"# ",nint(current_time*time_step_duration*1d12,8)," ps;"&
    ," temperature is ",nint(temperature)," K"
    !write(302,'(A,I10.4,A)')"# temperature is ",nint(temperature)," K"
    !open(305, file = 'd_test.txt')
    !open(306, file = 'e_test.txt')
    !open(307, file = 'log_test.txt')
    do  i_x= halfsize,-halfsize,-1
      do  i_y= halfsize,-halfsize,-1
          write(302,411) &
            i_x*5d-1*lattice_parameter,&
            i_y*5d-1*lattice_parameter,&
            c_iia(i_x,i_y)/carbon_concentration
      enddo
      write(302,*) " "
    enddo
    411 format(SP,3(1x,E11.4e2))
!401 format(300(1x,E14.5e3))
    close(302)
    !close(305)
    !close(307)
    !print*," "
    print'(A,$)',"wo ok;"
endsubroutine   wo_matrices

subroutine      wo_matrices_e
    use concentration_mod
    USE DEFORMATIONS_MOD
    integer i_y
    open(303, file = 'e_xx.txt')
    open(304, file = 'e_yy.txt')
    DO  I_Y=-HALFSIZE,HALFSIZE
    !DO  I_Y= HALFSIZE,-HALFSIZE,-1
        write(303,401) ,strain_xx(-halfsize:halfsize,i_y)
        write(304,401) ,strain_yy(-halfsize:halfsize,i_y)
        !write(302,401) ,d1_s(-halfsize:halfsize,i_y)
        !write(302,401) ,log10(d1_s(-halfsize:halfsize,i_y) )
        !WRITE(302,401) ,LOG10(1D-60+ABS(C_ADD(-HALFSIZE:HALFSIZE,I_Y)))!*&
        !C_ADD(-HALFSIZE:HALFSIZE,I_Y)/(ABS(C_ADD(-HALFSIZE:HALFSIZE,I_Y)))
    ENDDO
    401 format(300(1x,E16.5e4))
    close(303)
    close(304)
    print*," "
    print*,"wo epsilon ok"
endsubroutine   wo_matrices_e

subroutine      po_info
    use phys_constants_mod
    use comp_parameter_mod
    implicit none
    open( 311,file='info.txt')
    write(311,*),lattice_parameter,   " A    is lattice_parameter"
    write(311,*),temperature,         " K    is temperature"
    write(311,*),barrier             ," eV   is carbon jump barrier"
    write(311,*),diffusion_const     ," A2/s is D0"
    write(311,*),1d-16*diffusion_const*exp(-barrier/temperature/boltzmann_const),"cm2/s is D"
    write(311,*),carbon_concentration,"      is carbon initial concentration"
    write(311,*),time_step_duration  ," s    is duration of a single time step"
    write(311,*),time_step_duration*time_max  ," s    is total physical time"
    write(311,*),k1                  ," eV   is k1"
    write(311,*),k2                  ," eV   is k2"
    write(311,*),k3                  ," eV   is k3"
    write(311,*),h_x                 ," A    is h_x"
    write(311,*),h_y                 ," A    is h_y"
    close(311)
    print*,lattice_parameter,   " A    is lattice_parameter"
    print*,temperature,         " K    is temperature"
    print*,barrier             ," eV   is carbon jump barrier"
    print*,diffusion_const     ," A2/s is D0"
    print*,1d-16*diffusion_const*exp(-barrier/temperature/boltzmann_const),"cm2/s is D"
    print*,carbon_concentration,"      is carbon initial concentration"
    print*,time_step_duration  ," s    is duration of a single time step"
    print*,time_step_duration*time_max  ," s    is total physical time"
    print*,time_max  ,"      is total time steps"
    print*,k1                  ," eV   is k1"
    print*,k2                  ," eV   is k2"
    print*,k3                  ," eV   is k3"
    print*,h_x                 ," A    is h_x"
    print*,h_y                 ," A    is h_y"
endsubroutine   po_info

real(8) function d1(e11,e22,e33)
    use deformations_mod
    use phys_constants_mod
    real(8),intent(in) :: e11,e22,e33
    real(8) :: reci_kt, exp1, exp2
    reci_kt = 1d0/(temperature*boltzmann_const)
    exp1    = exp(-reci_kt*(k2*e22+k3*e33))
    exp2    = exp(-reci_kt*(k3*e22+k2*e33))
    d1      = diffusion_const*exp(-reci_kt*barrier)*5d-1*exp(-reci_kt*(k1*e11))*(exp1+exp2)
endfunction d1

real(8) function d2(e11,e22,e33)
    use deformations_mod
    use phys_constants_mod
    real(8),intent(in) :: e11,e22,e33
    real(8) :: reci_kt, exp1, exp2
    reci_kt = 1d0/(temperature*boltzmann_const)
    exp1    = exp(-reci_kt*(k2*e33+k3*e11))
    exp2    = exp(-reci_kt*(k3*e33+k2*e11))
    d2      = diffusion_const*exp(-reci_kt*barrier)*5d-1*exp(-reci_kt*(k1*e22))*(exp1+exp2)
endfunction d2

SUBROUTINE STRAIN_DAMPER(X,Y,E11,E22)
    REAL(8),INTENT(IN)      ::  X,Y
    REAL(8),INTENT(INOUT)   ::  E11,E22
    !HIRTH LOTHE FORMULA 3.42 PAGE 57
    E11=(-Y*( Y*Y+3*X*X)*((X*X+Y*Y+1d-3)**-2))*1d-3
    E22=( Y*(-Y*Y+  X*X)*((X*X+Y*Y+1d-3)**-2))*1d-3
ENDSUBROUTINE

