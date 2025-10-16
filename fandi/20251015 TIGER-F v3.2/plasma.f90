subroutine plasma()
use variables
implicit none

real*8 :: distance

!           Shyy et. al (2002)
!           y
            ! \
            !  
            !    a
            !  
!===========! /
            !------------------------------------------> x
!           \      =================/
!                        b

    !initialize zero point of Plasma
    PlasmaZc = az(1)
    PlasmaYc = ay(1)

    PlasmaZr = az(198)      !1.525522743313484E-002
    PlasmaYr = ay(198)      !2.080064932467730E-002   

    PlasmaZu = PlasmaZc - COS((-AOA)*PI/180.d0)*len_a
    PlasmaYu = PlasmaYc - SIN((-AOA)*PI/180.d0)*len_a


    !--- Tanssfer --------------------------------------------------------------------------
    PlasmaZc = (PlasmaZc-PlasmaZu)/0.02d0*0.03535d0 + PlasmaZu
    PlasmaYc = (PlasmaYc-PlasmaYu)/0.02d0*0.03535d0 + PlasmaYu
    !--- Tanssfer --------------------------------------------------------------------------
    

    !write(*,*) PlasmaZu, PlasmaYu 
    !write(*,*) PlasmaZc, PlasmaYc 


    PlasmaVy1 = PlasmaYr-PlasmaYc
    PlasmaVy2 = PlasmaYu-PlasmaYr
    PlasmaVy3 = PlasmaYc-PlasmaYu
    PlasmaVz1 = PlasmaZr-PlasmaZc
    PlasmaVz2 = PlasmaZu-PlasmaZr
    PlasmaVz3 = PlasmaZc-PlasmaZu



    !initialize parameter of Plasma
    do k=1,nz; do j=1,ny; do i=1,nx
        EE(i,j,k) = 1.D0
        edelta(i,j,k) = 0.D0
        F_tavex(i,j,k) = 0.D0
        F_tavey(i,j,k) = 0.D0
    end do; end do; end do

    !----------------------- Dimensional parameter calculation -----------------------!

    E0 = poto / ( len_d * reference_L )
    k1 = ( E0 - Eb ) / ( len_b * reference_L )
    k2 = ( E0 - Eb ) / ( len_a * reference_L )
    E0 = E0 + k2*0.01535d0*reference_L

    k1 = ( E0 - Eb ) / ( len_b/0.02d0*0.03535d0 * reference_L )
    k2 = ( E0 - Eb ) / ( len_a/0.02d0*0.03535d0 * reference_L )




    !$OMP PARALLEL DO PRIVATE(i,c1,c2,c3)
    do k=1,nz; do j=1,ny

        i=0

        c1 = PlasmaVz1*( PlasmaYr - Ys(j) ) - PlasmaVy1*( PlasmaZr - Zs(k) )
        c2 = PlasmaVz2*( PlasmaYu - Ys(j) ) - PlasmaVy2*( PlasmaZu - Zs(k) )
        c3 = PlasmaVz3*( PlasmaYc - Ys(j) ) - PlasmaVy3*( PlasmaZc - Zs(k) )


        !----------------------- Dimensional parameter calculation -----------------------!
        if ( c1 <= 0 .AND. c2 <= 0 .AND. c3 <= 0 ) then

            EE(i,j,k) = E0 - k1*ABS(SIN((AOA)*PI/180.d0)*(Zs(k)-PlasmaZc) + COS((AOA)*PI/180.d0)*(Ys(j)-PlasmaYc))*reference_L &
                           - k2*ABS(COS((AOA)*PI/180.d0)*(Zs(k)-PlasmaZc) - SIN((AOA)*PI/180.d0)*(Ys(j)-PlasmaYc))*reference_L

            edelta(i,j,k) = 1.d0

            F_tavex(i,j,k) = theta * alfa * roc * ec * edelta(i,j,k) * delta_t * ( EE(i,j,k)*k2 ) / (( k1**2.d0 + k2**2.d0 ) ** 0.5d0)
            F_tavey(i,j,k) = theta * alfa * roc * ec * edelta(i,j,k) * delta_t * ( EE(i,j,k)*k1 ) / (( k1**2.d0 + k2**2.d0 ) ** 0.5d0)

            F_tavex(i,j,k) = F_tavex(i,j,k) / (Freestream**2.d0) / 1.225d0 * reference_L * SQRT(1.d0-ETA(i,j,k))
            F_tavey(i,j,k) = F_tavey(i,j,k) / (Freestream**2.d0) / 1.225d0 * reference_L * SQRT(1.d0-ETA(i,j,k))

        else 

            edelta(i,j,k) = 0.d0
            F_tavex(i,j,k) = 0.d0
            F_tavey(i,j,k) = 0.d0

        endif


        do i=1,nx+1

            edelta(i,j,k) = edelta(0,j,k)     
            F_tavex(i,j,k) = F_tavex(0,j,k)     
            F_tavey(i,j,k) = F_tavey(0,j,k)   

        end do


    end do; end do
    !$OMP END PARALLEL DO

    unDBD_onoff=1

    

end subroutine plasma


subroutine unsteady_plasma()
use variables
implicit none

    unDBD_cycle = 1.d0/(unDBD+1.E-10)

    if ( MOD(time,unDBD_cycle) <= 0.1d0*unDBD_cycle ) then
        unDBD_onoff=1
    else
        unDBD_onoff=0
    endif

end subroutine