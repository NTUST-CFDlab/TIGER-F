! 12 Jan 2025 - FDS

subroutine discretisation_QUICK_centre()
  use variables
  implicit none

!!! u_star calculation (x component) !!!
!$acc data present(iDx,iDy,iDz,Dxs,Dys,Dzs,nut(:,:,istart:iend+1), &
!$acc	     u(:,:,istart-2:iend+2),     v(:,:,istart-2:iend+2),     w(:,:,istart-2:iend+2), &
!$acc	u_star(:,:,istart-2:iend+2),v_star(:,:,istart-2:iend+2),w_star(:,:,istart-2:iend+2))

  !$OMP PARALLEL
   !$OMP DO PRIVATE(i,j,u_tilde_x1,u_tilde_x2,v_tilde_x1,v_tilde_x2,w_tilde_x1,w_tilde_x2,ue,uw,un,us,uf,ub,idds,Viseff) collapse(nclps)
   !$acc parallel
   !$acc loop independent private(u_tilde_x1,u_tilde_x2,v_tilde_x1,v_tilde_x2,w_tilde_x1,w_tilde_x2,ue,uw,un,us,uf,ub,idds,Viseff) collapse(3) gang vector
   do k=istart,iend
   do j=1,ny; do i=1,nx

      u_tilde_x1 = 0.5d0 * (u(i+1,j,k)+u(i,j,k)) 
      u_tilde_x2 = 0.5d0 * (u(i-1,j,k)+u(i,j,k))
      v_tilde_x1 = 0.5d0 * (v(i,j,k)+v(i+1,j,k))
      v_tilde_x2 = 0.5d0 * (v(i,j-1,k)+v(i+1,j-1,k))
      w_tilde_x1 = 0.5d0 * (w(i+1,j,k)+w(i,j,k))
      w_tilde_x2 = 0.5d0 * (w(i,j,k-1)+w(i+1,j,k-1))

        if (u_tilde_x1 > 0) then 

        ue = 0.5d0*(u(i,j,k)     + u(i+1,j,k)  ) -0.125d0*iDx(i+1)*iDx(i+1)/Dxs(i)* &
				(  (u(i+1,j,k)   - u(i,j,k)    ) / iDx(i+1) &
				 - (u(i,j,k)     - u(i-1,j,k)  ) / iDx(i)   )
                
        else

        ue = 0.5d0*(u(i,j,k)     + u(i+1,j,k)  ) -0.125d0*iDx(i+1)*iDx(i+1)/Dxs(i+1)* &
				(  (u(i+2,j,k)   - u(i+1,j,k)  ) / iDx(i+2) &
				 - (u(i+1,j,k)   - u(i,j,k)    ) / iDx(i+1) )

        end if

        if (u_tilde_x2 > 0) then 

        uw = 0.5d0*(u(i-1,j,k)   + u(i,j,k)    ) -0.125d0*iDx(i)*iDx(i)/Dxs(i-1)* &
				(  (u(i,j,k)     - u(i-1,j,k)  ) / iDx(i)   &
				 - (u(i-1,j,k)   - u(i-2,j,k)  ) / iDx(i-1) ) 
                
        else

        uw = 0.5d0*(u(i-1,j,k)   + u(i,j,k)    ) -0.125d0*iDx(i)*iDx(i)/Dxs(i)* &
				(  (u(i+1,j,k)   - u(i,j,k)    ) / iDx(i+1) &
				 - (u(i,j,k)     - u(i-1,j,k)  ) / iDx(i)   ) 

        end if

        if (v_tilde_x1 > 0) then 

        un = 0.5d0*(u(i,j,k)     + u(i,j+1,k)  ) -0.125d0*Dys(j)*Dys(j)/iDy(j)* &
				(  (u(i,j+1,k)   - u(i,j,k)    ) / Dys(j)  &
				 - (u(i,j,k)     - u(i,j-1,k)  ) / Dys(j-1))  
                
        else

        un = 0.5d0*(u(i,j,k)     + u(i,j+1,k) )-0.125d0*Dys(j)*Dys(j)/iDy(j+1)* &
				(  (u(i,j+2,k)   - u(i,j+1,k)  ) / Dys(j+1)&
				 - (u(i,j+1,k)   - u(i,j,k)    ) / Dys(j)  )

        end if

        if (v_tilde_x2 > 0) then 

        us = 0.5d0*(u(i,j-1,k)   + u(i,j,k)    ) -0.125d0*Dys(j-1)*Dys(j-1)/iDy(j-1)* &
				(  (u(i,j,k)     - u(i,j-1,k)  ) / Dys(j-1)&
				 - (u(i,j-1,k)   - u(i,j-2,k)  ) / Dys(j-2)) 
                
        else

        us = 0.5d0*(u(i,j-1,k)   + u(i,j,k)    ) -0.125d0*Dys(j-1)*Dys(j-1)/iDy(j)* &
				(  (u(i,j+1,k)   - u(i,j,k)    ) / Dys(j)  &
				 - (u(i,j,k)     - u(i,j-1,k)  ) / Dys(j-1))  

        end if

        if (w_tilde_x1 > 0) then 

        uf = 0.5d0*(u(i,j,k)     + u(i,j,k+1)  ) -0.125d0*Dzs(k)*Dzs(k)/iDz(k)* &
				(  (u(i,j,k+1)   - u(i,j,k)    ) / Dzs(k)  &
				 - (u(i,j,k)     - u(i,j,k-1)  ) / Dzs(k-1))

            
        else

        uf = 0.5d0*(u(i,j,k)     + u(i,j,k+1)  ) -0.125d0*Dzs(k)*Dzs(k)/iDz(k+1)* &
				(  (u(i,j,k+2)   - u(i,j,k+1)  ) / Dzs(k+1)&
				 - (u(i,j,k+1)   - u(i,j,k)    ) / Dzs(k)  )


        end if

        if (w_tilde_x2 > 0) then 

        ub = 0.5d0*(u(i,j,k-1)   + u(i,j,k)    ) -0.125d0*Dzs(k-1)*Dzs(k-1)/iDz(k-1)* &
				(  (u(i,j,k)     - u(i,j,k-1)  ) / Dzs(k-1)&
				 - (u(i,j,k-1)   - u(i,j,k-2)  ) / Dzs(k-2))  
            
        else

        ub = 0.5d0*(u(i,j,k-1)   + u(i,j,k)    ) -0.125d0*Dzs(k-1)*Dzs(k-1)/iDz(k)* &
				(  (u(i,j,k+1)   - u(i,j,k)    ) / Dzs(k)  &
				 - (u(i,j,k)     - u(i,j,k-1)  ) / Dzs(k-1))

        end if

		idds = 0.5d0*iDx(i)/Dxs(i)
		Viseff = nu + LES*((1.d0 - idds)*nut(i,j,k) + idds*nut(i+1,j,k))	! (nu + nut) on east face

        u_star(i,j,k) = dt*(										&
						- (u_tilde_x1*ue-u_tilde_x2*uw) / Dxs(i) 	&
						- (v_tilde_x1*un-v_tilde_x2*us) / iDy(j) 	&
						- (w_tilde_x1*uf-w_tilde_x2*ub) / iDz(k) 	&

                        + Viseff*( (u(i+1,j,k)-u(i,j,k)) / iDx(i+1) - (u(i,j,k)-u(i-1,j,k)) / iDx(i) ) / Dxs(i) &
                        + Viseff*( (u(i,j+1,k)-u(i,j,k)) / Dys(j) - (u(i,j,k)-u(i,j-1,k)) / Dys(j-1) ) / iDy(j) &
                        + Viseff*( (u(i,j,k+1)-u(i,j,k)) / Dzs(k) - (u(i,j,k)-u(i,j,k-1)) / Dzs(k-1) ) / iDz(k) )
                        

  	end do;end do
  	end do
  	!$OMP END DO


!!! v_star calculation (y component) !!!

   !$OMP DO PRIVATE(i,j,u_tilde_y1,u_tilde_y2,v_tilde_y1,v_tilde_y2,w_tilde_y1,w_tilde_y2,ve,vw,vn,vs,vf,vb,idds,Viseff) collapse(nclps)
   !$acc loop independent private(u_tilde_y1,u_tilde_y2,v_tilde_y1,v_tilde_y2,w_tilde_y1,w_tilde_y2,ve,vw,vn,vs,vf,vb,idds,Viseff) collapse(3) gang vector
   do k=istart,iend
   do j=1,ny;do i=1,nx

      u_tilde_y1 = 0.5d0 * (u(i,j,k)+u(i,j+1,k))
      u_tilde_y2 = 0.5d0 * (u(i-1,j,k)+u(i-1,j+1,k))
      v_tilde_y1 = 0.5d0 * (v(i,j,k)+v(i,j+1,k))
      v_tilde_y2 = 0.5d0 * (v(i,j-1,k)+v(i,j,k))
      w_tilde_y1 = 0.5d0 * (w(i,j,k)+w(i,j+1,k))
      w_tilde_y2 = 0.5d0 * (w(i,j,k-1)+w(i,j+1,k-1))

        if (u_tilde_y1 > 0) then 
        
        ve = 0.5d0*(v(i,j,k)     + v(i+1,j,k)  ) -0.125d0*Dxs(i)*Dxs(i)/iDx(i)* &
				(  (v(i+1,j,k)   - v(i,j,k)    ) / Dxs(i)  &
				 - (v(i,j,k)     - v(i-1,j,k)  ) / Dxs(i-1))

        else

        ve = 0.5d0*(v(i,j,k)     + v(i+1,j,k)  ) -0.125d0*Dxs(i)*Dxs(i)/iDx(i+1)* &
				(  (v(i+2,j,k)   - v(i+1,j,k)  ) / Dxs(i+1)&
				 - (v(i+1,j,k)   - v(i,j,k)    ) / Dxs(i)  )  

        end if

        if (u_tilde_y2 > 0) then 

        vw = 0.5d0*(v(i-1,j,k)   + v(i,j,k)    ) -0.125d0*Dxs(i-1)*Dxs(i-1)/iDx(i-1)* &
				(  (v(i,j,k)     - v(i-1,j,k)  ) / Dxs(i-1)&
				 - (v(i-1,j,k)   - v(i-2,j,k)  ) / Dxs(i-2)) 

        else 

        vw = 0.5d0*(v(i-1,j,k)   + v(i,j,k)    ) -0.125d0*Dxs(i-1)*Dxs(i-1)/iDx(i)* &
				(  (v(i+1,j,k)   - v(i,j,k)    ) / Dxs(i)  &
				 - (v(i,j,k)     - v(i-1,j,k)  ) / Dxs(i-1))  

        end if

        if (v_tilde_y1 > 0) then 

        vn = 0.5d0*(v(i,j,k)     + v(i,j+1,k)  ) -0.125d0*iDy(j+1)*iDy(j+1)/Dys(j)* &
				(  (v(i,j+1,k)   - v(i,j,k)    ) / iDy(j+1) &
				 - (v(i,j,k)     - v(i,j-1,k)  ) / iDy(j)   ) 
            
        else

        vn = 0.5d0*(v(i,j,k)     + v(i,j+1,k)  ) -0.125d0*iDy(j+1)*iDy(j+1)/Dys(j+1)* &
				(  (v(i,j+2,k)   - v(i,j+1,k)  ) / iDy(j+2) &
				 - (v(i,j+1,k)   - v(i,j,k)    ) / iDy(j+1) ) 

        end if

        if (v_tilde_y2 > 0) then 

        vs = 0.5d0*(v(i,j-1,k)   + v(i,j,k)    ) -0.125d0*iDy(j)*iDy(j)/Dys(j-1)* &
				(  (v(i,j,k)     - v(i,j-1,k)  ) / iDy(j)   &
				 - (v(i,j-1,k)   - v(i,j-2,k)  ) / iDy(j-1) ) 
            
        else

		vs = 0.5d0*(v(i,j-1,k)   + v(i,j,k)    ) -0.125d0*iDy(j)*iDy(j)/Dys(j)* &
              	(  (v(i,j+1,k)   - v(i,j,k)    ) / iDy(j+1)&
                 - (v(i,j,k)     - v(i,j-1,k)  ) / iDy(j)  ) 

        end if


        if (w_tilde_y1 > 0) then 

		vf = 0.5d0*(v(i,j,k)     + v(i,j,k+1)  ) -0.125d0*Dzs(k)*Dzs(k)/iDz(k)* &
				(  (v(i,j,k+1)   - v(i,j,k)    ) / Dzs(k)  &
				 - (v(i,j,k)     - v(i,j,k-1)  ) / Dzs(k-1))
            
        else

		vf = 0.5d0*(v(i,j,k)     + v(i,j,k+1)  ) -0.125d0*Dzs(k)*Dzs(k)/iDz(k+1)* &
				(  (v(i,j,k+2)   - v(i,j,k+1)  ) / Dzs(k+1)&
				 - (v(i,j,k+1)   - v(i,j,k)    ) / Dzs(k)  )
        end if

        if (w_tilde_y2 > 0) then 

		vb = 0.5d0*(v(i,j,k-1)   + v(i,j,k)    ) -0.125d0*Dzs(k-1)*Dzs(k-1)/iDz(k-1)* &
				(  (v(i,j,k)     - v(i,j,k-1)  ) / Dzs(k-1)&
				 - (v(i,j,k-1)   - v(i,j,k-2)  ) / Dzs(k-2))
            
        else

		vb = 0.5d0*(v(i,j,k-1)   + v(i,j,k)    ) -0.125d0*Dzs(k-1)*Dzs(k-1)/iDz(k)* &
				(  (v(i,j,k+1)   - v(i,j,k)    ) / Dzs(k)  &
				 - (v(i,j,k)     - v(i,j,k-1)  ) / Dzs(k-1))

        end if


		idds = 0.5d0*iDy(j)/Dys(j)
		Viseff = nu + LES*((1.d0 - idds)*nut(i,j,k) + idds*nut(i,j+1,k))	! (nu + nut) on north face

        v_star(i,j,k) = dt*(										&
						- (u_tilde_y1*ve-u_tilde_y2*vw) / iDx(i)	&
                        - (v_tilde_y1*vn-v_tilde_y2*vs) / Dys(j)	&
                        - (w_tilde_y1*vf-w_tilde_y2*vb) / iDz(k)	&

                        + Viseff*( (v(i+1,j,k)-v(i,j,k)) / Dxs(i) - (v(i,j,k) - v(i-1,j,k)) / Dxs(i-1) ) / iDx(i) &
                        + Viseff*( (v(i,j+1,k)-v(i,j,k)) / iDy(j+1) - (v(i,j,k) - v(i,j-1,k)) / iDy(j) ) / Dys(j) &
                        + Viseff*( (v(i,j,k+1)-v(i,j,k)) / Dzs(k) - (v(i,j,k) - v(i,j,k-1)) / Dzs(k-1) ) / iDz(k) )
                        !+dt*0.5d0*( F_tavey(i,j,k)+F_tavey(i,j+1,k) )*unDBD_onoff
                        !-dt*( Viseff*ABS(F_tavey(i,j,k)) )
    
  end do;end do
  end do
  !$OMP END DO  
  

!!! w_star calculation (z component) !!!

   !$OMP DO PRIVATE(i,j,u_tilde_z1,u_tilde_z2,v_tilde_z1,v_tilde_z2,w_tilde_z1,w_tilde_z2,we,ww,wn,ws,wf,wb,idds,Viseff) collapse(nclps)
   !$acc loop independent private(u_tilde_z1,u_tilde_z2,v_tilde_z1,v_tilde_z2,w_tilde_z1,w_tilde_z2,we,ww,wn,ws,wf,wb,idds,Viseff) collapse(3) gang vector
   do k=istart,iend 
   do j=1,ny;do i=1,nx

        u_tilde_z1 = 0.5d0 * (u(i,j,k+1)+u(i,j,k)) 
        u_tilde_z2 = 0.5d0 * (u(i-1,j,k+1)+u(i-1,j,k))  
        v_tilde_z1 = 0.5d0 * (v(i,j,k)+v(i,j,k+1))
        v_tilde_z2 = 0.5d0 * (v(i,j-1,k+1)+v(i,j-1,k))
        w_tilde_z1 = 0.5d0 * (w(i,j,k+1)+w(i,j,k))
        w_tilde_z2 = 0.5d0 * (w(i,j,k)+w(i,j,k-1))

        if (u_tilde_z1 > 0) then 
    
		we = 0.5d0*(w(i,j,k)     + w(i+1,j,k)  ) -0.125d0*Dxs(i)*Dxs(i)/iDx(i)* &
				(  (w(i+1,j,k)   - w(i,j,k)    ) / Dxs(i)  &
				 - (w(i,j,k)     - w(i-1,j,k)  ) / Dxs(i-1)) 

        else

		we = 0.5d0*(w(i,j,k)     + w(i+1,j,k)  ) -0.125d0*Dxs(i)*Dxs(i)/iDx(i+1)* &
				(  (w(i+2,j,k)   - w(i+1,j,k)  ) / Dxs(i+1) &
				 - (w(i+1,j,k)   - w(i,j,k)    ) / Dxs(i)   )

        end if

        if (u_tilde_z2 > 0) then 

		ww = 0.5d0*(w(i-1,j,k)   + w(i,j,k)    ) -0.125d0*Dxs(i-1)*Dxs(i-1)/iDx(i-1)* &
				(  (w(i,j,k)     - w(i-1,j,k)  ) / Dxs(i-1) &
				 - (w(i-1,j,k)   - w(i-2,j,k)  ) / Dxs(i-2) )

        else 

		ww = 0.5d0*(w(i-1,j,k)   + w(i,j,k)    ) -0.125d0*Dxs(i-1)*Dxs(i-1)/iDx(i)* &
				(  (w(i+1,j,k)   - w(i,j,k)    ) / Dxs(i)  &
				 - (W(i,j,k)     - w(i-1,j,k)  ) / Dxs(i-1)) 

        end if

        if (v_tilde_z1 > 0) then 

		wn = 0.5d0*(w(i,j,k)     + w(i,j+1,k)  ) -0.125d0*Dys(j)*Dys(j)/iDy(j)* &
				(  (w(i,j+1,k)   - w(i,j,k)    ) / Dys(j)  &
                 - (w(i,j,k)     - w(i,j-1,k)  ) / Dys(j-1))
            
        else

		wn = 0.5d0*(w(i,j,k)     + w(i,j+1,k)  ) -0.125d0*Dys(j)*Dys(j)/iDy(j+1)* &
				(  (w(i,j+2,k)   - w(i,j+1,k)  ) / Dys(j+1)&
                 - (w(i,j+1,k)   - w(i,j,k)    ) / Dys(j)  ) 

        end if

        if (v_tilde_z2 > 0) then 

		ws = 0.5d0*(w(i,j-1,k)   + w(i,j,k)    ) -0.125d0*Dys(j-1)*Dys(j-1)/iDy(j-1)* &
				(  (w(i,j,k)     - w(i,j-1,k)  ) / Dys(j-1)&
                 - (w(i,j-1,k)   - w(i,j-2,k)  ) / Dys(j-2))
            
        else

		ws = 0.5d0*(w(i,j-1,k)   + w(i,j,k)    ) -0.125d0*Dys(j-1)*Dys(j-1)/iDy(j)* &
				(  (w(i,j+1,k)   - w(i,j,k)    ) / Dys(j)  &
                 - (w(i,j,k)     - w(i,j-1,k)  ) / Dys(j-1))

        end if


        if (w_tilde_z1 > 0) then 

		wf = 0.5d0*(w(i,j,k)     + w(i,j,k+1)  ) -0.125d0*iDz(k+1)*iDz(k+1)/Dzs(k)* &
				(  (w(i,j,k+1)   - w(i,j,k)    ) / iDz(k+1)&
                 - (w(i,j,k)     - w(i,j,k-1)  ) / iDz(k)  )
            
        else

		wf = 0.5d0*(w(i,j,k)     + w(i,j,k+1)  ) -0.125d0*iDz(k+1)*iDz(k+1)/Dzs(k+1)* &
				(  (w(i,j,k+2)   - w(i,j,k+1)  ) / iDz(k+2)&
                 - (w(i,j,k+1)   - w(i,j,k)    ) / iDz(k+1))
        end if

        if (w_tilde_z2 > 0) then 

		wb = 0.5d0*(w(i,j,k-1)   + w(i,j,k)    ) -0.125d0*iDz(k)*iDz(k)/Dzs(k-1)* &
				(  (w(i,j,k)     - w(i,j,k-1)  ) / iDz(k)  &
                 - (w(i,j,k-1)   - w(i,j,k-2)  ) / iDz(k-1))
            
        else

		wb = 0.5d0*(w(i,j,k-1)   + w(i,j,k)    ) -0.125d0*iDz(k)*iDz(k)/Dzs(k)* &
                (  (w(i,j,k+1)   - w(i,j,k)    ) / iDz(k+1)&
				 - (w(i,j,k)     - w(i,j,k-1)  ) / iDz(k)  )

        end if

		idds = 0.5d0*iDz(k)/Dzs(k)
		Viseff = nu + LES*((1.d0 - idds)*nut(i,j,k) + idds*nut(i,j,k+1))	! (nu + nut) on front face

        w_star(i,j,k) = dt*(										&
						- (u_tilde_z1*we-u_tilde_z2*ww) / iDx(i)	&
                        - (v_tilde_z1*wn-v_tilde_z2*ws) / iDy(j)	&
                        - (w_tilde_z1*wf-w_tilde_z2*wb) / Dzs(k)	&
                                
                        + Viseff*( ( w(i+1,j,k)-w(i,j,k) ) / Dxs(i) -  ( w(i,j,k)-w(i-1,j,k) ) / Dxs(i-1)  ) / iDx(i) &
                        + Viseff*( ( w(i,j+1,k)-w(i,j,k) ) / Dys(j) -  ( w(i,j,k)-w(i,j-1,k) ) / Dys(j-1)  ) / iDy(j) &
                        + Viseff*( ( w(i,j,k+1)-w(i,j,k) ) / iDz(k+1) - ( w(i,j,k)-w(i,j,k-1) ) / iDz(k)   ) / Dzs(k) )
                        !+ dt*0.5d0*( F_tavex(i,j,k)+F_tavex(i,j,k+1) )*unDBD_onoff
                        !+ dt*( Viseff*ABS(F_tavex(i,j,k)) )
  end do;end do
  end do
  !$acc end parallel
  !$OMP END DO
 !$OMP END PARALLEL
!$acc end data

end subroutine discretisation_QUICK_centre
