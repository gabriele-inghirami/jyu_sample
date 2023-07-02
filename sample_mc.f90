! *  *  *  *  *
! *
! *  JYU_SAMPLE
! *
! *  Version: 0.5
! *                          
! *  Date (DD/MM/YYYY): 08/10/2019
! *                  
! *  File: jyu_sample/sample_mc.f90 
! *                           
! *  Author: 
! *  Gabriele Inghirami (University of Jyvaskyla and Helsinki Institute of physics- Finland)
! *  E-mail: gabriele.g.inghirami@jyu.fi
! *  in collaboration with:
! *  Harri Niemi (University of Jyvaskyla and Helsinki Institute of physics- Finland)
! * 
! *  Copyright - Important attribution note:
! *               
! *  THIS PROGRAM CONTAINS CODE:
! *
! *  - INCLUDED IN ECHO-QGP v.1.0.x AND DEVELOPED BY:
! *  Valentina Rolando (INFN and University of Ferrara - Italy)             
! *  with the contribution of:
! *  Giuseppe Pagliara and Alessandro Drago (INFN and University of Ferrara - Italy)
! *  References: 
! *  Eur.Phys.J. C73 (2013) 2524 - arXiv: 1305.7052
! *  Eur.Phys.J. C75 (2015) no.9, 406, Erratum: Eur.Phys.J. C78 (2018) no.5, 354 - 1501.04468
! *                         
! *  License: GPL version 2.0 (Please, read the file LICENSE.TXT)       
! *                                      
! *  This program is free software; you can redistribute it and/or  
! *  modify it under the terms of the GNU General Public License   
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.       
! *                                                                
! *  This program is distributed in the hope that it will be useful, 
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of 
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! *  GNU General Public License for more details.                 
! *                                                                   
! *  You should have received a copy of the GNU General Public License 
! *  along with this program; if not, write to the Free Software     
! *  Foundation, Inc., 51 Franklin Street, Fifth Floor,             
! *  Boston, MA  02110-1301, USA.                                  
! *
! *  *  *  *  * 


	 module work 
	  
	  use constants
      use settings
	  implicit none
!  procedure(), pointer:: cooper_frye

	  real, dimension(0:3) :: xfo, ufo, dV
	  real, dimension(1:3) :: vfo
	  real rho, prex
      real, dimension(1:N_dth) :: dth_arr
      real, dimension(1:N_dphi-1) :: dphi_arr
      real, dimension(1:N_dp) :: dp_arr
      real, dimension(1:4, 1:N_dth, 1:N_dphi-1,1:N_dp) :: dmom_arr
      real, dimension(0:3,0:3) :: Lambda_boost, Lambda_boost_back, pi_shear, pi_shear_lrf, tmparr
      logical :: viscosity
      integer, parameter :: kt=0, kx=1, ky=2, kz=3
      real(8), allocatable, dimension (:) :: particle_density, particle_density_common
      real(8), dimension(0:3) :: part_coord, part_mom
	  contains
	  ! ! !     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

          subroutine fill_N_arrays()
          implicit none

          integer i,j,k
          real :: dx

          dx=PMAX_BOX/(N_dp-1)
          do i=1,N_dp
             dp_arr(i)=dx*(i-1)
          end do

          dx=2.d0*GREEK_PI/(N_dphi-1)
          do i=1,N_dphi-1
             dphi_arr(i)=dx*(i-1)
          end do
         
          dx=GREEK_PI/(N_dth-1)
          do i=1,N_dth
             dth_arr(i)=dx*(i-1)
          end do

          do k=1,N_dp
             do j=1,N_dphi-1
                do i=1,N_dth
                   dmom_arr(1,i,j,k)=dp_arr(k)*dp_arr(k)
                   dmom_arr(2,i,j,k)=dp_arr(k)*sin(dth_arr(i))*cos(dphi_arr(j))
                   dmom_arr(3,i,j,k)=dp_arr(k)*sin(dth_arr(i))*sin(dphi_arr(j))
                   dmom_arr(4,i,j,k)=dp_arr(k)*cos(dth_arr(i))
                 end do
              end do
           end do
         
          write(*,*) "Arrays for coarse loop filled"

          end subroutine fill_N_Arrays
	 subroutine work_choose_random_index(iii, MMM)
	! ! !  extracts integer random value between 0 and MMM
	    implicit none
	    integer, intent(out) :: iii !index
	    integer, intent(in)  :: MMM !max 
	    
	    real(8) :: r
	    
	    call random_number(r)
	    r=(MMM*1.0+1.0)*r
	    iii=int(r)
	    return
	 end subroutine work_choose_random_index
	 ! ! !     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	!  subroutine work_particles_generator(pa_indx, hy_indx, hyper_i, stt)
	  subroutine work_global_density(hyper_i, particle_density, global_density)
	!  this routine calculates N= Sum_i (N_i) in a particular cell
	!  where N_i = n_i u^mu dSigma_mu (fm^-3 * fm^3 = adimensional )
	  use constants
	  use common
      use eos
	  use numlib
      use settings
	  implicit none
	  
	 !   contains desity for each specie  
	!   hypersurface array for the ideal case    
	    real(8), intent(in), dimension(1:27) :: hyper_i 
            real(8), allocatable, dimension(:) :: particle_density

	!   this is the sum over the species of the particle density 
	!   N_i = sum_i (n_i u^mu dSigma_mu)
	    real(8), intent(out) :: global_density
	    
	!   u^mu d_Sigma_mu
	    real(8) u_dot_dsigma
	    ! counters
	    integer ipart, ipi
	    
	    real(8) gm2T, tausq
	    real(8) bes_arg, bessum, m, gspin, edens, mu, Tfo
	    real(8),save :: global_density_common=0.d0
!            logical, save :: first_time=.true.
	      
	    ! tau  or t 
	    xfo(0)=hyper_i(1)
	    ! x
	    xfo(1)=hyper_i(2)
	    ! y
	    xfo(2)=hyper_i(3)
	    ! eta or z
	    xfo(3)=hyper_i(4)
	    
	    tausq=xfo(0)*xfo(0)
	!     Volume element
	    dV(0)=hyper_i(5)*xfo(0)
	    dV(1)=hyper_i(6)*xfo(0)
	    dV(2)=hyper_i(7)*xfo(0)
	    dV(3)=hyper_i(8)*xfo(0)
	   
            Tfo=hyper_i(9)
            edens=hyper_i(10)

             
          
	!     velocities
	    vfo(1)=hyper_i(12)
	    vfo(2)=hyper_i(13)
	    vfo(3)=0    
	    
	  
	    !global_density=0.d0
	    !particle_density=0.d0
	    
	    ! We compute the scalar product  u^mu dSigma_mu 
	    ! in order to reject the negative contributions
	    call scalar_product(xfo, vfo, dV, u_dot_dsigma) !! u^mu dSigma_mu (x is needed for the metrics)
	    
	    ! check if this is a negative contribution and possibly reject it
	    ! if negative: global_density=0.0 -> no particle is produced 

!            if(first_time) then
	      particle_density_common=0.d0
	      global_density_common=0.d0
              !we start from 2, because 1 is the photon
	      do ipart=2,npart_main
                 m=pdata(ipart)%mass
                 gspin=pdata(ipart)%spin_degeneracy
                 call get_mu(ipart,edens,mu) 
                 if(m .eq. 0) then
                   cycle
                 end if
	         bes_arg=m/Tfo
		  ! if we are dealing with pions use the expansion to higher order 
		  ! (up to CUT_PI)
		  !if ((abs(pdg_number(ipart)) .eq. 211) .or. (pdg_number(ipart) .eq. 111) .or. ) then
		  if (m .lt. 0.6) then
		    bessum=0.0
		    do ipi=1, CUT_PI
		      bessum=bessum + (bessk2(1.0*ipi*bes_arg)*1.0/(1.0*ipi))*exp(ipi*mu/Tfo)
		    end do
		  else
		    bessum=bessk2(bes_arg)*exp(mu/Tfo)
		  endif
		  gm2T= gspin* m*m * Tfo
          !       write(*,*) "g factor: ", ipart, g(ipart)
          !		  if(mu .eq. 0) then
                  particle_density_common(ipart)=  gm2T * bessum
          !	          else
          !                 particle_density_common(ipart)=exp(mu/Tfo) * gm2T * bessum
          !		  end if
	!	  write(*,*) ipart, particle_density_common(ipart)
	!	  write(*,*) m(ipart), mu(ipart), Tfo, g(ipart), gm2T, bessum, exp(mu(ipart)/Tfo) * gm2T * bessum
		  global_density_common=global_density_common+particle_density_common(ipart)
              end do !!!npart      

!              first_time=.false.
!            end if

	    if (u_dot_dsigma .gt. 0.0) then    
	      global_density=global_density_common*NOR_DENS*u_dot_dsigma
	      particle_density(:)=particle_density_common(:)*NOR_DENS*u_dot_dsigma
            else
              global_density=0.
              particle_density(:)=0.
	    endif
	  return
	 end subroutine work_global_density 
	 ! !     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	 subroutine work_assign_features(ipart, ni, Ntot)
	  use eos, only: npart_main
	  implicit none
	  real, intent (in) :: Ntot !! Global particle density of this cell
	  real(8), intent (in), allocatable, dimension(:) :: ni ! particle density array (i for the part. specie)
	  integer, intent(out) :: ipart 
	  integer i
	  real harvest, rnd, nsum
	  
	  call random_number(harvest)      
	  rnd=harvest*Ntot ! we extract a random between 0 and the global density
          !we start from 2, because 1 is the photon
	  ipart=2      
	  nsum=ni(ipart)
	  do while (rnd > nsum ) 
	    ipart=ipart+1
	    nsum=nsum+ni(ipart)
	  end do
            

	 return
	 end subroutine work_assign_features
	 ! ! !     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  
	  subroutine scalar_product(x,v,dsigma, scprod)
	  implicit none 
	! !  computes the scalar product u^mu dSigma_mu
	  real, intent(in), dimension(0:3) :: x,dsigma
	  real, intent(in), dimension(1:3) :: v
	  real, intent(out) :: scprod
	  real betasq, tausq
	  integer i
	  real, dimension(0:3) :: u_mu
	  
	  tausq=x(0)*x(0)
	  
	  betasq=v(1)*v(1)+v(2)*v(2)+v(3)*v(3)*tausq
	    u_mu(0)=1.0/sqrt(1-betasq)
	    do i=1, 3
	      u_mu(i)=u_mu(0)*v(i)
	    end do
	  
	    scprod= u_mu(0)*dsigma(0) + u_mu(1)*dsigma(1) + &
		  & u_mu(2)*dsigma(2) + u_mu(3)*dsigma(3)
	    
	  return
	 end subroutine scalar_product
	 ! ! !     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	  subroutine work_momentum_generator(ipart, hyper_i, location, ppp, feedback)
          use eos
 implicit none
  real, intent(in), dimension(1:27) :: hyper_i  
  integer,intent(inout):: ipart
  real, dimension(0:3) ::  ppp, temp_mom, location

  integer i, j
  real(8) cftrial, CoopFrye, massimo, r, delta_kr, pmuds, pmutot
!  integer, parameter :: maxtrials=1000000000
  integer, parameter :: maxtrials=20000000
  integer counttrials
  integer feedback
  real(8) :: tau,x,y,eta,m,m2,ch,sh,p0,pt2,p,phi,theta,glf,betasq,px,py,pz,pmaxm,ds_tau,ds_x,ds_y,ds_eta,mu,bf
  real(8), dimension(1:3) :: vb
  real(8), dimension(0:3) :: um,ub,dsigma_hyp,dsigma_hyp_lrf,um_back
  real(8) :: distr_max, pguess, mom_shear_prod
  real(8) :: edens, press, mt, ptau, peta, rap, Tfo, yrap, yrap_minus_eta
  real(8) :: ptt_bj, ptx_bj, pty_bj, ptz_bj, pxt_bj, pxz_bj, pyt_bj, pyz_bj, pzt_bj, pzx_bj, pzy_bj, pzz_bj

  counttrials=0
  feedback=0

  m=pdata(ipart)%mass
  m2=m**2

  tau=hyper_i(1)
  x=hyper_i(2)
  y=hyper_i(3)
  if(D2) then
    call random_number(r)
    eta=(2.d0*r-1.d0)*eta_pseudorap_max !we add the random pseudorapidity
  else  
    eta=hyper_i(4)
  end if
  ds_tau=hyper_i(5)*tau
  ds_x=hyper_i(6)*tau
  ds_y=hyper_i(7)*tau
  ds_eta=hyper_i(8)*tau
  Tfo=hyper_i(9)
  vb(1)=hyper_i(12)
  vb(2)=hyper_i(13)
  vb(3)=0
  ch=cosh(eta)
  sh=sinh(eta)
  location(0)=tau*ch
  location(1)=x
  location(2)=y
  location(3)=tau*sh
 
  edens=hyper_i(10)
  press=hyper_i(11)
  call get_mu(ipart, edens, mu)
  bf=real(pdata(ipart)%bf)
  
  ptt_bj=hyper_i(15)
  ptx_bj=hyper_i(16)
  pty_bj=hyper_i(17)
  ptz_bj=0.
  pxt_bj=ptx_bj
  pi_shear(kx,kx)=hyper_i(18)
  pi_shear(kx,ky)=hyper_i(20)
  pxz_bj=0.
  pyt_bj=pty_bj
  pi_shear(ky,kx)=pi_shear(kx,ky)
  pi_shear(ky,ky)=hyper_i(19)
  pyz_bj=0.
  pzt_bj=ptz_bj
  pzx_bj=pxz_bj
  pzy_bj=pyz_bj
  pzz_bj=hyper_i(21)/(tau**2) !we restore the pi^{eta eta} form
  
  !first, we transform the tensor in Minkowski coordinates
  !TO RECHECK BETTER!!!!
  pi_shear(kt,kt)=ptt_bj*ch**2+pzz_bj*(tau*sh)**2
  pi_shear(kt,kx)=ptx_bj*ch
  pi_shear(kx,kt)=pi_shear(kt,kx)
  pi_shear(kt,ky)=pty_bj*ch
  pi_shear(ky,kt)=pi_shear(kt,ky)
  pi_shear(kt,kz)=ptt_bj*ch*sh+pzz_bj*tau*ch*tau*sh
  pi_shear(kz,kt)=pi_shear(kt,kz)
  pi_shear(kx,kz)=ptx_bj*sh+pxz_bj*tau*ch
  pi_shear(kz,kx)=pi_shear(kx,kz)
  pi_shear(ky,kz)=pty_bj*sh+pyz_bj*tau*ch
  pi_shear(kz,ky)=pi_shear(ky,kz)
  pi_shear(kz,kz)=ptt_bj*sh**2+pzz_bj*(tau*ch)**2
    
  !the division by tau comes from the Jacobian
!  ds_t=(ds_tau*ch-sh*ds_eta/tau)/tau
!  ds_x=ds_x/tau
!  ds_y=ds_y/tau
!  ds_z=(-sh*ds_tau+ch*ds_eta/tau)/tau
  dsigma_hyp(0)=(ds_tau*ch-sh*ds_eta/tau)
  dsigma_hyp(1)=ds_x
  dsigma_hyp(2)=ds_y
  dsigma_hyp(3)=(-sh*ds_tau+ch*ds_eta/tau)

  betasq=vb(1)*vb(1)+vb(2)*vb(2)+vb(3)*vb(3)*tau*tau
  ub(0)=1.0/sqrt(1-betasq)
  do i=1, 3
    ub(i)=ub(0)*vb(i)    
  end do 

  um(0)=ub(0)*ch+tau*sh*ub(3)
  um(1:2)=ub(1:2) 
  um(3)=ub(0)*sh+tau*ch*ub(3)
 ! um_back(0)=um(0)
 ! um_back(1:3)=-um(1:3)

  Lambda_boost=0.d0
!  Lambda_boost(0,:)=um_back(:)
!  Lambda_boost(1:3,0)=um_back(1:3)
  Lambda_boost(0,:)=um(:)
  Lambda_boost(1:3,0)=um(1:3)

   do i=1,3
     do j=1,3
        if(i .ne. j) then 
          delta_kr=0
        else
          delta_kr=1.d0
        end if  
        Lambda_boost(i,j)=delta_kr+um(i)*um(j)/(um(0)+1.d0)
     end do
   end do

!now we create the matrix to boost back contravariant objects to the LRF
  Lambda_boost_back=0.d0
  um_back(0)=um(0)
  um_back(1:3)=-um(1:3)
  Lambda_boost_back(0,:)=um_back(:)
  Lambda_boost_back(1:3,0)=um_back(1:3)
 
   do i=1,3
     do j=1,3
        if(i .ne. j) then 
          delta_kr=0
        else
          delta_kr=1.d0
        end if 
        Lambda_boost_back(i,j)=delta_kr+um_back(i)*um_back(j)/(um_back(0)+1.d0)        
     end do
   end do

!now we boost the shear tensor back in the LRF
   tmparr=matmul(pi_shear,transpose(Lambda_boost_back))
   pi_shear_lrf=matmul(Lambda_boost_back,tmparr)
     
!we also boost the hypersurface in the LRF
!however, being a covariant object, the boost back matrix is just the forward boost matrix for the contravariant objects

   dsigma_hyp_lrf=matmul(Lambda_boost,dsigma_hyp) 
      
   ! we extract roughly the maximum value of the cooper_frye formula in this cell
   if(viscosity) then
     call coarse_loop_visc(ipart, dsigma_hyp_lrf, massimo, Tfo, edens, press)
   else
     call coarse_loop(ipart, dsigma_hyp_lrf, massimo)
   end if
   if(massimo .le. 0) then
     feedback=1
     return
   end if
  
  distr_max=ffmax(m,mu,Tfo,bf)
  do while (.true.)
   do while (.true.)
     if (counttrials > maxtrials) then 
        print *, "Too many trials... Particle discarded."
        feedback=1
        return
     end if
     counttrials=counttrials+1
     ! extract a random momentum in a box   		---
     call random_number(r)
     p=r*PMAX_BOX !p
     pguess=p*p*ff(p,m,mu,Tfo,bf)
     call random_number(r)
     if(pguess .gt. r*distr_max*1.01) then
        exit
     end if
   end do
   p0=sqrt(m2+p*p)
   !now we sample phi and theta or phi and rapidity 
   call random_number(r)
   phi=r*2.d0*GREEK_PI
   call random_number(r)
   theta=acos(-1+2*r)
   px=p*sin(theta)*cos(phi)
   py=p*sin(theta)*sin(phi)
   pz=p*cos(theta)
   temp_mom(0)=p0
   temp_mom(1)=px
   temp_mom(2)=py
   temp_mom(3)=pz

   pmuds=(p0*dsigma_hyp_lrf(0)+px*dsigma_hyp_lrf(1)+py*dsigma_hyp_lrf(2)+pz*dsigma_hyp_lrf(3))/p0

   if(viscosity) then
    !the minus sign in the next equation comes from the covariant p_{\mu} p_{\nu} components
    mom_shear_prod=p0*(p0*pi_shear_lrf(kt,kt)-px*pi_shear_lrf(kt,kx)-py*pi_shear_lrf(kt,ky)-pz*pi_shear_lrf(kt,kz))&
    -px*(p0*pi_shear_lrf(kx,kt)-px*pi_shear_lrf(kx,kx)-py*pi_shear_lrf(kx,ky)-pz*pi_shear_lrf(kx,kz))&
    -py*(p0*pi_shear_lrf(ky,kt)-px*pi_shear_lrf(ky,kx)-py*pi_shear_lrf(ky,ky)-pz*pi_shear_lrf(ky,kz))&
    -pz*(p0*pi_shear_lrf(kz,kt)-px*pi_shear_lrf(kz,kx)-py*pi_shear_lrf(kz,ky)-pz*pi_shear_lrf(kz,kz))
    pmutot=pmuds*(1+(1+bf*ff(p,m,mu,Tfo,bf))*mom_shear_prod/(2*Tfo**2*(edens+press)))
    if(pmutot .gt. 0) then
      call random_number(r)
      CoopFrye=r*massimo
      if(pmutot .gt. massimo) then
         write(*,*)
         write(*,*) "Error in finding the maximum of the final acceptance distribution"
         write(*,*)
         call exit(7)
      end if
      if (pmutot .gt. CoopFrye) then   !everything is okay, we store the momentum
        feedback=0
        ppp=matmul(Lambda_boost,temp_mom) 
        if(D2) then
          yrap=0.5d0*log((ppp(0)+ppp(3))/(ppp(0)-ppp(3)))
          yrap_minus_eta=yrap+eta !we add the random pseudorapidity
          mt=sqrt(m2+ppp(1)**2+ppp(2)**2)
          ppp(0)=mt*cosh(yrap_minus_eta)
          ppp(3)=mt*sinh(yrap_minus_eta) 
        end if
        return
      end if
     end if  !end do pmuds > 0
   else !ideal case
    if(pmuds .gt. 0) then
     call random_number(r)
     CoopFrye=r*massimo
     if(pmuds .gt. massimo) then
        write(*,*)
        write(*,*) "Error in finding the maximum of the final acceptance distribution"
        write(*,*)
        call exit(7)
     end if
     if (pmuds .gt. CoopFrye) then   !everything is okay, we store the momentum
       feedback=0
       ppp=matmul(Lambda_boost,temp_mom) 
       if(D2) then
         yrap=0.5d0*log((ppp(0)+ppp(3))/(ppp(0)-ppp(3)))
         yrap_minus_eta=yrap+eta !we add the random pseudorapidity
         mt=sqrt(m2+ppp(1)**2+ppp(2)**2)
         ppp(0)=mt*cosh(yrap_minus_eta)
         ppp(3)=mt*sinh(yrap_minus_eta) 
       end if
       return
     end if
    end if  !end do pmutot > 0
   end if !end if viscosity
 end do !main loop
 end subroutine work_momentum_generator
 ! !     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   
  subroutine coarse_loop( ipart, hyper, massimo)
! ! it finds the Maximum of the  formula 
! ! performing a coarse loop oven the momentum space 
! ! for a given dV, x^mu and u^mu
  use eos  
  implicit none
  integer, intent (in):: ipart 		    ! particle identifyier 
  real, dimension(1:4), intent(in) :: hyper     ! f.o. hyp elements
  real, intent(out) :: massimo ! cooper-frye maximum value
  real eta,tau,ch,sh,x,y,t,z,ds_tau,ds_t,ds_x,ds_y,ds_eta,ds_z,mass,p0,m2
  real dpt,dphi,dcosth
  real pmuds
  integer i,j,k
 
  m2=(pdata(ipart)%mass)**2

!  tau=hyper(1)
!  x=hyper(2)
!  y=hyper(3)
!  eta=hyper(4)
!  ds_tau=hyper(10)
!  ds_x=hyper(11)
!  ds_y=hyper(12)
!  ds_eta=hyper(13)
  ds_t=hyper(1)
  ds_x=hyper(2)
  ds_y=hyper(3)
  ds_z=hyper(4)

!  ch=cosh(eta)
!  sh=sinh(eta)

  !the division by tau comes from the Jacobian
  !ds_t=(ds_tau*ch-sh*ds_eta/tau)/tau
  !ds_x=ds_x/tau
  !ds_y=ds_y/tau
  !ds_z=(-sh*ds_tau+ch*ds_eta/tau)/tau
!  ds_t=(ds_tau*ch-sh*ds_eta/tau)
!  ds_x=ds_x
!  ds_y=ds_y
!  ds_z=(-sh*ds_tau+ch*ds_eta/tau)

  !we convert the dsigma_mu in cartesian coordinates
  !dmom_arr(1,:,:,:) is p^2
  massimo=-1.d10
 
  do k=1,N_dp
     do j=1,N_dphi-1
        do i=1,N_dth
           p0=sqrt(m2+dmom_arr(1,i,j,k))
           pmuds=(p0*ds_t+dmom_arr(2,i,j,k)*ds_x+dmom_arr(3,i,j,k)*ds_y+dmom_arr(4,i,j,k)*ds_z)/p0
           if (pmuds > massimo) then
              massimo=pmuds
           endif
        end do
     end do
  end do

!  write(*,*) "massimo:", massimo
  
  massimo=massimo*MAXfactor
  return
 end subroutine coarse_loop 

  subroutine coarse_loop_visc( ipart, hyper_lrf, massimo, temp, edens, press)
! ! it finds the Maximum of the  formula 
! ! performing a coarse loop oven the momentum space 
! ! for a given dV, x^mu and u^mu
  use eos 
  implicit none
  integer, intent (in):: ipart 		    ! particle identifyier 
  real, dimension(1:4), intent(in) :: hyper_lrf     ! f.o. hyp elements in the LRF
  real, intent(out) :: massimo ! cooper-frye maximum value
  real eta,tau,ch,sh,x,y,t,z,ds_tau,ds_t,ds_x,ds_y,ds_eta,ds_z,mass,p0,m2,ptau,peta,mt,rap,pt2
  real dpt,dphi,dcosth, stat_fac
  real pmuds, pmutot, mom_shear_prod, press, temp, edens, mu_chempot
  integer i,j,k
  real, dimension(1:27) :: hyp !full f.o. hypersurface
  real :: px, py, pz, pmod
     
  m2=(pdata(ipart)%mass)**2
  stat_fac=pdata(ipart)%bf
  call get_mu(ipart, edens, mu_chempot)

  ds_t=hyper_lrf(1)
  ds_x=hyper_lrf(2)
  ds_y=hyper_lrf(3)
  ds_z=hyper_lrf(4)
  
  tau=hyp(1)
  eta=hyp(4)
  
    
  massimo=-1.d10
 
  do k=1,N_dp
     do j=1,N_dphi-1
        do i=1,N_dth
           p0=sqrt(m2+dmom_arr(1,i,j,k))
           px=dmom_arr(2,i,j,k)
           py=dmom_arr(3,i,j,k)
           pz=dmom_arr(4,i,j,k)
           pmuds=(p0*ds_t+px*ds_x+py*ds_y+pz*ds_z)/p0
           mom_shear_prod=p0*(p0*pi_shear_lrf(kt,kt)-px*pi_shear_lrf(kt,kx)-py*pi_shear_lrf(kt,ky)-pz*pi_shear_lrf(kt,kz))&
           -px*(p0*pi_shear_lrf(kx,kt)-px*pi_shear_lrf(kx,kx)-py*pi_shear_lrf(kx,ky)-pz*pi_shear_lrf(kx,kz))&
           -py*(p0*pi_shear_lrf(ky,kt)-px*pi_shear_lrf(ky,kx)-py*pi_shear_lrf(ky,ky)-pz*pi_shear_lrf(ky,kz))&
           -pz*(p0*pi_shear_lrf(kz,kt)-px*pi_shear_lrf(kz,kx)-py*pi_shear_lrf(kz,ky)-pz*pi_shear_lrf(kz,kz))       
           pmutot=pmuds*(1+(1+stat_fac/(exp((p0-mu_chempot)/temp)-stat_fac))*mom_shear_prod/(2*temp**2*(edens+press)))
           if (pmutot > massimo) then
              massimo=pmutot
           endif
        end do
     end do
  end do

!  write(*,*) "massimo:", massimo
  
  massimo=massimo*MAXfactor
  return
 end subroutine coarse_loop_visc

real function ff(p,m,mu,T,bf)
implicit none
real, intent(in) :: p,m,mu,T,bf

ff=1.d0/(exp((sqrt(p*p+m*m)-mu)/T)-bf)

end function ff

real(8) function ffmax_section(m,mu,T,bf)
implicit none
real(8), intent(in) :: m,mu,T,bf
real(8), parameter :: dp_lim=1.d-5
real(8) pleft,p1,p2,pright,fleft,f1,f2,fright,dp
integer steps

steps=0
pleft=0.
fleft=0.
p1=PMAX_BOX/3.
f1=p1*p1*ff(p1,m,mu,T,bf)
p2=PMAX_BOX*2/3.
f2=p2*p2*ff(p2,m,mu,T,bf)
pright=PMAX_BOX
fright=PMAX_BOX**2*ff(PMAX_BOX,m,mu,t,bf)
dp=PMAX_BOX

do while (dp>dp_lim)
if(f1 .eq.f2) then !since we have only one maximum, it must be in the middle
    pleft=p1
    pright=p2
end if

if(f1 .gt. f2) then
    pright=p2
else
    pleft=p1
end if

dp=pright-pleft
fleft=pleft**2*ff(pleft,m,mu,T,bf)
p1=pleft+dp/3.
f1=p1*p1*ff(p1,m,mu,T,bf)
p2=pright-dp/3.
f2=p2*p2*ff(p2,m,mu,T,bf)
fright=pright**2*ff(pright,m,mu,T,bf)
steps=steps+1
end do

ffmax_section=(pleft+dp/2)**2 * ff(pleft+dp/2,m,mu,T,bf)

end function ffmax_section

real(8) function ffmax(m,mu,T,bf)
implicit none
real(8), intent(in) :: m,mu,T,bf
!the accuracy at which to stop, the bracketing interval for the momentum
real(8), parameter :: dp_lim=1.d-6, pmin=0.d0, pmax=3.d0
real(8) :: f,f1,p,dp,pl,ph,dpold,pold
real(8) :: E2,K,Z,E
integer steps
integer, parameter :: maxsteps=100

!we want to return the maximum of g=p^2*f(Boltzmann/Dirac/Bose)
!we proceed by finding p such that g'=0. Since we can factor a p, actually we find the root of f such that g'=p*f, i.e. f=g'/p (p>0)
!f1=f' is the derivative of f
!we use a mixed Newton-Raphson / bisection method, very close to that described
!in the Numerical Recipes (rtsafe)

steps=1
!write(*,*) "m,mu,T and bf:",m,mu,T,bf
!we know that the maximum is for p=0 and that the function is positive there,
!while at pmax is negative
ph=pmin
pl=pmax
p=(pl+ph)/2.d0
dpold=pmax-pmin
dp=dpold
E2=m**2+p**2
E=sqrt(E2)
K=exp((E-mu)/T)
Z=K-bf
f=1.d0/Z*(2.d0-p**2*K/(Z*E*T))
f1=-p*K*(-2*p**2*K*E + p**2*Z*E -T*Z*p**2 + 4*E2*K*T - 4*E2*T*bf)/(E2*E*T**2 * Z**3)
do steps=1,maxsteps
!   write(*,*) "Step: ", steps, "p=",p, "f=",f, "f1=",f1,"dp=",dp
   dpold=dp
   if(( (((p-ph)*f1-f) * ((p-pl)*f1-f)) .gt. 0) .or. (abs(2.d0*f) .gt. (dpold*f1))) then !Newton's method failing or slow convergence
     dp=(ph-pl)/2.d0
     p=pl+dp
     if(pl .eq. p) exit !the increment is negligible
   else
     dp=f/f1
     pold=p
     p=p-dp
     if(pold .eq. p) exit !again, we consider the case in which the finite machine precision can lead to a negligble increment
   end if
   if(abs(dp) .lt. dp_lim) exit
   E2=m**2+p**2
   E=sqrt(E2)
   K=exp((E-mu)/T)
   Z=K-bf
   f=1/Z*(2.d0-p**2*K/(Z*E*T))
   f1=-p*K*(-2*p**2*K*E + p**2*Z*E -T*Z*p**2 + 4*E2*K*T - 4*E2*T*bf)/(E2*E*T**2 * Z**3)
   if(f .lt. 0) then
     pl=p
   else
     ph=p
   end if
end do

!we compute the maximum corresponding to p
ffmax=p**2*ff(p,m,mu,T,bf)

end function ffmax

    subroutine init_random_seed()
    implicit none
    integer, allocatable :: seed(:)
    integer :: n, pid,statflag, tt(2), xx, dates(8), times, i
    integer(8) :: clock_val
    logical :: clock_used
#ifdef __INTEL_COMPILER
    integer, external :: getpid
#endif
    call random_seed(size = n)
    allocate(seed(n))
    open(75, file='/dev/urandom', access='stream', form='UNFORMATTED', action="read", status="old", iostat=statflag)
    if( statflag .eq. 0) then
        read(75) seed
        close(75)
        call random_seed(put=seed)
        write(*,*) "Seed from /dev/urandom:",seed(:)
    else
        call system_clock(clock_val)
        if(clock_val .ne. 0) then
            tt = transfer(clock_val, tt)
            clock_used=.true.
        else
            clock_used=.false.
            call date_and_time(values=dates)
            times = (dates(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                  + dates(2) * 31_8 * 24 * 60 * 60 * 1000 &
                  + dates(3) * 24 * 60 * 60 * 60 * 1000 &
                  + dates(5) * 60 * 60 * 1000 &
                  + dates(6) * 60 * 1000 &
                  + dates(7) * 1000 &
                  + dates(8)
            tt = transfer(times, tt)
        end if
        !warning, the effective randomness of this procedure has not been carefully evaluated
        !hopefully, it should be used only when /dev/urandom is not available
        xx = ieor(tt(1), tt(2))
        pid = getpid() + 104412563 
        xx = ieor(xx, pid)
        if (n .ge. 3) then
            seed(1) = tt(1) + 375412561
            seed(2) = tt(2) + 475414183
            seed(3) = pid
            if (n .gt. 3) then
              seed(4:) = xx + 131 * (/ (i, i = 0, n - 4) /)
            end if
        else
            seed = xx + 131 * (/ (i, i = 0, n - 1 ) /)
        end if
        call random_seed(put=seed)
        if(clock_used) then
          write(*,*) "Seed from system_clock call:",seed(:)
          write(*,*) "System clock value:",clock_val
        else
          write(*,*) "Seed from date_and_time call:",seed(:)
          write(*,*) "Built time value:",times
        end if
      end if

    end subroutine init_random_seed   



 end module work

 
! ! ! !  ***********************************************************************
! ! ! !  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ! ! !  -----------------------------------------------------------------------
! ! ! !  §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
! ! ! !  ----------------------------------------------------------------------- 
! ! ! !  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ! ! !  ***********************************************************************
 
program jyu_sample    
    use common
    use numlib, only: poisson
    use eos
    use work

    
    implicit none
  
   
    real harvest
    ! counter of particle produced in the local cell
    integer produced 
    !Global particle density in the cell (adimensional)
    real global_density 
    !particle density in the cell (adimensional) for each specie  
    ! counter of produced particles
    integer pindex 
    
    ! error flags 
    integer AllocateStatus, eof
    
!   hypersurface array info:
!   IDEAL (hysu_ide)
!   1-tau		2-x	3-y	4-eta
!   5-rho		6-vx	7-vy	8-veta		9-prex
!   10-dV0	11-dV1		12-dV2		13-dV3
!   VISOCUS (hysu_visco)
!   1-bulk
!   2-pi^xy    3-pi^xz    4-pi^yz    5-pi^xx    6-pi^yy    7-pi^zz
!   8-pi^tt    9-pi^tx    10-pi^ty    11-pitz
!   12-fracden=1.0/(temp*temp*(eng+prex))
    real, dimension(:,:), allocatable ::  hysu
        
    integer, parameter :: samples = 25
    
    character*1024 bitbucket
    
    integer i,j,k 
    
    integer ih 		! index of frozen cells
    integer iout 	! index for time-step
    integer cit 	! total amount of frozen cells that we read. We can check it against frozen_cells
    integer feedback ! it checks if the momentum has been assigned correctly
    integer num_arguments
    character(len=32) outfile
    character(len=160) infile
    character visclabel
    integer iterations
    integer(kind=8) :: totcells,partcells
    real, dimension(1:27) :: tmp_input_arr
    integer :: indx_produced,tot_iterations
    integer, parameter :: max_samp_part=20000
    integer, dimension(1:max_samp_part) :: index_array
    real, dimension(1:max_samp_part,1:4) :: pos_array, mom_array
    integer,dimension(:), allocatable :: nspecie

    ih=0
    iout=1
    cit=0
    
001 format (I3,ES14.7)      
    !***** START ******!
        
 
    num_arguments=command_argument_count()
    if (num_arguments .eq. 3) then
       call get_command_argument(1, infile)
       write(*,*) "I will use the input file: ",infile
       LID_in=index(infile, ' ')-1
       call get_command_argument(2, outfile)
       write(*,*) "I will use the output file: ",outfile
       LID_out=index(outfile, ' ')-1
       call get_command_argument(3, visclabel)
       if(visclabel == "0") then
         viscosity=.false.
       else if (visclabel == "1") then
         viscosity=.true.
       else
         write(*,*) "There is something wrong in the arguments..."
         write(*,*) "Syntax: ./sample inputfile output_directory 0/1 ( 0 = no viscosity, 1 = with viscosity )" 
         call exit(1)
       end if 
    else
       write(*,*) "There is something wrong in the arguments..."
       write(*,*) "Syntax: ./sample inputfile output_directory 0/1 ( 0 = no viscosity, 1 = with viscosity )" 
       call exit(1)
    end if

    
    
    ! READ THE HYPERSURFACE FILE  and see if it's T or e -----
    !     open(unit=14,status='old',file=input1(1:li)//'.dat',form='unformatted', iostat=filerror, access='stream')
    open(unit=22,status='old',file=infile(1:LID_in),form='formatted', iostat=filerror)
    call check_file(filerror, infile(1:LID_in))
    write(*,*) "Counting freezeout file line numbers"
    DO WHILE (.true.)
       read(22,*,iostat=eof)
       if(eof == 0) then
         frozen_cells=frozen_cells+1
       else
         exit
       end if
    end do
    frozen_cells=frozen_cells-2 !we remove the two lines of comments
    write(*,*) "There are ",frozen_cells," cells"
    close(22)

    allocate(hysu(frozen_cells,1:27), STAT = AllocateStatus)
    if (AllocateStatus .ne. 0) then
       write(*,*) "Unable to allocate the hysu array. I quit."
       call exit(3)
    end if
    hysu=0.0

    call count_particles()
    call allocate_arrays()
    call read_particles()
    call read_chemical_potential_file()
  
    
    allocate( particle_density(1:npart_main), particle_density_common(1:npart_main), nspecie(1:npart_main),STAT = AllocateStatus)
    if (AllocateStatus .ne. 0) then
       write(*,*) "Unable to allocate the particle_density and/or the particle_density_common and/or the nspecie arrays. I quit."
       call exit(3)
    end if
    hysu=0.0


       print *, "Reading freeze-out hypersurface data"
       open(unit=22,status='old',file=infile(1:LID_in),form='formatted', iostat=filerror)
       read(22,*,iostat=eof)
       read(22,*,iostat=eof)
       ih=1
       DO WHILE (.true.) ! until we find the end of file 
            read(22, *,  iostat=eof) (tmp_input_arr(j), j=1, 27)
            if(eof .eq. 0) then
              hysu(ih,:)=tmp_input_arr(:)
              ih=ih+1
            else
              ih=ih-1
              exit
            end if 
       END DO
       close(22)
       if(ih .ne. frozen_cells) then
         write(*,*) "Hey, there is something weird going on with the number of cells in the f.o. hypersurface..."
         call exit(2)
       end if
    print *,"Done"
    print *,"Initializing the random number generator"
    call init_random_seed()
    print *,"Opening output file ", outfile(1:LID_out)
    open(unit=23,status='replace',file=outfile(1:LID_out),form='formatted', iostat=filerror)
    write(23,*) "Samples: ", samples
    
    print *,"°°°°°°°°°°°°°°°°°°°°°       here we go!      °°°°°°°°°°°°°°°°°°°°°"   
 
    nspecie(:)=0
    call fill_N_arrays()
    if(D2) then
      tot_iterations=int(2*eta_pseudorap_max)*samples
    else
      tot_iterations=samples
    end if
    totcells=int(frozen_cells,8)*tot_iterations
    do iterations=1,tot_iterations
     write(*,*) "Sampling number ",iterations, "( of ", tot_iterations," )"
     ih=0
     indx_produced=1
     index_array=0.
     pos_array=0.
     mom_array=0.
     do while(.true.)
      ih=ih+1
      partcells=int(ih,8)+int(frozen_cells,8)*int(iterations-1 ,8)
      if(mod(partcells,1000) .eq. 0) then 
          write(*,*) "Cells done so far: ", partcells, " of ", totcells,&
          &" , i.e. ",float(partcells)*100/(totcells)," %"
      end if
      produced = 0          
!      do while (produced == 0 ) ! we try until we find a particle
	! choose a random position over the hypersurface
	!call work_choose_random_index(ih, frozen_cells)	
	! calculate the global particle density in it 	
	call work_global_density(hysu(ih,:), particle_density, global_density)	
	! and discriminate whether a particle is produced or not	 
	! if the global density is above 0.01, we sample it through a Poisson distribution
	if (global_density > 0.01) then 
	  produced=poisson(global_density)   
	! if the global density is below 0.01, we sample it straightforwardly
	else if ((global_density .le. 0.01) .and. (global_density>0.0)) then
	  call random_number(harvest)
	  if (harvest .le. global_density) produced=1	  
	else  
	  produced=0	  
	endif 	
!      end do
!       print *, "here", produced
      do j=1, produced	
!       WARNING the index j_p of the produced particles is incremented in assign_features
	call work_assign_features(pindex,particle_density, global_density)
!	print *, j_p, name(particle_index(j_p))
	call work_momentum_generator(pindex, hysu(ih,:), part_coord(:), part_mom(:), feedback)
        if(feedback .eq. 0) then
          nspecie(pindex)=nspecie(pindex)+1
          index_array(indx_produced)=pindex
          pos_array(indx_produced,:)=part_coord(:)
          mom_array(indx_produced,:)=part_mom(:)
          indx_produced=indx_produced+1
          if(indx_produced>max_samp_part) then
            write(*,*) "Not enough space to store informations about sampled particles..."
            write(*,*) "Please, increase the max_samp_part parameter."
            call exit(3)
         end if
        end if
      end do
      if(ih .eq. frozen_cells) exit
    END DO
    do j=1,indx_produced-1
          write(23,*) index_array(j), pdata(index_array(j))%PDG_id, pdata(index_array(j))%name, pos_array(j,:), mom_array(j,:)
    end do
   end do 

   close(23)
   print *,"Output file closed"
    
    do j=2,npart_main
        write(*,*) pdata(j)%name, nspecie(j)
    end do  
    print *,"All done"
 end program jyu_sample
