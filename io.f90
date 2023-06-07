! *  *  *  *  *
! *
! *  JYU_SAMPLE
! *
! *  Version: 0.5
! *                          
! *  Date (DD/MM/YYYY): 08/10/2019
! *                  
! *  File: jyu_sample/io.f90 
! *                           
! *  Author: 
! *  Gabriele Inghirami+ (University of Jyvaskyla and Helsinki Institute of physics- Finland)
! *  E-mail: gabriele.g.inghirami@jyu.fi
! *  in collaboration with:
! *  Harri Niemi+ (University of Jyvaskyla and Helsinki Institute of physics- Finland)
! *  +: actually, given the amount of reused code, the real author is V. Rolando, see later note
! * 
! *  Copyright - Important attribution note:
! *               
! *  THIS PROGRAM HEAVILY REUSES CODE:
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

module io
  use constants
  implicit none 
  integer, dimension(:), allocatable :: id_list
  integer ID_start, ID_stop, los, part_list
  integer,parameter :: antibaryons=1
  integer, parameter :: chempot=0
  contains
! !******************************************************************
! !******************************************************************
! !***********************   INPUT  ROUTINES   **********************
! !******************************************************************  
! 
  subroutine io_read_particlelist()
	use common
	use common_thermal
	use init
	implicit none
	integer filerror, i, AllocateStatus, ipt, ipart
	integer j, k
	integer dummy_int
	
	real temperature(2)
	real temp_mu(2,maxpar)
	real weight
	integer check_mu
	real Tfo_MeV
        real dpt_cl
!   
        real, dimension(maxpar) :: R_m, R_mu, R_width, R_isospin
        integer,dimension(maxpar) :: R_pdg_number , R_bf, R_g
        integer, dimension(maxpar) :: R_baryon, R_strange, R_bottom, R_charme, R_charge, R_n_decays
!       
        integer R_mc_dec, R_daughters, R_mc1,  R_mc2 , R_mc3  ,R_mc4  ,R_mc5, R_npart
        real R_ratio
	integer hypercharge

	character*24 R_name(maxpar) !particle name      

	R_m=0.d0
        R_mu=0.d0
	R_pdg_number=0
        R_bf=0
        R_g=0
	R_baryon=0
        R_strange=0
        R_charme=0
        R_charge=0

	Tfo_MeV=1000.0*Tfo
	check_mu=0

	open(unit=30,status='old',file='../eos_data/pdglist.txt',form='formatted', iostat=filerror)
	call check_file(filerror, '../eos_data/pdglist.txt')	

! 
! ! skip the first 35 lines
   do i=1, 35
	read (30, *) 
	end do
	i=1
	do while (i.lt.maxpar) 
	  read (30, *, end=166) R_pdg_number(i)
! 	  PRINT *, R_pdg_number(i)
	  READ (30, *, END=166) R_name(i)
          READ(30,*,END=166) R_m(i), R_width(i), R_g(i), R_baryon(i), R_strange(i), R_charme(i), R_bottom(i), &
                           & R_isospin(i), R_charge(i), R_n_decays(i)
          do j=1,R_n_decays(i)
            read(30,*) R_mc_dec, R_daughters, R_ratio, R_mc1,  R_mc2 , R_mc3,  R_mc4, R_mc5
!           PRINT *, R_mc_dec, R_daughters, R_ratio, R_mc1,  R_mc2 , R_mc3,  R_mc4, R_mc5
          end do
! 	  PRINT *, R_name(i)
	  if (R_baryon(i).NE.0) then
	      R_bf(i)=-1
	    else 
	      R_bf(i)=1
	  endif    	  
	  hypercharge= R_strange(i)+R_charme(i)+ R_baryon(i) !top and bottom quarks are out of range
	  i=i+1
	  end do
	  print *, "you have more particle than expected, change maxpar in the common file" 
	  call exit
 166    continue
	R_npart=i-1
! 	WRiTE(*,2001) "INDEX","PDG","NAME","MASS","DEG", "B", "e", "S", "b/f", "mu"
! 	PRINT *, "---------------------------------------------------"
	CLOSE(30)

	if (chempot .eq. 1) then 
	  open(unit=31,status='old',file='../eos_data/chemical_potential.txt',form='formatted', iostat=filerror)	
	    if(0.ne.filerror) then 
	      print *, "*** I am forced to quit *** "
	      print *, "cannot find the file ../eos_data/chemical_potential.txt"
	      print *, "*** I am forced to quit *** "
	      call exit(8)
	    end if
  ! ! 	print *, Tfo
	  R_mu=0.0
	  check_mu=0
	  read (31,*) temperature(2), (temp_mu(2, j), j=1, R_npart)
	  do while (check_mu.lt.1) 
	    read (31,*, end=167) temperature(1), (temp_mu(1, j), j=1, R_npart)
	    if ((temperature(1)-Tfo_MeV)*(Tfo_MeV-temperature(2)) .GE. 0.0) then 
	      check_mu=2
	      weight=(temperature(2)-Tfo_MeV)/(temperature(2)-temperature(1))
	      if (weight<0.0) then 
		print *, "error (weight<0.0)  in io-3D.f08 subroutine read_particlelist ", weight
		call exit (5)
	      endif
	      do j=1, R_npart
		R_mu(j)=temp_mu(2,j)-weight*(temp_mu(2,j)-temp_mu(1,j))
	      end do
	    endif
	    do j=1, R_npart
	      temp_mu(2, j)=temp_mu(1, j)
	    end do
	    temperature(2)=temperature(1)
	  end do
167	  continue
	  CLOSE(31)
	else if (chempot .eq. 0 ) then 
	  R_mu(:)=0.0
	  check_mu=2
	else 
	  print *, "There is something wrong somewhere, I found chempot=",chempot, " but it can only be 0 or 1"
	  print *, "I am going to quit now"
	  call exit	  
      	endif 

	if (check_mu==0) then
	  print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	  print *, "Hey, something is wrong in routine read_particlelist,"
	  print *, "I cannot interpolate the chemical_potential."
	  print *, "please check again the file ../eos_data/chemical_potential.txt"
	  print *, "I AM ABORTING THE SIMULATION"
	  call exit (5)
	endif 
! ! 
	k=1
	do i=1, R_npart
	   m(k)=R_m(i)
           mu(k)=R_mu(i)
           pdg_number(k)=R_pdg_number(i)
           bf(k)=R_bf(i)
           g(k)=R_g(i)
           baryon(k)=R_baryon(i)
           strange(k)=R_strange(i)
           charme(k)=R_charme(i)
           charge(k)=R_charge(i)
           name(k)=R_name(i)
           mu(k)=R_mu(i)
           if( k .eq. maxpar) then
               write(*,*) "Sorry, but you to allocate more space for the particle arrays by increasing maxpar in common.f90"
               call exit(2)
           end if
           k=k+1
           if(((antibaryons .eq. 1) .or. ( caller .eq. 1 ) ) .and. (R_baryon(i) /= 0)) then
             m(k)=          m(k-1)
             pdg_number(k)=-pdg_number(k-1)
             bf(k)=         bf(k-1)
             g(k)=          g(k-1)
             baryon(k)=    -baryon(k-1)
             strange(k)=   -strange(k-1)
             bottom(k)=    -bottom(k-1)
             charme(k)=    -charme(k-1)
             charge(k)=    -charge(k-1)
             n_decays(k)=  -n_decays(k-1)
             name(k)=       'Anti_'//name(k-1)
             mu(k)=        -mu(k-1)
             if( k .eq. maxpar) then
                write(*,*) "Sorry, but you to allocate more space for the particle arrays by increasing maxpar in common.f90"
                call exit(2)
             end if
             k=k+1
           endif
        end do
        npart=k-1

 2000 FORMAT(I6,I10,2X, A25,F10.6,I5,4I3,3x, F10.6)
 2001 FORMAT(A6,A10,3x,A4,22x,A4,7x, A3, 1x, A1, 2(2x, A1),2x, A3, 2x, A3)
 2002 FORMAT(A10,A7,3x,A4,22x,A4,6x, A3, 1x, A3, 3x, A5, 2x, A3)
 2003 FORMAT(I6,I10,2X, A25,F10.6,I5,I3,3x, F10.6)

       	WRiTE(*,*) "Particle list and relative features"
      	WRiTE(*,2002) "INDEX","PDG","NAME","MASS","DEG", "b/f", "mu"
	do i=1, npart	
	write(*,2003) i,pdg_number(i),name(i),m(i),g(i), bf(i), mu(i)
	end do
	PRINT *, "---------------------------------------------------"

  return
  end subroutine io_read_particlelist
! !******************************************************************




! !******************************************************************
! !******************************************************************
! !***********************   OUTPUT  ROUTINES   *********************
! !******************************************************************  
  subroutine io_build_particle_filename(idx_part, string, fname, lix)
    use common, only: outdir, LID_out, name, pdg_number
    
    integer, intent (in) :: idx_part
    character*32, intent(inout) :: string
    character*128, intent(out) :: fname
    integer, intent(out) :: lix
    
    character*32 particle  
    character*3 segno
    character*11 pdg_n
    integer lin, lis, lip
    
    lis=index(string, ' ')-1
    
    particle=name(idx_part)
    lip=index(particle, ' ')-1
    
    if (pdg_number(idx_part)<0) then
      segno='_-_'
    else if (pdg_number(idx_part)>0) then
      segno='_+_'
    else 
      segno='_0_'
    endif    
    
    if (abs(pdg_number(idx_part))<100) then !A
      write(pdg_n, '(6I1, I2,A3)')  0,0,0,0,0,0,abs(pdg_number(idx_part)), segno(1:3)
    else if (abs(pdg_number(idx_part)) .ge. 100 .and. abs(pdg_number(idx_part)) <1000 )  then
      write(pdg_n, '(5I1, I3,A3)')  0,0,0,0,0,  abs(pdg_number(idx_part)), segno(1:3)
    else 
      print *, ' '
      if   (abs(pdg_number(idx_part)) .ge. 1000 .and. abs(pdg_number(idx_part)) <10000 ) then !B
	write(pdg_n, '(4I1, I4,A3)')  0,0,0,0,    abs(pdg_number(idx_part)), segno(1:3)
      else if  (abs(pdg_number(idx_part)) .ge. 10000 .and. abs(pdg_number(idx_part)) <100000 ) then
	write(pdg_n, '(3I1, I5,A3)')  0,0,0,      abs(pdg_number(idx_part)), segno(1:3)
      else 
	print *, ' '
	if (abs(pdg_number(idx_part)) .ge. 100000 .and. abs(pdg_number(idx_part)) <1000000 )  then !C
	  write(pdg_n, '(2I1,I6,A3)') 0,0,         abs(pdg_number(idx_part)), segno(1:3)
	else if (abs(pdg_number(idx_part)) .ge. 1000000 .and. abs(pdg_number(idx_part)) <10000000 )  then
	  write(pdg_n, '(I1,I7,A3)')  0,        abs(pdg_number(idx_part)), segno(1:3)
	else 
	  if (abs(pdg_number(idx_part)) .ge. 10000000 .and. abs(pdg_number(idx_part)) <100000000 )  then !D
	    write(pdg_n, '(I8,A3)')  0,        abs(pdg_number(idx_part)), segno(1:3)
	  else   
	    print *, 'ooops... there was an error in the print_spectra routine.'
	    print *, 'I am sorry, I have to quit. Please contact the developer.'	
! 	    remember that pdg_n is 8+3 here. if you have a bigger pdg number remember to change the strings too
	    call exit (2)
	  endif !D
	endif  !C
      endif	 !B
    endif 	 !A
        
!      128 = 32 // 1 // 32 // 1 // 8+3 // 32 // 4 (=113 < 128)
    fname=outdir(1:LID_out)//'/'//string(1:lis)//'_'//pdg_n(1:11)//particle(1:lip)//'.txt'
    lix=index(fname, ' ')-1
    
    string='                                '
    return
  end subroutine io_build_particle_filename

! !******************************************************************

subroutine io_print_mc_output(totproduced, energy_integral, T0mudSigmamu, oversample,outputfile)
  use common
  use common_thermal, only: baryon, strange, charme, bottom, isospin, charge  
  use common_MC
  

  
implicit none
  integer, intent(in) :: totproduced !! produced particles
  real, intent(in) :: energy_integral, T0mudSigmamu, oversample
  integer filerror, AllocateStatus
  integer j_p, i, k, ipart, irap, iphi
  
  ! -- histogram related
  real c_pt, c_y, c_phi !bin central value
  real low_pt, up_pt, step_pt 
  real low_logpt,up_logpt, step_logpt
  real low_phi, up_phi, step_phi
  real low_y, up_y, step_y  
  real hi_nor
  integer ipt
  
  real deltaphi, deltapt, deltay
!   -- histograms
  integer, dimension (:,:), allocatable :: dN_dy
  integer, dimension (:,:,:), allocatable :: dN_ptdptdy_yconst
  
  ! -- 
  integer lif_pth , lif
  character*32 stringa32,outputfile
  character*128 filename_pt, filename
  !!-- quantum numbers
  integer strangeness, baryon_ness, charmness, charge_ness  
  
  real pc_m(0:3)
  real pm_m(0:3)
  real th, ta, ch3, sh3, pb0, pb3
    
    produced_particles=totproduced 
    print *, "Done... producing output files"
    print *, "WARNING, this subroutine is temporarily not working"
    print *, "This is due to the different dimensions of the part_coord and part_mom arrays"
    !open(unit=13,status='replace',&
    !& file=(outdir(1:LID_out)//'produced_particles_B'//'.txt'),  form='formatted', iostat=filerror)
    !call check_file(filerror, outdir(1:LID_out)//'produced_particles_B'//'.txt')  
    !write (13, *) "oversample", oversample
    !do j_p=1, totproduced
    !  write (13, *) pdg_number(particle_index(j_p)), &
    !  &  name(particle_index(j_p)), &      
    !  & (part_coord(j_p,k), k=0,3), &
    !  & (part_mom(j_p,k), k=0,3)
    !end do
    !close(13) 
    if(outputfile .eq. "") then
      open(unit=13,status='replace',&
      & file=(outdir(1:LID_out)//'produced_particles_M'//'.txt'),  form='formatted', iostat=filerror)
      call check_file(filerror, outdir(1:LID_out)//'produced_particles_M'//'.txt')  
    else
      open(unit=13,status='replace',&
      & file=(outdir(1:LID_out)//outputfile),  form='formatted', iostat=filerror)
      call check_file(filerror, outdir(1:LID_out)//outputfile//'.txt')  
    end if
    write (13, *) "oversample", oversample
    do j_p=1, totproduced
!      th=tanh(part_mom(j_p,3)-part_coord(j_p,3))
!      ch3=cosh(part_coord(j_p,3))
!      sh3=sinh(part_coord(j_p,3))
!      pb0=part_mom(j_p,0)
!      ta=part_coord(j_p,0)
!      pb3=pb0*th/ta
!      pc_m(0)=ta*ch3
!      pc_m(1)=part_coord(j_p,1)
!      pc_m(2)=part_coord(j_p,2)
!      pc_m(3)=ta*sh3
!      pm_m(0)=pb0*ch3 + pb3*ta*sh3
!      pm_m(1)=part_mom(j_p,1)
!      pm_m(2)=part_mom(j_p,2)
!      pm_m(3)=pb0*sh3 + pb3*ta*ch3
      
!      write (13, *) pdg_number(particle_index(j_p)), &
!      & name(particle_index(j_p)),       &
!      & part_coord(j_p,:), 	 &
!      & part_mom(j_p,:)
    end do
    close(13) 

    
    charge_ness=0
    baryon_ness=0
    strangeness=0
    charmness=0
 !   do i=1, totproduced
 !     charge_ness=charge_ness+charge(particle_index(i)) 
 !     baryon_ness=baryon_ness+baryon(particle_index(i)) 
 !     strangeness=strangeness+strange(particle_index(i)) 
 !     charmness=charmness+charme(particle_index(i)) 
 !   end do
    
    write(*,*)  "Cells on the hypersurface:", frozen_cells
    write(*,*)  "Energy of the fluid:", T0mudSigmamu
    write(*,*)  "Energy of the particles: ", energy_integral
    write(*,*)  "Particles on the hypersurface:", totproduced
!   write(*,*)  "Of which"
!    do ipart=1, npart
!      write(*,*)  ipart, name(ipart), "#  ", nspecie(ipart), " with energy", enespecie_m(ipart), "and p0", enespecie_b(ipart)
!    end do
!    write(*,*)  "---------------------------------------------"
!    write(*,*)  "----------      VIOLATIONS       ------------"
    write(*,*)  "Charge:        ", charge_ness
    write(*,*)  "Baryon number: ", baryon_ness
    write(*,*)  "Strangeness:   ", strangeness
    write(*,*)  "Charmness:     ", charmness  
    
    open(unit=36,status='replace' ,file=(outdir(1:LID_out)//'utils_merge'//'.txt'),  form='formatted', iostat=filerror)
    call check_file(filerror, outdir(1:	LID_out)//'utils_merge'//'.txt')
    write(36,*)  "Cells on the hypersurface:", frozen_cells
    write(36,*)  "Energy of the fluid:", T0mudSigmamu
    write(36,*)  "Energy of the particles: ", energy_integral
    write(36,*)  "Particles on the hypersurface:", totproduced
!    write(36,*)  "Of which"
!    do ipart=1, npart
!      write(36,*)  ipart, name(ipart), "#  ", nspecie(ipart), " with energy", enespecie_m(ipart), "and p0", enespecie_b(ipart)
!    end do
    write(36,*)  "---------------------------------------------"
    write(36,*)  "----------      VIOLATIONS       ------------"
    write(36,*)  "Charge:        ", charge_ness
    write(36,*)  "Baryon number: ", baryon_ness
    write(36,*)  "Strangeness:   ", strangeness
    write(36,*)  "Charmness:     ", charmness  
    close(36)
    
return
end subroutine io_print_mc_output

end module io
