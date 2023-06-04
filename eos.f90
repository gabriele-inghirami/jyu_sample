! *  *  *  *  *
! *
! *  JYU_SAMPLE
! *
! *  Version: 0.5
! *                          
! *  Date (DD/MM/YYYY): 08/10/2019
! *                  
! *  File: jyu_sample/eos.f90 
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
! *  - INSPIRED BY A C++ PROGRAM WRITTEN IN 2013 BY:
! *  Hannu Holopainen (Franfkurt University and FIAS - Germany)
! *  (some parts of the present code are actually just Fortran translations of that code)
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

module eos
implicit none

real(8), parameter :: lifetime_to_be_frozen=1.d1
real(8), parameter :: hbarc=0.197327054064
integer, parameter :: max_daughters=5
logical, parameter :: rescale_branching_ratios_to_max_3_part_decay=.true. !it removes the 4 and 5 particles decays

integer :: npart_main, npart_decays, nbaryons_in_the_list, npart_main_in_the_list, npart_decays_in_the_list, nstable, nfrozen
integer :: filerror, alloc_error, readerror
!character(len=*),parameter :: particle_file="pdg05.dat", chempot_file="s95p-PCE175-v1_pichem1.dat"
character(len=*),parameter :: particle_file="particle_and_eos_data/pdg05.dat"
character(len=*),parameter :: chempot_file="particle_and_eos_data/s95p-PCE175-v1_pichem1.dat"
type particle_infos
integer :: PDG_id
integer :: spin_degeneracy
integer :: baryon_number
integer :: strangeness
integer :: charmness
integer :: bottomness
integer :: isospin
integer :: charge
integer :: ndecays
integer :: decays_offset_index
integer :: bf !(-1 for bosons and +1 for fermions)
character(len=27) :: name
real(8) :: mass
real(8) :: width
logical :: frozen
logical :: stable
end type particle_infos

type(particle_infos), allocatable, dimension(:) :: pdata

integer, allocatable, dimension(:) :: decay_PDG_id_parent, decay_id_parent, num_daughters
real(8), allocatable, dimension(:) :: branching_ratio
integer, allocatable, dimension(:,:) :: mc, mc_pid

real(8) :: min_edens_in_table, delta_edens, max_edens_in_table
integer :: edens_table_npoints, frozen_particles_to_read

real(8), allocatable, dimension(:,:) :: chempot_table
real(8), allocatable, dimension(:) :: edens_array

contains

subroutine count_particles()
implicit none
character(len=22) :: pname
real(8) :: mass, width
integer :: i, PDG_id, deg, nb, ns, nc, nbot, isospin, charge, ndecays
integer :: PDG_id_parent, n_daug, mc1, mc2, mc3, mc4, mc5
real(8) :: dec_frac
!111 format (i8,2x,a22,1x,f8.0,2x,f8.0,8(1x,i2))
!222 format (i8,1x,i2,2x,f5.0,7x,5(i8))
111 format (i8,2x,a22)

npart_main=0
npart_main_in_the_list=0
npart_decays=0
npart_decays_in_the_list=0
nbaryons_in_the_list=0
open(unit=20,file=particle_file,form='formatted', IOSTAT=filerror, action="READ")
if(filerror .ne. 0) then
    write(*,*) "Problem in opening the file "//particle_file
    write(*,*) "I quit"
    call exit(2)
end if

do while(readerror .eq. 0)
   read(20,111,IOSTAT=readerror,advance='no') PDG_id, pname
!DDD   write(*,*) PDG_id, pname, mass,  width,  deg,  nb,  ns,  nc,  nbot,  isospin,  charge,  ndecays
   if(PDG_id .eq. 7)  then
     close(20)
     return
   else
     read(20,*)  mass,  width,  deg,  nb,  ns,  nc,  nbot,  isospin,  charge,  ndecays
     npart_main=npart_main+1
     npart_decays=npart_decays+ndecays
     npart_main_in_the_list=npart_main_in_the_list+1
     npart_decays_in_the_list=npart_decays_in_the_list+ndecays
     if(nb .eq. 1) then
       nbaryons_in_the_list=nbaryons_in_the_list+1
       npart_main=npart_main+1 !we create room for an additional anti-baryon
       npart_decays=npart_decays+ndecays
     end if
     if(rescale_branching_ratios_to_max_3_part_decay) then !with this option, we do not count the decays in more than 3 particles
      do i=1,ndecays
!        read(20,222) PDG_id_parent, n_daug, dec_frac, mc1, mc2, mc3, mc4, mc5
        read(20,*) PDG_id_parent, n_daug, dec_frac, mc1, mc2, mc3, mc4, mc5
!DDD        write(*,*)   PDG_id_parent, n_daug, dec_frac, mc1, mc2, mc3, mc4, mc5
        if(PDG_id_parent .ne. PDG_id) then
          write(*,*) "Particle PDG ID mismatch between daughter's parent", PDG_id, " and ", PDG_id_parent
          write(*,*) "Sorry, but I cannot continue if you don't fix the error."
          call exit(2)
        end if
        if(n_daug .gt. 3) then
          npart_decays=npart_decays-1
          npart_decays_in_the_list=npart_decays_in_the_list-1
          if(nb .eq. 1) then !if the particle is a baryon, we remove one decay also from its antibaryon
            npart_decays=npart_decays-1
          end if
        end if
      end do     
     else
      do i=1,ndecays
!        read(20,222) PDG_id_parent, n_daug, dec_frac, mc1, mc2, mc3, mc4, mc5
        read(20,*) PDG_id_parent, n_daug, dec_frac, mc1, mc2, mc3, mc4, mc5
!DDD        write(*,*)   PDG_id_parent, n_daug, dec_frac, mc1, mc2, mc3, mc4, mc5
        if(PDG_id_parent .ne. PDG_id) then
          write(*,*) "Particle PDG ID mismatch between daughter's parent", PDG_id, " and ", PDG_id_parent
          write(*,*) "Sorry, but I cannot continue if you don't fix the error."
          call exit(2)
        end if
      end do     
     end if!end of rescale_branching_ratios_to_max_3_part_decay condition
   end if
end do

write(*,*) "Error in reading "//particle_file
close(20)
call exit(3)
  
end subroutine count_particles

subroutine read_particles()
implicit none
integer :: i, j, k, offset, pid, did, pid_bar, did_bar
character(len=22) :: pname, tmp_name
real(8) :: lifetime
!111 format (i8,2x,a22,1x,f8.0,2x,f8.0,8(1x,i2))
111 format (i8,2x,a22)
!222 format (i8,1x,i2,2x,f5.0,6x,5(1x,i7))
!222 format (i8,1x,i2,2x,f5.0,7x,5(i8))
integer :: tmp_decay_PDG_id_parent, tmp_num_daughters, ndecays_to_eliminate
real(8) :: tmp_branching_ratio, sum_br, sum_verification
real(8), dimension(1:5) :: tmp_mc


open(unit=20,file=particle_file,form='formatted', IOSTAT=filerror, action="READ")
if(filerror .ne. 0) then
    write(*,*) "Problem in opening the file "//particle_file
    write(*,*) "I quit"
    call exit(2)
end if

offset=0
pid=0   
did=0
nfrozen=0
nstable=0
do i=1,npart_main_in_the_list
   pid=pid+1
   read(20,111,IOSTAT=readerror,advance='no') pdata(pid)%PDG_id, pdata(pid)%name
   read(20,*) pdata(pid)%mass, pdata(pid)%width, pdata(pid)%spin_degeneracy, pdata(pid)%baryon_number, pdata(pid)%strangeness,&
   & pdata(pid)%charmness, pdata(pid)%bottomness, pdata(pid)%isospin,pdata(pid)%charge, pdata(pid)%ndecays
   !we replace the spaces with underscores
   do k=1,len(trim(adjustl(pdata(pid)%name)))
      if(pdata(pid)%name(k:k) .eq. ' ') pdata(pid)%name(k:k)="_"
   end do

   if(readerror .ne. 0)  then
     close(20)
     write(*,*) "Error in reading "//particle_file
     close(20)
     call exit(3)
   end if 
   if(pdata(pid)%width .gt. 0) then
     lifetime = hbarc/pdata(pid)%width
   else
     lifetime=1.d100
   end if
   if(lifetime .ge. lifetime_to_be_frozen) then
      pdata(pid)%frozen = .true.
      nfrozen=nfrozen+1
   else 
      pdata(pid)%frozen = .false.
   end if 
   if (pdata(pid)%baryon_number .ne. 0) then
       pdata(pid)%bf=-1
   else
       pdata(pid)%bf=1
   endif
   pdata(pid)%decays_offset_index=offset
   offset=offset+pdata(pid)%ndecays

   if(rescale_branching_ratios_to_max_3_part_decay) then !with this option, we do not count the decays in more than 3 particles
    ndecays_to_eliminate=0
    sum_br=0
    do j=1,pdata(pid)%ndecays 
      did=pdata(pid)%decays_offset_index+j-ndecays_to_eliminate
      read(20,*) tmp_decay_PDG_id_parent, tmp_num_daughters, tmp_branching_ratio, tmp_mc(1:5)
      if(tmp_num_daughters .gt. 3) then
        offset=offset-1
        ndecays_to_eliminate=ndecays_to_eliminate+1
        cycle
      end if
      decay_PDG_id_parent(did)=tmp_decay_PDG_id_parent
      num_daughters(did)=tmp_num_daughters
      branching_ratio(did)=tmp_branching_ratio
      sum_br=sum_br+tmp_branching_ratio
      mc(1:5,did)=tmp_mc(1:5)
      decay_id_parent(did)=pid
      num_daughters(did)=abs(num_daughters(did)) !sometimes in the PDG table this value is negative...
      if(decay_PDG_id_parent(did) .ne. pdata(pid)%PDG_id) then
         write(*,*) "PDG id mismatch between ",decay_PDG_id_parent(did)," and ", pdata(pid)%PDG_id
         write(*,*) "Quitting."
         call exit(3)
      end if
     end do
     pdata(pid)%ndecays=pdata(pid)%ndecays-ndecays_to_eliminate
     if(sum_br .gt. 0) then
       sum_verification=0.d0 !now we rescale the branching ratios and we want to be sure that their sum is 1
       do j=1,pdata(pid)%ndecays 
          did=pdata(pid)%decays_offset_index+j
          branching_ratio(did)=branching_ratio(did)/sum_br !we rescale the branching ratios to so get 1 as their sum
          sum_verification=sum_verification+branching_ratio(did)
       end do
       if(abs(sum_verification-1.d0) .gt. 1.d-9) then
         write(*,*) "Something went wrong when rescaling the branching ratios after removing 4 and 5 particle decays"
         write(*,*) "I quit"
         call exit(3)
       end if 
     end if
    else
     do j=1,pdata(pid)%ndecays 
      did=pdata(pid)%decays_offset_index+j
      decay_id_parent(did)=pid
!      read(20,222) decay_PDG_id_parent(did), num_daughters(did), branching_ratio(did), mc(1:5,did)
      read(20,*) decay_PDG_id_parent(did), num_daughters(did), branching_ratio(did), mc(1:5,did)
      num_daughters(did)=abs(num_daughters(did)) !sometimes in the PDG table this value is negative...
!DDD      write(*,*) decay_PDG_id_parent(did), num_daughters(did), branching_ratio(did), mc(1:5,did)
      if(decay_PDG_id_parent(did) .ne. pdata(pid)%PDG_id) then
         write(*,*) "PDG id mismatch between ",decay_PDG_id_parent(did)," and ", pdata(pid)%PDG_id
         write(*,*) "Quitting."
         call exit(3)
      end if
     end do     
   end if !end of rescale_branching_ratios_to_max_3_part_decay condition
   
   if((pdata(pid)%ndecays .eq. 1) .and. (num_daughters(did) .eq. 1) .and. (mc(1,did) .eq. decay_PDG_id_parent(did))) then
     pdata(pid)%stable=.true.
     nstable=nstable+1
   end if
   if(pdata(pid)%baryon_number .eq. 1) then !we create the antibaryon
     pid_bar=pid
     pid=pid+1
     pdata(pid)%PDG_id=-pdata(pid_bar)%PDG_id
     pdata(pid)%frozen=pdata(pid_bar)%frozen
     if(pdata(pid)%frozen) nfrozen=nfrozen+1
     pdata(pid)%stable=pdata(pid_bar)%stable
     if(pdata(pid)%stable) nstable=nstable+1
     pdata(pid)%name='anti_'//pdata(pid_bar)%name
     pdata(pid)%mass=pdata(pid_bar)%mass
     pdata(pid)%width=pdata(pid_bar)%width
     pdata(pid)%spin_degeneracy=pdata(pid_bar)%spin_degeneracy
     pdata(pid)%baryon_number=-pdata(pid_bar)%baryon_number
     pdata(pid)%strangeness=-pdata(pid_bar)%strangeness
     pdata(pid)%charmness=-pdata(pid_bar)%charmness
     pdata(pid)%bottomness=-pdata(pid_bar)%bottomness
     pdata(pid)%isospin=-pdata(pid_bar)%isospin
     pdata(pid)%charge=-pdata(pid_bar)%charge
     pdata(pid)%ndecays=pdata(pid_bar)%ndecays
     pdata(pid)%bf=pdata(pid_bar)%bf
     pdata(pid)%decays_offset_index=offset
     offset=offset+pdata(pid)%ndecays
     do j=1,pdata(pid)%ndecays 
        did=pdata(pid)%decays_offset_index+j
        did_bar=pdata(pid_bar)%decays_offset_index+j
        decay_id_parent(did)=pid
        decay_PDG_id_parent(did)=-decay_PDG_id_parent(did_bar)
        num_daughters(did)=num_daughters(did_bar)
        branching_ratio(did)=branching_ratio(did_bar)
        mc(:,did)=-mc(:,did_bar)
     end do     
   end if
end do

close(20)

call remap_mc_mc_pid()

end subroutine read_particles

subroutine remap_mc_mc_pid
implicit none
integer :: i,j

do j=1,npart_decays
   do i=1,5
      if(mc(i,j) .ne. 0) then
         mc_pid(i,j)=convert_PDG_id_pid(mc(i,j))
      end if
   end do
end do

end subroutine remap_mc_mc_pid

integer function convert_PDG_id_pid(pid)
implicit none
integer i,pid

convert_PDG_id_pid=0
do i=1,npart_main
   if(pdata(i)%PDG_id .eq. pid) then
     convert_PDG_id_pid=i
     return
   end if
   !we have this conditions because in the case of antibaryons we flipped all pdg id of decays,
   ! but there are bosons which have only positive id, like the photon or the neutral pion
   if(pdata(i)%PDG_id .eq. abs(pid)) then
     convert_PDG_id_pid=i
   end if
end do

if(convert_PDG_id_pid .ne. 0) return

write(*,*) "Unable to convert PDG id into internal id. I am forced to quit."
write(*,*) "Given PDG id was: ",pid
call exit(3)

end function


subroutine read_chemical_potential_file()
implicit none
integer :: i,j,k, rindx
real(8), allocatable, dimension(:) :: tmp_chempot

open(unit=20,file=chempot_file,form='formatted', IOSTAT=filerror, action="READ")
if(filerror .ne. 0) then
    write(*,*) "Problem in opening the file "//chempot_file
    write(*,*) "I quit"
    call exit(2)
end if

read(20,*)  min_edens_in_table
read(20,*)  delta_edens, edens_table_npoints
read(20,*)  frozen_particles_to_read

max_edens_in_table=min_edens_in_table+delta_edens*edens_table_npoints

allocate(chempot_table(1:edens_table_npoints,1:npart_main),tmp_chempot(1:frozen_particles_to_read),STAT=alloc_error)
if(alloc_error .ne. 0) then
   write(*,*) "Error in allocating memory for the chempot_table or tmp_chempot arrays"
   call exit(2)
end if

chempot_table=0.

do rindx=0,edens_table_npoints-1
   i=edens_table_npoints-rindx
   read(20,*,IOSTAT=readerror) tmp_chempot(1:frozen_particles_to_read)
   if(readerror .ne. 0) then
      write(*,*) "Error in reading "//chempot_file
      call exit(2)
   end if
   k=1
   do j=2,npart_main
      if(pdata(j)%frozen) then
        chempot_table(i,j)=tmp_chempot(k)
        k=k+1
        if(k .gt. frozen_particles_to_read) exit
      else
        chempot_table(i,j)=0.d0
      end if
   end do
end do

allocate(edens_array(1:edens_table_npoints), STAT=alloc_error)
if(alloc_error .ne. 0) then
   write(*,*) "Error in allocating memory for edens_array "
   call exit(2)
end if

do i=1,edens_table_npoints   
   edens_array(i)=(i-1)*delta_edens+min_edens_in_table
end do

close(20)

write(*,*) chempot_file//" read."
end subroutine read_chemical_potential_file

subroutine allocate_arrays()

implicit none

allocate(pdata(1:npart_main),STAT=alloc_error)
if(alloc_error .ne. 0) then
  write(*,*) "Error in allocating the pdata array of structures containing particle infos. Quitting."
  call exit(2)
end if

allocate(decay_PDG_id_parent(1:npart_decays), decay_id_parent(1:npart_decays), num_daughters(1:npart_decays), mc(5,1:npart_decays),&
   & mc_pid(1:5,1:npart_decays), branching_ratio(1:npart_decays), STAT=alloc_error)
if(alloc_error .ne. 0) then
   write(*,*) "Error in allocating the arrays containing particle decays infos. Quitting."
   call exit(2)
end if

end subroutine allocate_arrays


recursive subroutine get_mu(pid, edens, mu)
implicit none

integer, intent(in) :: pid
real(8), intent(out) :: mu 
real(8), intent(inout) :: edens
real(8) :: x, mu_i, mu_daug
integer :: index_low, index_high, i, j, k, indx, dec_indx

if ((edens .lt. min_edens_in_table) .or. (edens .ge. max_edens_in_table)) then
   write(*,*) "Error, I cannot compute the chemical potential with the value of the energy density that you provided,"
   write(*,*) "as it is outside the boundaries (or exactly at the top) of the PCE EOS table. I quit."
   !it would be possible to include the case in which edens .eq. max_edens_in_table, but, since it is unlikely to use it,
   !probably it is better to avoid to include a dedicated if condition
   call exit(3)       
end if

index_low = int(floor((edens - min_edens_in_table)/delta_edens))+1 !our first array index is 1
index_high = index_low + 1

x = (edens - edens_array(index_low))/delta_edens

if(pdata(pid)%frozen) then
  mu=(1.0 - x)*chempot_table(index_low,pid) + x*chempot_table(index_high,pid)
!DDD  write(*,*) edens, x, index_low, pid, index_high, mu
else
  mu = 0.d0
  do i=1,pdata(pid)%ndecays
     indx=pdata(pid)%decays_offset_index+i
     do j=1,num_daughters(indx)
        dec_indx = mc_pid(j,indx)
        if(pdata(dec_indx)%frozen) then
           mu_i = (1.0 - x)*chempot_table(index_low,dec_indx) + x*chempot_table(index_high,dec_indx)
           mu = mu + mu_i*branching_ratio(indx)
        else
           call get_mu(dec_indx,edens,mu_daug)
           mu = mu + mu_daug*branching_ratio(indx)
        end if
     end do
   end do
end if

end subroutine get_mu


end module eos
