! *  *  *  *  *
! *
! *  JYU_SAMPLE
! *
! *  Version: 0.5
! *                          
! *  Date (DD/MM/YYYY): 08/10/2019
! *                  
! *  File: jyu_sample/decays.f90 
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

module utilities

implicit none
real(8), parameter :: PIGREEK=3.14159265358979323846
real(8), parameter :: HBARC=0.197326

contains

pure function lambda(l0,l1,l2)
implicit none
real(8) lambda
real(8), intent(in) :: l0,l1,l2

lambda=l0**4 + l1**4 + l2**4 - 2.d0*(l0**2)*(l1**2) - 2.d0*(l0**2)*(l2**2) - 2.d0*(l1**2)*(l2**2)
return
end function lambda

subroutine Lboost(unboosted,boosted,vel)
implicit none
real(8), intent(in), dimension(0:3) :: unboosted
real(8), intent(out), dimension(0:3) :: boosted
real(8), dimension(0:3) :: u_mu
real(8), intent(in), dimension(1:3) :: vel
real(8) :: betasq, delta_kr
real(8), dimension(0:3,0:3) :: Lambda_boost
integer :: i,j,k

betasq=vel(1)**2+vel(2)**2+vel(3)**2
u_mu(0)=1.0/sqrt(1-betasq)
do i=1, 3
   u_mu(i)=u_mu(0)*vel(i)
end do

Lambda_boost(0,:)=u_mu(:)
Lambda_boost(1:3,0)=u_mu(1:3)

do i=1,3
   do j=1,3
      if(i .ne. j) then 
        delta_kr=0
      else
        delta_kr=1.d0
      end if  
      Lambda_boost(i,j)=delta_kr+u_mu(i)*u_mu(j)/(u_mu(0)+1.d0)
   end do
end do

boosted=matmul(Lambda_boost,unboosted)

end subroutine Lboost
     
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

subroutine check_file(fstatus, filename)
  implicit none
  integer, intent(in) :: fstatus
  character(len=*), intent(in) :: filename
  if (fstatus .ne. 0) then
    write(*,*) "******  ERROR  ******"
    write(*,*) "the file named:    " , filename
    write(*,*) "cannot be opened. I am forced to quit!"

    call exit(1)
  end if
end subroutine check_file

end module utilities

module decay_routines
    use eos
    implicit none
    integer, dimension(:), allocatable :: pid_in, pid_out
    real(8), dimension(:,:), allocatable :: pos_array_in, mom_array_in, pos_array_out, mom_array_out
    integer :: ipart_in, ipart_out

contains
 
    subroutine compute_decays(pid,location,mom)
    implicit none
    integer :: j,did, pid
    real(8) :: sum_of_branching_ratios, probability_br
    real(8), dimension(1:4) :: location, mom

    if(pdata(pid)%stable) then !the particle is stable
       ipart_out=ipart_out+1
       pid_out(ipart_out)=pid
       pos_array_out(ipart_out,:)=location(:)
       mom_array_out(ipart_out,:)=mom(:)
    else
       !we decide the decay channel
       sum_of_branching_ratios=0
       call random_number(probability_br)
       do j=1,pdata(pid)%ndecays
          did=pdata(pid)%decays_offset_index+j
          sum_of_branching_ratios=sum_of_branching_ratios+branching_ratio(did)
          if(probability_br .lt. sum_of_branching_ratios) then
            exit
          end if
       end do
       if(num_daughters(did) .eq. 2) then
         call compute_2_part_decay(pid,location, mom, did)
       else if(num_daughters(did) .eq. 3) then
         call compute_3_part_decay(pid,location, mom, did)
       else
         write(*,*) "Sorry, but currently I am not able to deal with decays with more than 3 daughter particles"
         write(*,*) "Plese, check if you compiled the code with the parameter inside eos.f90:"
         write(*,*) "rescale_branching_ratios_to_max_3_part_decay=.true."
         write(*,*) "Sorry, I prefer to stop here"
         call exit(3)
       end if
    end if
    end subroutine compute_decays

    subroutine compute_part_decay(m0,m1,m2,mom1,mom2)
    use utilities
    implicit none
    real(8) :: m0,m1,m2,p,phi,r,costheta,theta
    real(8), dimension(1:4) :: mom1, mom2
    
    p=sqrt(lambda(m0,m1,m2))/(2*m0)
    call random_number(r)
    phi=2.d0*PIGREEK*r
    call random_number(r)
    costheta = -1.d0 + 2.d0*r
    theta = acos(costheta)
    mom1(2) = p*sin(theta)*cos(phi) !px of particle 1
    mom1(3) = p*sin(theta)*sin(phi) !py of particle 2
    mom1(4) = p*cos(theta) !pz of particle 3
    mom2(2:4)=-mom1(2:4)
    mom1(1)=sqrt(m1**2+sum(mom1(2:4)**2))
    mom2(1)=sqrt(m2**2+sum(mom2(2:4)**2))
 


    end subroutine compute_part_decay

    subroutine compute_2_part_decay(pid,location, mom0, did)
    use utilities
    implicit none
    integer :: id1,id2
    integer :: did, pid
    real(8), dimension(1:4) :: location, mom0, mom1, mom2, b_mom1, b_mom2
    !mom0 is the four momentum of the decaying particle, mom1 and mom2 the four momenta of the daughters in the decaying part LRF
    !b_mom1 and b_mom2 are the four momenta of the daughters in the computational ref. frame, i.e. mom1 and mom2 after the L. boost
    real(8) :: m0,m1,m2
    !m0, m1 and m2 are the rest masses of the decaying particle and its daughters
    real(8), dimension(1:3) :: vel !3 velocity of the decaying particle

    id1=mc_pid(1,did)
    id2=mc_pid(2,did) 

    m0=pdata(pid)%mass
    m1=pdata(id1)%mass
    m2=pdata(id2)%mass
    if((m1+m2) .gt. m0) then !here in the future we might improve this condition
      m0 = m1 + m2 + 0.002
      mom0(1)=sqrt(m0**2+mom0(2)**2+mom0(3)**2+mom0(4)**2) !mom0(1)=energy,mom0(2)=px,mom0(3)=py,mom0(4)=pz of the decaying part.
    end if
    vel(1:3)=mom0(2:4)/mom0(1)
    call compute_part_decay(m0,m1,m2,mom1,mom2) !here we get the four momenta of the daughter particles in the decaying part LRF
    call Lboost(mom1,b_mom1,vel)
    call Lboost(mom2,b_mom2,vel)
    call compute_decays(id1,location,b_mom1) !here we let the daughters decays, if they are not stable
    call compute_decays(id2,location,b_mom2) !the position does not change, this might be changed when using an afterburner
    !DDD temporary section for debugging purposes
    if((mom1(1) .le. mom1(4)) .or. (mom2(1) .le. mom2(4))) then
      write(*,*) "Error in 2 part decay"
      write(*,*) "Last ipart_out:", ipart_out
      write(*,*) mom0(:)
      write(*,*) mom1(:)
      write(*,*) mom2(:)
      call exit(4)
    end if
    !DDD end of temporary debugging section
    
    end subroutine compute_2_part_decay

    subroutine compute_3_part_decay(pid,location, mom0, did)
    use utilities
    implicit none
    integer :: id1,id2,id3
    integer :: did, pid
    real(8), dimension(1:4) :: location, mom0, mom1, mom2, mom3, mom23, b_mom1, b_mom2_rel23, b_mom2, b_mom3_rel23, b_mom3
    real(8) :: m0,m1,m2,m3, m23, m23_max_sq, m23_min_sq
    real(8) :: acceptance_value, sample_prob, r, distromax
    real(8), dimension(1:3) :: vel, vel_23
    !m23 is the invariant mass of the particles 2-3, m23_max_sq its max value squared, i.e. (m0-m1)^2
    ! m23_min_sq its min value squared, i.e. (m2+m3)^2
    !for the other variables, please, read compute_2_part_decay comments

    id1=mc_pid(1,did)
    id2=mc_pid(2,did) 
    id3=mc_pid(3,did) 

    m0=pdata(pid)%mass
    m1=pdata(id1)%mass
    m2=pdata(id2)%mass
    m3=pdata(id3)%mass

    if((m1+m2+m3) .gt. m0) then !here in the future we might improve this condition
      m0 = m1 + m2 + m3 + 0.002
      mom0(1)=sqrt(m0**2+mom0(2)**2+mom0(3)**2+mom0(4)**2) !mom0(1)=energy,mom0(2)=px,mom0(3)=py,mom0(4)=pz of the decaying part.
    end if
    vel(1:3)=mom0(2:4)/mom0(1)

    m23_max_sq=(m0-m1)**2
    m23_min_sq=(m2+m3)**2
 
    distromax= sqrt((2.d0*m0**4 + 2.d0*m1**4 + 6.d0*((m0*m1)**2))*(m0**4 + m1**4 + m2**4 + m3**4 + 6.d0*((m0*m1)**2)))/&
             &(m0*m0*(m2*m2 + m3*m3))

    acceptance_value=0.d0
    sample_prob=1.d0
    do while (sample_prob .gt. acceptance_value)
      call random_number(r)
      m23=sqrt(m23_min_sq+r*(m23_max_sq-m23_min_sq))
      acceptance_value=sqrt(lambda(m0,m1,m23)*lambda(m23,m2,m3))/((m0*m23)**2)
      call random_number(r)
      sample_prob=r*distromax 
    end do
    
    call compute_part_decay(m0,m1,m23,mom1,mom23) !here we get the four momenta of the daug part 1 and 23 in the dec part 0 LRF
    call compute_part_decay(m23,m2,m3,mom2,mom3) !here we get the four momenta of the daug part 2 and 3 in the dec part 23 LRF
    call Lboost(mom1,b_mom1,vel) ! we boost the particle 1 in the computational ref frame
    vel_23(1:3)=mom23(2:4)/mom23(1)
    call Lboost(mom2,b_mom2_rel23,vel_23) !we boost the particle 2 in the 1-23 RF
    call Lboost(mom3,b_mom3_rel23,vel_23) !we boost the particle 2 in the 1-23 RF
    call Lboost(b_mom2_rel23,b_mom2,vel) !we boost the particle 2 in the computational ref frame 
    call Lboost(b_mom3_rel23,b_mom3,vel) !we boost the particle 2 in the computational ref frame 
    call compute_decays(id1,location,b_mom1) !here we let the daughters decays, if they are not stable
    call compute_decays(id2,location,b_mom2) !the position does not change, this might be changed when using an afterburner
    call compute_decays(id3,location,b_mom3) !the position does not change, this might be changed when using an afterburner
     
    !DDD temporary section for debugging purposes
    if((mom1(1) .le. mom1(4)) .or. (mom2(1) .le. mom2(4)) .or. (mom3(1) .le. mom3(4))) then
      write(*,*) "Error in 3 part decay"
      write(*,*) "Last ipart_out:", ipart_out
      write(*,*) mom0(:)
      write(*,*) mom1(:)
      write(*,*) mom2(:)
      write(*,*) mom3(:)
      call exit(4)
    end if
    !DDD end of temporary debugging section

    end subroutine compute_3_part_decay


    subroutine write_results(hs,out_filename,end_signal)
    implicit none
    character(len=*) :: out_filename
    integer :: j
    integer :: filerror
    logical, save :: first_time=.true.
    logical :: end_signal
    character(len=120) :: hs
    integer, save :: wt

    if(first_time) then
      first_time=.false.
      open(unit=23,status='replace',file=out_filename,form='formatted', iostat=filerror)
      write(23,'(120a)') hs
      wt=0
    end if

    wt=wt+1
    write(*,*) "I am writing the ouput for the ",wt,"time"
    do j=1,ipart_out
          write(23,*) pid_out(j), pdata(pid_out(j))%PDG_id, pdata(pid_out(j))%name,&
         &            pos_array_out(j,:), mom_array_out(j,:)
    end do

    if(end_signal) then
      close(23) !if it is the last writing, we close the output file
    else
      ipart_out=0 !if there are still particle to write, we reset the index of the output file
    end if
    end subroutine 


end module decay_routines

 
program decay    
    use utilities
    use decay_routines
    
    implicit none
    integer i,j,k,events
!   maximum number of particles that are read from the input file before starting to process them
!   please, tune it according to the memory available in your computer
    integer, parameter :: max_block_of_particle_to_read=5000000
!   maximum number of processed particles that are kept in memory before write them in the output file
    integer :: max_block_of_particle_to_write=max_block_of_particle_to_read*3
    integer :: numpart
    integer feedback ! it checks if the momentum has been assigned correctly
    integer num_arguments
    character(len=500) outfile
    character(len=500) infile
    character(len=132) header_string
    character modelabel
    character(10) :: eventlabel
    real(8) harvest
    integer :: alloc_status,eof
    integer :: bitbucket_PDGid
    integer :: LID_in, LID_out
    character(len=100) :: bitbucket_name

 
    num_arguments=command_argument_count()
    if (num_arguments .eq. 2) then
       call get_command_argument(1, infile)
       write(*,*) "I will read the input particle list from file: ",infile
       LID_in=index(infile, ' ')-1
       call get_command_argument(2, outfile)
       write(*,*) "I will write the output particle list in file: ",outfile
       LID_out=index(outfile, ' ')-1
!       call get_command_argument(3, modelabel)
!       read(modelabel,"(I1)") mode
!       write(*,*) "I will use decay mode ",mode
!       mode=0
!       write(*,*) "Decay mode set 0, as it is the only option available at the moment"
!       call get_command_argument(4, eventlabel)
!       read(eventlabel,*) nevents
!       write(*,*) "I will perform ",nevents," decays"
    else
       write(*,*) "There is something wrong in the arguments..."
       write(*,*) "Syntax: ./decay.exe inputfile outputfile" 
       call exit(1)
    end if

    max_block_of_particle_to_write=max_block_of_particle_to_read

!   we allocate the arrays to contain the particle data
    allocate(pid_in(1:max_block_of_particle_to_read), pid_out(1:max_block_of_particle_to_write),&
   &pos_array_in(1:max_block_of_particle_to_read,1:4), mom_array_in(1:max_block_of_particle_to_read,1:4),&
   &pos_array_out(1:max_block_of_particle_to_write,1:4), mom_array_out(1:max_block_of_particle_to_write,1:4),STAT=alloc_status)
    if(alloc_status .ne. 0) then
      write(*,*) "Sorry, but I am not able to allocate the arrays to store the informations about the particles."
      write(*,*) "I suggest you either to reduce the number of decay repetitions when invoking the program or reduce the parameter"
      write(*,*) "max_block_of_particle_to_read in decay.f90, recompile the code and try again"
      call exit(2)
    end if

    pid_in=0
    pos_array_in=0.d0
    mom_array_in=0.d0


!   we get informations about the particles    
    call count_particles()
    call allocate_arrays()
    call read_particles()
    call read_chemical_potential_file()
    call init_random_seed()
    
    open(unit=22,status='old',file=infile(1:LID_in),form='formatted', iostat=filerror)
    call check_file(filerror, infile(1:LID_in))
    write(*,*) "Reading input datafile:", infile(1:LID_in)
!   we assume that the first lines of the file is a comment with the number of f.o. hyp. samplings
!   we will save it in a string and then we will rewrite it in the output file
    read(22,'(132a)',iostat=eof) header_string
    ipart_in=1  !ipart_in is already at the correct index, it is increased after each reading
    ipart_out=0 !ipart_out must be increased before writing
    DO WHILE (.true.)
       read(22,*,iostat=eof) pid_in(ipart_in), bitbucket_PDGid, bitbucket_name,&
      &pos_array_in(ipart_in,:), mom_array_in(ipart_in,:)
       if(eof == 0) then
         ipart_in=ipart_in+1
         if(ipart_in .gt. max_block_of_particle_to_read) then
           do i=1,ipart_in-1 
              call compute_decays(pid_in(i),pos_array_in(i,:),mom_array_in(i,:))
           end do
           call write_results(header_string,outfile(1:LID_out),.false.)
           pid_in=0
           pos_array_in=0.d0
           mom_array_in=0.d0
           ipart_in=1
         end if
       else
         ipart_in=ipart_in-1
         do i=1,ipart_in
            call compute_decays(pid_in(i),pos_array_in(i,:),mom_array_in(i,:))
         end do
         call write_results(header_string,outfile(1:LID_out),.true.)
         exit
       end if
    end do
    close(22)
    write(*,*) "Program decay ended."
  
 end program decay
