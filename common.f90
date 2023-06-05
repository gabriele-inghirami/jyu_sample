! *  *  *  *  *
! *
! *  JYU_SAMPLE
! *
! *  Version: 0.5
! *                          
! *  Date (DD/MM/YYYY): 08/10/2019
! *                  
! *  File: jyu_sample/common.f90 
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
! *  - INSPIRED BY A C++ PROGRAM WRITTEN IN 2013 BY:
! *  Hannu Holopainen (Franfkurt University and FIAS - Germany)
! *
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

  module constants
  implicit none
      real(8), parameter :: PIGRECO=3.14159265358979323846
      real(8), parameter :: HBARC=0.197326
      real(8), parameter :: NOR_DENS=1.0/(2.0*PIGRECO*PIGRECO*HBARC*HBARC*HBARC)
      real(8), parameter :: TWOPIhbarc_cube = (2.0*PIGRECO*HBARC)**3
      !the D2 option selects boost invariant sampling.
      !If D2=.true. then the pseudorapidity of the particle is randomly sampled from a uniform distribution from
      !-eta_pseudorap_max to +eta_pseudorap_max
      !If D2=.false. then we perform a 3D sampling and the pseudorapidity is simply that of the freeze-out hypersurface cell
      logical, parameter :: D2=.false.
      real(8), parameter :: eta_pseudorap_max=2.d0
      !1=CORNELIUS, 2=ECHO-QGP 2D+1 ideal, 3=ECHO-QGP 2D+1 visco, 4=ECHO-QGP 3D+1 ideal, 5=ECHO-QGP 3D+1 visco
      integer, parameter :: hyp_format=1
  end module constants
 ! !------------------------------------------------------------------  
  module common ! OVERALL
    use  constants
    implicit none

      integer dimension_flag
! 	directory and file names for the i/o
      character*32 inputdir,file,input,input1,outdir,output      
      integer LID_in, LID_out
      integer  frozen_cells,  produced_particles
      real Energy_int,particle_energy
    
      integer npart
      integer seed_settings, trueseed
      
contains
 
! !------------------------------------------------------------------  
 subroutine check_file(Status, filename)
  implicit none
  integer, intent(in) :: Status
  character(len=*), intent(in) :: filename
  if (Status .ne. 0) then
    print *, "******  ERROR  ******"
    print *, "the file named:    " , filename
    print *, "cannot be opened. I am forced to quit!"
    
    print *, "filename test", filename
    
    call exit(1)
  end if
 end subroutine check_file
 ! !------------------------------------------------------------------   
 end module common
