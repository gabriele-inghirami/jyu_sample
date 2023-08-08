! *  *  *  *  *
! *
! *  JYU_SAMPLE
! *
! *  Version: 0.5
! *                          
! *  Date (DD/MM/YYYY): 08/08/2023
! *                  
! *  File: jyu_sample/settings.f90 
! *                           
! *  Author: 
! *  Gabriele Inghirami (GSI Helmholtzzentrum für Schwerionenforschung GmbH)
! *  E-mail: g.inghirami@gsi.de
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

module settings

  implicit none

  ! non technical user options

  ! format of the hypersurface
  ! 1=CORNELIUS, 2=ECHO-QGP 2D+1 ideal, 3=ECHO-QGP 2D+1 visco, 4=ECHO-QGP 3D+1 ideal, 5=ECHO-QGP 3D+1 visco
  integer, parameter :: hyp_format=1

  ! boost invariant sampling.
  ! if D2=.true. then the pseudorapidity of the particle is randomly sampled from a uniform distribution from
  ! -eta_pseudorap_max to +eta_pseudorap_max
  ! if D2=.false. then we perform a 3D sampling and the pseudorapidity is simply that of the freeze-out hypersurface cell
  logical, parameter :: D2=.false.

  ! maximum absolute value of the pseudorapidity added to the sampled particle momentum
  real(8), parameter :: eta_pseudorap_max=2.d0

  ! tecnical expert options

  ! the minimum value of the Cooper-Frye maximum probability to sample a particle
  ! very small values allow to try to sample a particle even if it is very unlikely that it will be created,
  ! however, if too small, the marginal improved accuracy comes at the expense of overall longer sampling times
  real(8), parameter :: minCF=1.e-12

  ! maximum value of the momentum to be used in the integral of the distrbution function for the sampling
  ! the higher this value, the longer the integration time; however, if it is too small, it might introduce biases
  ! when sampling the momentum of heavy hadrons
  real(8), parameter :: PMAX_BOX=12

  ! below this particle rest mass value [GeV] we use a high order expansion
  ! with Bessel functions to compute the particle density
  real(8), parameter :: MASS_LIMIT_FOR_HIGHER_EXPANSION=0.6

  ! the number of terms in the high order expansion with Bessel functions 
  integer, parameter :: CUT_PI=10

  ! the resolution of the coarse grained grid to get a quick,
  ! but raw estimate the maximum of the distribution to sample a particle
  ! N_dp = num. of cells of momentum
  ! N_dphi = num. of cells of the azimuthal angle phi (0 - 2\pi)
  ! N_dth = num. of cells of the polar angle theta (0 - \pi)
  integer, parameter :: N_dp=121, N_dphi=81, N_dth=41

  ! parameter to increase the estimated value of the maximum of the sampling distribution found with the coarse method
  real(8), parameter :: MAXfactor=1.2

end module settings
