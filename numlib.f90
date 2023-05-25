! *  *  *  *  *
! *
! *  JYU_SAMPLE
! *
! *  Version: 0.5
! *                          
! *  Date (DD/MM/YYYY): 08/10/2019
! *                  
! *  File: jyu_sample/numlib.f90 
! *                           
! *  Author: 
! *  Gabriele Inghirami (University of Jyvaskyla and Helsinki Institute of physics- Finland)
! *  E-mail: gabriele.g.inghirami@jyu.fi
! * 
! *  Note: the formulas for the Bessel functions have been taken from
! *  "Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables",
! *  by Milton Abramowitz and Irene Stegun, NBS (1964)
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

module numlib

implicit none
 
contains

       real function I0(x)
        implicit none
        real, intent(in) :: x
        real :: t

        t=x/3.75
        if (x .le. 3.75) then
           I0=1.+3.5156229*t**2.+3.0899424*t**4.+1.2067492*t**6.+0.2659732*t**8.+0.0360768*t**10.+0.0045813*t**12.
        else
           I0=(0.39894228+0.01328592/t+0.00225319/(t**2.)-0.00157565/(t**3.)+000916281/(t**4.)-0.2057706/(t**5.)&
                  &+0.02635537/(t**6.)-0.01647644/(t**7.)+0.00392377/(t**8.))/(sqrt(x)*exp(-x))
        end if
      end function I0

       
      real function I1(x)
        implicit none
        real, intent(in) :: x
        real :: t

        t=x/3.75
        if (x .le. 3.75) then
           I1=(0.5+0.87890594*t**2.+0.51498869*t**4.+0.15084934*t**6+0.02658733*t**8+0.00301532*t**10+0.00032411*t**12)*x
        else
           I1=(0.39894228-0.03988024/t-0.00362018/(t**2.)+0.00163801/(t**3.)-0.01031555/(t**4.)+0.02282967/(t**5.)&
                  &-0.2895312/(t**6.)+0.01787654/(t**7.)-0.00420059/(t**8.))/(sqrt(x)*exp(-x))
        end if
      end function I1

      real function bessk0(xinput)
        implicit none
        real, intent(in) :: xinput
        real :: x
 
        if (xinput .le. 2) then
           x=0.5*xinput
           bessk0=-log(x)*I0(xinput)-0.57721566+0.42278420*x**2.+0.23069756*x**4.+0.03488590*x**6.+0.00262698*x**8.+&
                  &0.00010750*x**10.+0.00000740*x**12.
        else
           x=2./xinput
           bessk0=(1.25331414-0.07832358*x+0.02189568*x**2.-0.01062446*x**3.+0.00587872*x**4.-0.00251540*x**5.+0.00053208*x**6.)&
                  &/(sqrt(xinput)*exp(xinput))
        end if
      end function bessk0

      real function bessk1(xinput)
        implicit none
        real, intent(in) :: xinput
        real :: x

        if (xinput .le. 2) then
           x=0.5*xinput
           bessk1=(xinput*log(x)*I1(xinput)+1.+0.15443144*x**2.-0.67278579*x**4.-0.18156897*x**6.-0.01919402*x**8.&
                  &-0.00110404*x**10.-0.00004686*x**12.)/xinput
        else
           x=2./xinput
           bessk1=(1.25331414+0.23498619*x-0.03655620*x**2.+0.01504268*x**3.-0.00780353*x**4.+0.00325614*x**5.-0.00068245*x**6.)&
                  &/(sqrt(xinput)*exp(xinput))
        end if
      end function bessk1

      real function bessk2(xinput)
        implicit none
        real, intent(in) :: xinput

        bessk2=bessk0(xinput)+2./xinput*bessk1(xinput)

      end function bessk2

      
      function poisson(L)
      implicit none
      integer :: poisson
      real(8) :: ex, r, p
      real(8), intent(in) :: L

      poisson=-1
      p=1.d0
      ex=exp(-L)
      do while(.true.)
         poisson=poisson+1
         call random_number(r)
         p=p*r
         if(p .lt. ex) exit
      end do
      return 
      end function poisson

end module
