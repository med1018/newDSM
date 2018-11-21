module vec_sph_harmonics

    implicit none

contains

!########################################################################################
    subroutine comp_vec_sph_harmonics(R,S,T,l_max,m_max,theta,phi)

        integer, intent(in) :: l_max, m_max
        double precision, intent(in) :: theta, phi
        complex(kind(0d0)), intent(out), dimension(3,-m_max:m_max,0:l_max) :: R, S, T

        integer :: l, m
        double precision :: P(0:m_max,0:l_max), dP, c
        complex(kind(0d0)) :: z

        double precision, parameter :: threshold_theta = 1d-5
        complex(kind(0d0)), parameter :: im = (0d0, 1d0)
        double precision, parameter :: pi = 3.1415926535897932d0

        R = 0d0
        S = 0d0
        T = 0d0

        if ( theta > threshold_theta ) then

            call assoc_Legendre_plys(P, l_max, m_max, cos(theta) )

            R(1,0,0) = 1d0 / sqrt( 4d0 * pi )

            do l = 1, l_max

                c = sqrt( dble(2*l+1) / 4d0 / pi )

                do m = 0, min(l, m_max)

                    if ( m /= 0 ) c = - c / sqrt( dble( (l-m+1)*(l+m) ) )

                    dP = dble(l) * P(m,l) / tan(theta) - dble(l+m) * P(m,l-1) / sin(theta)

                    R(1,m,l) = c * P(m,l)
                    S(2,m,l) = c * dP / sqrt( dble( l*(l+1) ) )
                    S(3,m,l) = im * dble(m) * R(1,m,l) / sin(theta) / sqrt( dble( l*(l+1) ) )

                enddo

            enddo

        elseif ( theta >= 0 .and. theta <= threshold_theta ) then

            do l = 0, l_max

                c = sqrt( dble(2*l+1) / 4.0d0 / pi )

                R(1,0,l) =   c * (1.0d0 - dble( l*(l+1) ) * theta * theta / 4d0 )
                S(2,0,l) = - c * sqrt( dble( l*(l+1) ) ) * theta / 2d0


                do m = 1, min(l, m_max)

                    c = - c * sqrt( dble( (l-m+1)*(l+m) ) / dble( 2*m ) )

                    R(1,m,l) = c * theta ** m
                    S(2,m,l) = c * ( theta ** (m-1) ) * dble(m) / sqrt( dble( l*(l+1) ) )
                    S(3,m,l) = im * S(2,m,l)

                enddo

            enddo

        endif

        do m = 1, m_max

            z = exp( im * dble(m) * phi )

            R(:,m,:) = R(:,m,:) * z
            S(:,m,:) = S(:,m,:) * z

            R(:,-m,:) = dble( (-1) ** m ) * conjg( R(:,m,:) )
            S(:,-m,:) = dble( (-1) ** m ) * conjg( S(:,m,:) )

        enddo

        T(2,:,:) =   S(3,:,:)
        T(3,:,:) = - S(2,:,:)

    end subroutine

!########################################################################################
    subroutine assoc_Legendre_plys(P,l_max,m_max,x)

        integer, intent(in) :: l_max, m_max
        double precision, intent(in)  :: x
        double precision, intent(out) :: P(0:m_max,0:l_max)

        integer :: l, m
        double precision :: s

        P = 0d0

        P(0,0) = 1d0

        s = sqrt(1d0 - x * x)

        do m = 0, m_max-1

            P(m,m+1)   = dble(2*m+1) * x * P(m,m)
            P(m+1,m+1) = dble(2*m+1) * s * P(m,m)

        enddo

        do m = 0, m_max
        do l = m+1, l_max-1

            P(m,l+1) = ( dble(2*l+1) * x * P(m,l) - dble(l+m) * P(m,l-1) ) / dble(l-m+1)

        enddo
        enddo

    end subroutine
!########################################################################################

end module