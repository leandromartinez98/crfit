!
! Subroutine compute_enbls: Computes the non-bonded energy with short-range modification
! 
! L. Martinez, Dez 2, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine compute_enbls(x,etot)
  
  use force_field
  implicit none
  integer :: inb, kbb, ix, iy, iz, jx, jy, jz
  double precision :: etot, elec_energy, vdw_energy, x(*), d2, d, p6, p12

  ! Compute interactions for all atoms 

  if ( ffcomp(7) ) then

    vdw_energy = 0.d0
    elec_energy = 0.d0
    do inb = 1, n_nonbonded

      ix = (ijnonbonded(inb,1)-1)*3 + 1
      iy = ix + 1
      iz = ix + 2

      jx = (ijnonbonded(inb,2)-1)*3 + 1
      jy = jx + 1
      jz = jx + 2

      d2 = (x(ix) - x(jx))**2 + & 
           (x(iy) - x(jy))**2 + &
           (x(iz) - x(jz))**2

      if ( d2 > cutoff2 ) cycle

      d = dsqrt(d2)

      if ( d > cutnb(inb,1) ) then
      
        elec_energy = elec_energy + qq(inb)/d

        p6 = ss6(inb) / d2**3
        p12 = p6*p6
        vdw_energy = vdw_energy + epseps(inb)*( p12 - 2.d0*p6 )

      else

        ! The energy of the short-range modification is added to VDW

        vdw_energy = vdw_energy + cutnb(inb,2)*d + cutnb(inb,3)

      end if

    end do
    etot = vdw_energy + elec_energy

  ! Compute non-bonded interactions only for backbone atoms

  else

    vdw_energy = 0.d0
    elec_energy = 0.d0
    do kbb = 1, n_nbbackbone

      inb = inbbackbone(kbb)

      ix = (ijnonbonded(inb,1)-1)*3 + 1
      iy = ix + 1
      iz = ix + 2

      jx = (ijnonbonded(inb,2)-1)*3 + 1
      jy = jx + 1
      jz = jx + 2

      d2 = (x(ix) - x(jx))**2 + & 
           (x(iy) - x(jy))**2 + &
           (x(iz) - x(jz))**2

      if ( d2 > cutoff2 ) cycle

      d = dsqrt(d2)

      if ( d > cutnb(inb,1) ) then
      
        elec_energy = elec_energy + qq(inb)/d

        p6 = ss6(inb) / d2**3
        p12 = p6*p6
        vdw_energy = vdw_energy + epseps(inb)*( p12 - 2.d0*p6 )

      else

        ! The energy of the short-range modification is added to VDW

        vdw_energy = vdw_energy + cutnb(inb,2)*d + cutnb(inb,3)

      end if

    end do
    etot = vdw_energy + elec_energy

  end if

  return
end subroutine compute_enbls

