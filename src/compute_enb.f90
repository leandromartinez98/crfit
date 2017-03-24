!
! Subroutine compute_enb: Computes the non-bonded energy terms
! 
! L. Martinez, Dez 2, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine compute_enb(x,vdw_energy,elec_energy)
  
  use force_field
  implicit none
  integer :: inb, kbb, ix, iy, iz, jx, jy, jz
  double precision :: x(*), d2, d, p6, p12, vdw_energy, elec_energy

  vdw_energy = 0.d0
  elec_energy = 0.d0

  ! If all non-bonded interactions will be computed (not only backbone ones)

  if ( ffcomp(7) ) then

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
    
      elec_energy = elec_energy + qq(inb)/d

      p6 = ss6(inb) / d2**3
      p12 = p6*p6
      vdw_energy = vdw_energy + epseps(inb)*( p12 - 2.d0*p6 )

    end do

  ! Compute interactions only for backbone atoms

  else

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
    
      elec_energy = elec_energy + qq(inb)/d

      p6 = ss6(inb) / d2**3
      p12 = p6*p6
      vdw_energy = vdw_energy + epseps(inb)*( p12 - 2.d0*p6 )

    end do

  end if

  return
end subroutine compute_enb

