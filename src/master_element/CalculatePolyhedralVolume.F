      subroutine polyhedralVolumeByFaces( ncoords, volcoords,
     .  ntriangles, triangleFaceTable, volume )
c
c     This function works by considering a polyhedron to be bonded
c     by a collection of triangular facets. It then uses the gauss
c     divergence theorem to compute the volume as 1/3 the
c     integral of the divergence of the coordinate vector :
c     3V = integral{ (x1,x2,x3) dot A}.
c
c     number of vertices in the polyhedron
c     ====================================
c     integer ncoords
c
c     coordinates of vertices
c     =======================
c     double volcoords(n_v3d, ncoords)
c
c     number of triangular facets composing polyhedron
c     ================================================
c     integer ntrianges
c
c     table defining vertices of each triangle
c     ========================================
c     integer triangleFaceTable(3,ntriangles)
c
c     volume of polyherdon
c     ====================
c     double volume
c
      implicit none
      integer ncoords, ntriangles, triangleFaceTable
      double precision volcoords, volume
      double precision xface
      
      integer itriangle, ip, iq, ir, k

      dimension volcoords(3, ncoords)
      dimension triangleFaceTable(3,ntriangles)
      dimension xface(3)

      ! initialize volume
      volume = 0.0d0

      ! loop over each triangular facet
      do itriangle = 1,ntriangles
        ! c-index ordering is used in the table, so change to fortran
        ip = triangleFaceTable(1,itriangle)+1
        iq = triangleFaceTable(2,itriangle)+1
        ir = triangleFaceTable(3,itriangle)+1
        ! set spatial coordinate of integration point
        do k = 1,3
          xface(k) = volcoords(k,ip) + volcoords(k,iq) + volcoords(k,ir)
        end do
        ! calculate contribution of triangular face to volume
        volume = volume
     .    + xface(1)*( ( volcoords(2,iq)-volcoords(2,ip) )*
     .                 ( volcoords(3,ir)-volcoords(3,ip) )
     .               - ( volcoords(2,ir)-volcoords(2,ip) )*
     .                 ( volcoords(3,iq)-volcoords(3,ip) ) )
     .    - xface(2)*( ( volcoords(1,iq)-volcoords(1,ip) )*
     .                 ( volcoords(3,ir)-volcoords(3,ip) )
     .               - ( volcoords(1,ir)-volcoords(1,ip) )*
     .                 ( volcoords(3,iq)-volcoords(3,ip) ) )
     .    + xface(3)*( ( volcoords(1,iq)-volcoords(1,ip) )*
     .                 ( volcoords(2,ir)-volcoords(2,ip) )
     .               - ( volcoords(1,ir)-volcoords(1,ip) )*
     .                 ( volcoords(2,iq)-volcoords(2,ip) ) )
      end do

      ! apply constants that were factored out for calculation of
      ! the integration point, the area, and the gauss divergence
      ! theorem.
      volume = volume/18.0d0

      end
