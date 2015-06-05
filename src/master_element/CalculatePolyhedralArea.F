      subroutine polyhedralAreaByFaces( ncoords, areacoords, ntriangles,
     .           triangleFaceTable, area)
c     This subroutine calculates the sum of area vectors composed of
c     triangular facets approximating a curved surface (works exactly)
c     for bilinear interpolation).
c
c     number of vertices composing surface
c     ====================================
c     integer ncoords
c
c     coordinates of vertices
c     =======================
c     double areacoords(n_v3d, ncoords)
c
c     number of triangles composing surface
c     =====================================
c     integer ntriangles
c
c     table defining triangles
c     ========================
c     integer triangleFaceTable
c
c     area output
c     ===========
c     double area(3)
c
      implicit none

      integer ncoords, ntriangles, triangleFaceTable
      integer itriangle, k, ip, iq, ir
      double precision areacoords, area
      double precision r1, r2
     
      dimension areacoords(3, ncoords)
      dimension triangleFaceTable(3,ntriangles)
      dimension area(3)
      dimension r1(3), r2(3)
     
      do k = 1,3
        area(k) = 0.0d0
      end do

      ! loop over triangles
      do itriangle = 1, ntriangles
        ip = triangleFaceTable(1,itriangle)
        iq = triangleFaceTable(2,itriangle)
        ir = triangleFaceTable(3,itriangle)
        ! construct vectors with common beginning point
        do k = 1,3
          r1(k) = areacoords(k,iq)-areacoords(k,ip)
          r2(k) = areacoords(k,ir)-areacoords(k,ip)
        end do
        ! cross product is twice the area vector
        ! triangularFaceTable should be constructed such
        ! that these vectors have the same convention (right hand rule)
        area(1) = area(1) + r1(2)*r2(3) - r2(2)*r1(3)
        area(2) = area(2) + r1(3)*r2(1) - r2(3)*r1(1)
        area(3) = area(3) + r1(1)*r2(2) - r2(1)*r1(2)
      end do

      ! apply the 1/2 that was pulled out
      do k = 1, 3
        area(k) = 0.5d0*area(k)
      end do

      end
