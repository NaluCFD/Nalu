      subroutine quadAreaByTriangleFacets( areacoords, area )
c     This subroutine computes the area of a quadrilateral surface by
c     decomposing it into four triangles and summing the contribution
c     of each triangle.
c
c     3D coordinates of each vertex of the face
c     =========================================
c     double areacoords(n_v3d,4)
c
c     area output
c     ===========
c     double area(3)
c
      implicit none
      
      double precision areacoords, area
      double precision r1, r2, xmid
      double precision one4th, half

      dimension areacoords(3,4)
      dimension area(3)

      dimension r1(3), r2(3)
      dimension xmid(3)

      integer triangularFacetTable(3,4)
      
      integer k, itriangle, ntriangles, iq
      
      ! this table defines the triangles composing the face with
      ! the convention that all cross products will point the
      ! same direction by convention
      data triangularFacetTable /
     .  5, 1, 2, 5, 2, 3, 5, 3, 4, 5, 4, 1 /

      one4th = 1.d0/4.d0
      half = 1.d0/2.d0

      do k = 1,3
        xmid(k) = one4th*( areacoords(k,1) + areacoords(k,2)
     .    + areacoords(k,3) + areacoords(k,4) )
      end do

      ntriangles = 4

      do k = 1,3
        area(k) = 0.0d0
        r2(k) = areacoords(k,1) - xmid(k)
      end do

      ! loop over triangles
      do itriangle = 1, ntriangles
        iq = triangularFacetTable(3,itriangle)
        ! construct vectors with common beginning point
        r1 = r2
        do k = 1,3
          r2(k) = areacoords(k,iq)-xmid(k)
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
        area(k) = half*area(k)
      end do

      end
