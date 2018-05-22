// @licstart revoropt
// This file is part of Revoropt, a library for the computation and 
// optimization of restricted Voronoi diagrams.
//
// Copyright (C) 2013 Vincent Nivoliers <vincent.nivoliers@univ-lyon1.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// @licend revoropt
#ifndef _REVOROPT_RVD_MEASURE_HPP_
#define _REVOROPT_RVD_MEASURE_HPP_

namespace Revoropt {

template<typename _Triangulation>
class NoAction {
  public :

  /* Types and template data */
  typedef _Triangulation Triangulation ;

  void operator()( unsigned int site, unsigned int triangle,
                   const RVDPolygon<Triangulation>& polygon
                 ) 
  {} ;
} ;

template<typename Triangulation>
class BarycenterAction {
  public :

  /* Types and template data */
  enum { Dim = Triangulation::VertexDim } ;
  typedef typename Triangulation::Scalar Scalar ;
  typedef Eigen::Matrix<Scalar,Dim,1> Vector ;

  /* Construction */
  BarycenterAction(unsigned int nsites) : 
    areas(nsites, 0), 
    barycenters(Dim * nsites, 0) 
  {} ;

  void operator()( unsigned int site, unsigned int triangle,
                   const RVDPolygon<Triangulation>& polygon
                 ) {
    //size of the polygon
    unsigned int size = polygon.size() ;
    if(size<3) return ;

    //triangulate the polygon
    const Vector& base_vertex = polygon[0].vertex ;
    for(unsigned int i = 1; i < polygon.size() - 1; ++i) {
      const Vector& v1 = polygon[ i        ].vertex ;
      const Vector& v2 = polygon[(i+1)%size].vertex ;

      Scalar area = triangle_area<3>(
          base_vertex.data(),
          v1.data(),
          v2.data()
          ) ;

      Eigen::Map<Vector> barycenter(barycenters.data() + Dim*site) ;
      barycenter += area * (base_vertex + v1 + v2) / 3 ;
      areas[site] += area ;
    }
  };

  std::vector<Scalar> areas ;
  std::vector<Scalar> barycenters ;
} ;

} //end of namespace Revoropt

#endif
