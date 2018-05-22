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
#ifndef _REVOROPT_LAZY_NN_FWD_HPP_
#define _REVOROPT_LAZY_NN_FWD_HPP_

#include "aabox_fwd.hpp"

#include <memory>

namespace Revoropt {

template< int _Dim = 3, typename _Scalar = double >
class LazyNN : public AABox<_Dim, _Scalar> {

  public :

  /* Typedefs */
  enum { Dim = _Dim } ;
  typedef _Scalar Scalar ;
  typedef Eigen::Matrix<Scalar, Dim, 1> Vector ;

  static const unsigned int NO_INDEX ;

  /* Construction */
  LazyNN( const Scalar* points, unsigned int size ) : size_(0), parent_(nullptr), last_node_(nullptr)
  {
    set_sites(points, size) ;
  }

  LazyNN() : size_(0), parent_(nullptr), last_node_(nullptr)
  {
  }

  /* Update */
  void set_sites( const Scalar* points, unsigned int size ) ;

  /* Nearest neighbor query */
  void knnSearch( const Scalar* query, unsigned int k, unsigned int* indices, Scalar* distances ) ;

  private :

  /* points in the node */
  Scalar* points_ ;
  std::vector<Scalar> points_owner_ ;
  unsigned int* indices_ ;
  std::vector<unsigned int> indices_owner_ ;
  unsigned int size_ ;

  /* partitionning data */
  //sorting along dimensions
  unsigned int* sorted_indices_[Dim] ;
  std::vector<unsigned int> sorted_indices_owner[Dim] ;
  bool dim_sorted_[Dim] ;
  //buffer for cut optimization
  Scalar* buffer_ ;
  std::vector<Scalar> buffer_owner_ ;

  /* child nodes if any */
  std::unique_ptr< LazyNN<Dim, Scalar> > child_[2] ;
  LazyNN<Dim, Scalar>* parent_ ;

  /* Reset bbox to the setup points */
  void reset() { this->init_to(points_, size_) ; }

  /* Nearest neighbor query */
  void knn_down( const Scalar* query, unsigned int k, unsigned int* indices, Scalar* distances ) ;
  void knn_up( const Scalar* query, unsigned int k, unsigned int* indices, Scalar* distances ) ;
  LazyNN<Dim, Scalar>* last_node_ ;

  /* Point query */
  void handle_pt( const Scalar* query, unsigned int pt_index, unsigned int k, unsigned int* indices, Scalar* distances ) ;

  /* Splitting */
  void split_at_query( const Scalar* query ) ;
  void init_from_parent( LazyNN<Dim, Scalar>* parent, unsigned int offset, unsigned int size ) ;

  /* Partitionning along a given axis */
  void partition( unsigned int dim ) ;
  void swap_pt( unsigned int i, unsigned int j ) ;
  Scalar swap_buffer_[Dim] ;

  /* splitting characteristics */
  unsigned int split_axis_ ;
  Scalar split_pos_ ;

} ;

template< int Dim, typename Scalar >
const unsigned int LazyNN<Dim, Scalar>::NO_INDEX = (unsigned int) -1 ;

} //end of namespace Revoropt

#endif
