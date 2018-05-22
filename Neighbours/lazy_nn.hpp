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
#ifndef _REVOROPT_LAZY_NN_HPP_
#define _REVOROPT_LAZY_NN_HPP_

#include "lazy_nn_fwd.hpp"

#include <limits>

namespace Revoropt {

template< int Dim, typename Scalar >
void LazyNN<Dim, Scalar>::set_sites( const Scalar* points, unsigned int size )
{
  //delete hierarchy
  child_[0].reset() ;
  child_[1].reset() ;

  //reset memory if necessary
  if(size != size_) {
    //reset size
    size_ = size ;
    //replace previous memory
    points_owner_.resize(Dim*size) ;
    buffer_owner_.resize(size) ;
    indices_owner_.resize(size) ;
    //get the allocated memory
    points_ = points_owner_.data() ;
    buffer_ = buffer_owner_.data() ;
    indices_ = indices_owner_.data() ;
  }

  //copy the point data
  std::copy(points, points + Dim*size, points_) ;

  //setup indices
  for(unsigned int i = 0; i < size; ++i) {
    indices_[i] = i ;
  }

  //initialize bounding box
  this->init_to(points_, size_) ;
}

template< int Dim, typename Scalar >
void LazyNN<Dim, Scalar>::knnSearch( 
    const Scalar* query, 
    unsigned int k,
    unsigned int* indices, 
    Scalar* distances 
    )
{
  //initialize distances
  for(unsigned int i = 0; i < k; ++i) {
    /* distances[i] = std::numeric_limits<Scalar>::max() ; */
    distances[i] = Scalar(std::numeric_limits<double>::max()) ; // FIXME (maybe): why std::numeric_limits<Scalar>::max() does not work?
  }

  //handle k bigger than size
  k = k > size_ ? size_ : k ;

  LazyNN<Dim, Scalar>* start = this ;

  //use hint
  if(last_node_ != nullptr) { //a previous search was done
    start = last_node_ ;
  } 

  unsigned int init_size = start->size_ > k ? k : start->size_ ;
  for(unsigned int i = 0; i < init_size; ++i) {
    start->handle_pt(query, i, k, indices, distances) ;
  }

  //knn search
  //start->knn_up(query, k, indices, distances) ;
  knn_down(query, k, indices, distances) ;
}

template< int Dim, typename Scalar >
void LazyNN<Dim, Scalar>::knn_up( 
    const Scalar* query, 
    unsigned int k,
    unsigned int* indices, 
    Scalar* distances 
    )
{
  if(!this->parent_) {
    return knn_down(query, k, indices, distances) ;
  }

  Scalar d = this->sq_distance_to(query) ;
  if(d < 0 && -d > distances[k-1]) {
    return knn_down(query, k, indices, distances) ;
  }

  this->parent_->knn_up(query, k, indices, distances) ;
}

template< int Dim, typename Scalar >
void LazyNN<Dim, Scalar>::knn_down( 
    const Scalar* query, 
    unsigned int k,
    unsigned int* indices, 
    Scalar* distances 
    )
{
  //fast return
  if(!size_) return ;
  if(this->sq_distance_to(query) > distances[k-1]) return ;

  //children exploration
  if(child_[0]) {
    if(query[split_axis_] < split_pos_) {
      if(child_[0]) child_[0]->knn_down(query, k, indices, distances) ;
      if(child_[1]) child_[1]->knn_down(query, k, indices, distances) ;
      return ;
    } else {
      if(child_[1]) child_[1]->knn_down(query, k, indices, distances) ;
      if(child_[0]) child_[0]->knn_down(query, k, indices, distances) ;
      return ;
    }
  }

  //linear search
  for(unsigned int i = 0; i < size_; ++i) {
    handle_pt(query, i, k, indices, distances) ;
  }

  //split
  if(size_ > 50) { //FIXME benchmark this parameter
    split_at_query(query) ;
  }
}

template< int Dim, typename Scalar >
void LazyNN<Dim, Scalar>::handle_pt( 
    const Scalar* query,
    unsigned int pt_index,
    unsigned int k,
    unsigned int* indices, 
    Scalar* distances 
    )
{
  Eigen::Map<const Vector> q(query) ;
  Eigen::Map<const Vector> p(points_ + Dim * pt_index) ;
  Scalar d = (p-q).squaredNorm() ;
  unsigned int pos = k ;
  while(pos > 0) {
    if(d < distances[pos - 1]) {
      if(pos < k) {
        distances[pos] = distances[pos-1] ;
        indices[pos] = indices[pos-1] ;
      }
      --pos ;
      if(pos == 0) {
        last_node_ = this ;
      }
    } else {
      break ;
    }
  }
  if(pos < k) {
    distances[pos] = d ;
    indices[pos] = indices_[pt_index] ;
  }
}

template< int Dim, typename Scalar >
void LazyNN<Dim, Scalar>::split_at_query( 
    const Scalar* query
    )
{
  //allocate space for new nodes
  child_[0].reset(new LazyNN<Dim, Scalar>) ;
  child_[1].reset(new LazyNN<Dim, Scalar>) ;

  Scalar max_dist = 0 ;
  unsigned int best_d = 0 ;

  /*
  //test the distances of the split boxes wrt. the query on each axis
  for(unsigned int d = 0; d < Dim; ++d) {
    this->split_at_ratio(d, 0.5, *child_[0], *child_[1]) ;
    Scalar dist = child_[0]->sq_distance_to(query) ;
    if(dist > max_dist) {
      max_dist = d ;
      best_d = d ;
      continue ;
    }
    dist = child_[1]->sq_distance_to(query) ;
    if(dist > max_dist) {
      max_dist = dist ;
      best_d = d ;
    }
  }
  */

  /*
  if(parent_) {
    best_d = (parent_->split_axis_ + 1) % Dim ;
  }
  */

  for(unsigned int d = 0; d < Dim; ++d) {
    Scalar dist = this->bounds_[2*d + 1] - this->bounds_[2*d] ;
    if(dist > max_dist) {
      max_dist = dist ;
      best_d = d ;
    }
  }

  //partition
  partition(best_d) ;
}

template< int Dim, typename Scalar >
void LazyNN<Dim, Scalar>::init_from_parent( 
    LazyNN<Dim, Scalar>* parent, 
    unsigned int offset, 
    unsigned int size 
    ) 
{
  parent_ = parent ;
  points_ = parent->points_ + Dim * offset ;
  indices_ = parent->indices_ + offset ;
  //revert_indices_ = parent->revert_indices_ ;
  //nodes_ = parent->nodes_ + offset ;
  size_ = size ;
  //offset_ = parent->offset_ + offset ;
  this->init_to(points_, size_) ;
  //std::fill(nodes_, nodes_ + size, this) ;
}

template< int Dim, typename Scalar >
void LazyNN<Dim, Scalar>::partition( 
    unsigned int d
    )
{
  //struct for easy sorting
  struct PtDimCompare {
    PtDimCompare(const Scalar* pts, unsigned int d) : pts_(pts), d_(d) {} ;
    bool operator() (unsigned int i, unsigned int j) {
      return pts_[Dim*i + d_] < pts_[Dim*j + d_] ;
    }
    const Scalar* pts_ ;
    const unsigned int d_ ;
  } ;

  //split position at midpoint
  Scalar mid = 0.5 * (this->bounds_[2*d] + this->bounds_[2*d + 1]) ;
  split_axis_ = d ;
  split_pos_ = mid ;
  //index of first point above the split
  unsigned int first_after = 0 ;
  for(unsigned int i = 0; i < size_; ++i) {
    //put points below the split are in the first portion of the array
    if(points_[Dim * i + d] < mid) {
      swap_pt(i, first_after) ;
      ++first_after ;
    }
  }
  //setup children
  child_[0]->init_from_parent(this, 0, first_after) ;
  child_[1]->init_from_parent(this, first_after, size_ - first_after) ;
}

template< int Dim, typename Scalar >
void LazyNN<Dim, Scalar>::swap_pt( 
    unsigned int i,
    unsigned int j
    )
{
  //swap coordinates
  std::copy(points_ + i*Dim, points_ + i*Dim + Dim, swap_buffer_) ;
  std::copy(points_ + j*Dim, points_ + j*Dim + Dim, points_ + i*Dim) ;
  std::copy(swap_buffer_, swap_buffer_ + Dim, points_ + j*Dim) ;
  //swap indices
  unsigned int tmp = indices_[i] ;
  indices_[i] = indices_[j] ;
  indices_[j] = tmp ;
  ////setup revert_indices
  //revert_indices_[indices_[j]] = offset_ + j ;
  //revert_indices_[indices_[i]] = offset_ + i ;
}


} //end of namespace Revoropt

#endif
