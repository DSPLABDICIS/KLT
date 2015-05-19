/*
 * boost_incl.hpp
 *
 *  Created on: Jan 25, 2013
 *      Author: luzdora
 */

#ifndef BOOST_INCL_HPP_
#define BOOST_INCL_HPP_

#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"
#include "boost/numeric/ublas/io.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/symmetric.hpp"

namespace boost_incl{

typedef boost::numeric::ublas::matrix<double> mat;
typedef boost::numeric::ublas::vector<double> vec;
typedef boost::numeric::ublas::bounded_vector<double,3> vec3;
typedef boost::numeric::ublas::symmetric_matrix<double> sym_mat;
typedef boost::numeric::ublas::vector_range<vec> vec_range;
typedef boost::numeric::ublas::matrix_range<mat> mat_range;
typedef boost::numeric::ublas::identity_matrix<double> identity_mat;
typedef boost::numeric::ublas::matrix_vector_range<mat> mvec_range;

}//end namespace

#endif /* BOOST_INCL_HPP_ */
