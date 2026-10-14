/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1998-2018 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Adapted from CQMP/Maxent src/maxent_matrix_def.hpp at
// e727650a1f0a427e95efebc409d907b4c84011ec; see ../README.md.

#pragma once

#include <alps/config.h> // needed to set up correct bindings
#include <Eigen/Core>

typedef Eigen::MatrixXd matrix_type;
typedef Eigen::MatrixXcd complex_matrix_type;
typedef Eigen::VectorXd vector_type;
typedef Eigen::MatrixXcd complex_vector_type;
typedef std::pair<vector_type, complex_vector_type> omega_complex_type;

///matrix-vector multiplication. 
inline vector_type maxent_prec_prod(const matrix_type &p, const vector_type &q) {
	return p*q;
}
///matrix-vector multiplication of transpose of matrix. 
inline vector_type maxent_prec_prod_trans(const matrix_type &p, const vector_type &q) {
	return p.transpose()*q;
}
///matrix-matrix multiplication. 
inline matrix_type maxent_prec_prod(const matrix_type &p, const matrix_type &q) {
	return p*q;
}
///matrix-matrix multiplication of transpose(p) with q. 
inline matrix_type maxent_prec_prod_trans(const matrix_type &p, const matrix_type &q) {
	return p.transpose()*q; 
}
