/*****************************************************************************
*
* Copyright (C) 2005 -2015  Emanuel Gull <gullc@comp-phys.org>
*
*****************************************************************************/

#pragma once
#include "alps/numeric/matrix.hpp"

///enum for green's function shape
enum shape_t {diagonal, blockdiagonal, nondiagonal};

///matrix inversion for non-diagonal matrices
void matrix_inversion(alps::numeric::matrix<std::complex<double> >&A, alps::numeric::matrix<std::complex<double> >&A_inv);
