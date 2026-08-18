/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Lukas Gamper <gamperl@gmail.com>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/python/numpy_array.hpp>
#include <iostream>

namespace alps {
    namespace python {
        namespace numpy {

           alps::python::numpy::array from_pyobject(boost::python::object const & source)
           {
               #if defined(ALPS_HAVE_BOOST_NUMPY)
               return boost::python::numpy::array(source);
               #else
               return boost::python::numeric::array(source);
               #endif
           }

            void convert(boost::python::object const & source, std::vector<double> & target) {
                import_numpy();
                target.resize(PyArray_Size(source.ptr()));
                PyArrayObject * ptr = (PyArrayObject *)source.ptr();
                memcpy(&target.front(), static_cast<double *>(PyArray_DATA(ptr)), PyArray_ITEMSIZE(ptr) * target.size());
            }

            alps::python::numpy::array convert(double source) {
                return convert(std::vector<double>(1, source));
            }

            std::vector<double> convert(boost::python::object const & source) {
               std::vector<double> target;
               convert(source, target);
               return target;
            }

            alps::python::numpy::array convert(std::vector<double> const & source) {
                import_numpy();
                npy_intp size = source.size();
                boost::python::object obj(boost::python::handle<>(PyArray_SimpleNew(1, &size, NPY_DOUBLE)));
                void * ptr = PyArray_DATA((PyArrayObject*) obj.ptr());
                memcpy(ptr, &source.front(), PyArray_ITEMSIZE((PyArrayObject*) obj.ptr()) * size);
                return boost::python::extract<alps::python::numpy::array>(obj);
            }

            alps::python::numpy::array convert(std::vector<std::vector<double> > const & source) {
                import_numpy();
                npy_intp size[2] = {static_cast<npy_intp>(source.size()), static_cast<npy_intp>(source[0].size()) };
                boost::python::object obj(boost::python::handle<>(PyArray_SimpleNew(2, size, NPY_DOUBLE)));
                void * ptr = PyArray_DATA((PyArrayObject*) obj.ptr());
                for (std::size_t i = 0; i < source.size(); ++i)
                    memcpy(static_cast<double *>(ptr) + i * size[1], &source[i].front(), PyArray_ITEMSIZE((PyArrayObject*) obj.ptr()) * size[1]);
                return boost::python::extract<alps::python::numpy::array>(obj);
            }

            alps::python::numpy::array convert(std::vector<std::vector<std::vector<double> > > const & source) {
                import_numpy();
                npy_intp size[3] = {
                      static_cast<npy_intp>(source.size())
                    , static_cast<npy_intp>(source[0].size())
                    , static_cast<npy_intp>(source[0][0].size())
                };
                boost::python::object obj(boost::python::handle<>(PyArray_SimpleNew(3, size, NPY_DOUBLE)));
                void * ptr = PyArray_DATA((PyArrayObject*) obj.ptr());
                for (std::size_t i = 0; i < source.size(); ++i)
                    for (std::size_t j = 0; j < source[i].size(); ++j)
                        memcpy(static_cast<double *>(ptr) + i * size[1] * size[2] + j * size[2], &source[i][j].front(), PyArray_ITEMSIZE((PyArrayObject*) obj.ptr()) * size[2]);
                return boost::python::extract<alps::python::numpy::array>(obj);
            }
        }
    }
}
