/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Matthias Troyer <troyer@comp-phys.org>             *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/ngs/mcobservable.hpp>

#include <alps/ngs/boost_python.hpp>
#include <alps/python/numpy_import.hpp>

#include <alps/alea/detailedbinning.h>

#include <boost/make_shared.hpp>

#include <valarray>

namespace alps {
    namespace detail {

        void observable_append(alps::mcobservable & self, boost::python::object const & data) {
            import_numpy();
            if (false);
            #define NGS_PYTHON_HDF5_CHECK_SCALAR(N)                                                                                                 \
                else if (std::string(data.ptr()->ob_type->tp_name) == N)                                                                            \
                    self << boost::python::extract< double >(data)();
            NGS_PYTHON_HDF5_CHECK_SCALAR("int")
            NGS_PYTHON_HDF5_CHECK_SCALAR("long")
            NGS_PYTHON_HDF5_CHECK_SCALAR("float")
            NGS_PYTHON_HDF5_CHECK_SCALAR("numpy.float64")
            else if (std::string(data.ptr()->ob_type->tp_name) == "numpy.ndarray" && PyArray_Check(data.ptr())) {
                PyArrayObject * ptr = (PyArrayObject *)data.ptr();
                if (!PyArray_ISNOTSWAPPED(ptr))
                    throw std::runtime_error("numpy array is not native" + ALPS_STACKTRACE);
                else if (!(ptr = PyArray_GETCONTIGUOUS(ptr)))
                    throw std::runtime_error("numpy array cannot be converted to continous array" + ALPS_STACKTRACE);
                self << std::valarray< double >(static_cast< double const *>(PyArray_DATA(ptr)), *PyArray_DIMS(ptr));
                Py_DECREF((PyObject *)ptr);
            } else
                throw std::runtime_error("unsupported type");
        }

        void observable_load(alps::mcobservable & self, alps::hdf5::archive & ar, std::string const & path) {
            std::string current = ar.get_context();
            ar.set_context(path);
            self.load(ar);
            ar.set_context(current);
        }

        alps::mcobservable create_RealObservable_export(std::string name) {
            return alps::mcobservable(boost::make_shared<alps::RealObservable>(name).get());
        }

        alps::mcobservable create_RealVectorObservable_export(std::string name) {
            return alps::mcobservable(boost::make_shared<alps::RealVectorObservable>(name).get());
        }
    }
}


BOOST_PYTHON_MODULE(pyngsobservable_c) {

    boost::python::def("createRealObservable", &alps::detail::create_RealObservable_export);
    boost::python::def("createRealVectorObservable", &alps::detail::create_RealVectorObservable_export);

    boost::python::class_<alps::mcobservable>(
        "observable",
        boost::python::no_init
    )
        .def("append", &alps::detail::observable_append)
        .def("merge", &alps::mcobservable::merge)
        .def("save", &alps::mcobservable::save)
        .def("load", &alps::detail::observable_load)
        .def("addToObservable", &alps::detail::observable_load)
    ;

}

