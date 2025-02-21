/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Institute for Theoretical Physics, ETH Zurich
 *               2013-2016 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <boost/filesystem.hpp>

#include <iterator>
#include <iostream>


using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/block_matrix/detail/alps.hpp"
typedef alps::numeric::matrix<double> matrix;

#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/models/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/mp_tensors/mps.h"

/// use NU1 symmetry
typedef NU1 grp;


std::string result_archive;


struct TestConfig {
    TestConfig()
    {
        result_archive = boost::filesystem::unique_path().native();
        result_archive += ".h5";
    }
    ~TestConfig()
    {
        if (boost::filesystem::exists(result_archive))
            boost::filesystem::remove(result_archive);
    }
};

BOOST_GLOBAL_FIXTURE( TestConfig );


struct ArgsConfig {
    std::string reference_archive;
    std::string custom_lattice_library;
    ArgsConfig()
    {
        BOOST_REQUIRE_MESSAGE ( boost::unit_test::framework::master_test_suite().argc == 3, "You miss some arguments: <custom_lattice_library> <reference_archive>"  );
        custom_lattice_library = boost::unit_test::framework::master_test_suite().argv[1];
        reference_archive = boost::unit_test::framework::master_test_suite().argv[2];
    }
};





void check_result_archive(std::string const& result_archive, std::string const& reference_archive, std::string const& rpath)
{
    alps::hdf5::archive ref_ar(reference_archive);
    alps::hdf5::archive res_ar(result_archive);

    std::vector<std::string> children = ref_ar.list_children(rpath);
    for (int i=0; i<children.size(); ++i) {
        if (!res_ar.is_group(rpath+"/"+children[i])) {
            BOOST_ERROR(children[i]+" not found in results.");
            continue;
        }
        
        
        if (ref_ar.dimensions(rpath+"/"+children[i]+"/mean/value") == 1) {
            std::vector<matrix::value_type> ref_meas, res_meas;
            ref_ar[rpath+"/"+children[i]+"/mean/value"] >> ref_meas;
            res_ar[rpath+"/"+children[i]+"/mean/value"] >> res_meas;
            BOOST_CHECK_EQUAL(ref_meas.size(), res_meas.size());
            if (ref_meas.size() != res_meas.size()) continue;
            
            for (int k=0; k<ref_meas.size(); ++k) {
                BOOST_CHECK_CLOSE(ref_meas[k], res_meas[k], 1e-9);
            }
        } else {
            std::vector<std::vector<matrix::value_type> > ref_meas, res_meas;
            ref_ar[rpath+"/"+children[i]+"/mean/value"] >> ref_meas;
            res_ar[rpath+"/"+children[i]+"/mean/value"] >> res_meas;
            BOOST_CHECK_EQUAL(ref_meas[0].size(), res_meas[0].size());
            if (ref_meas[0].size() != res_meas[0].size()) continue;
            
            for (int k=0; k<ref_meas[0].size(); ++k) {
                BOOST_CHECK_CLOSE(ref_meas[0][k], res_meas[0][k], 1e-9);
            }
        }
        BOOST_TEST_CHECKPOINT("Values of "+children[i]+": Finished");
        
        if (ref_ar.is_data(rpath+"/"+children[i]+"/labels")) {
            std::vector<std::string> ref_labels, res_labels;
            ref_ar[rpath+"/"+children[i]+"/labels"] >> ref_labels;
            res_ar[rpath+"/"+children[i]+"/labels"] >> res_labels;
            for (int k=0; k<ref_labels.size(); ++k) {
                BOOST_CHECK_EQUAL(ref_labels[k], res_labels[k]);
            }
            BOOST_TEST_CHECKPOINT("Labels of "+children[i]+": Finished");
        }
        BOOST_TEST_CHECKPOINT("Check of "+children[i]+": Finished");
    }
}

template <class Matrix, class SymmGroup>
void run_test(DmrgParameters & parms, std::string const& result_archive, std::string const& reference_archive) {
    /// Build lattice and model
    Lattice lattice(parms);
    BOOST_TEST_CHECKPOINT("Lattice initialized");
    Model<Matrix, SymmGroup> model(lattice, parms);
    BOOST_TEST_CHECKPOINT("Model initialized");

    typedef boost::ptr_vector<measurement<Matrix, SymmGroup> > measurements_type;
    measurements_type measurements = measurements::parse_and_create(lattice, model, parms);
    BOOST_TEST_CHECKPOINT("Measurements initialized");

    MPS<Matrix, SymmGroup> mps(lattice.size(), *(model.initializer(lattice, parms)));
    BOOST_TEST_CHECKPOINT("MPS initialized");

    std::string const & rpath = boost::unit_test::framework::current_test_case().p_name;
    std::for_each(measurements.begin(), measurements.end(), measure_and_save<Matrix, SymmGroup>(result_archive, rpath, mps));
    BOOST_TEST_CHECKPOINT("Measurements finished");

    check_result_archive(result_archive, reference_archive, rpath);
    BOOST_TEST_CHECKPOINT("Measurements compared");
}


BOOST_FIXTURE_TEST_SUITE( measurements, ArgsConfig )

BOOST_AUTO_TEST_CASE( fermion_hubbard )
{
    DmrgParameters parms;
    parms.set("max_bond_dimension", 2);
    parms.set("LATTICE", "open chain lattice");
    parms.set("L",       10                  );
    parms.set("MODEL",   "fermion Hubbard"   );
    parms.set("CONSERVED_QUANTUMNUMBERS", "Nup,Ndown");
    parms.set("Nup_total",   5);
    parms.set("Ndown_total", 5);
    parms.set("MEASURE_LOCAL[Local density up]",    "n_up");
    parms.set("MEASURE_LOCAL[Local exchance]",      "exchange_xy");
    parms.set("MEASURE_LOCAL[Local hopping]",       "fermion_hop");
    parms.set("MEASURE_LOCAL[Local double occ]",    "double_occupancy");
    parms.set("MEASURE_LOCAL[Local biquadratic]",   "biquadratic");
    parms.set("MEASURE_AVERAGE[Double occ]",        "double_occupancy");
    parms.set("MEASURE_AVERAGE[Hopping]",           "fermion_hop");
    parms.set("MEASURE_AVERAGE[Density]",           "n");
    parms.set("MEASURE_CORRELATIONS[dens-dens]",    "n:n");
    parms.set("MEASURE_HALF_CORRELATIONS[cdag-c]",  "cdag_up:c_up");
    parms.set("MEASURE_LOCAL_AT[Local at meas]",    "cdag_up:n_down:c_up|(0,1,2),(1,2,4)");

    run_test<matrix, grp>(parms, result_archive, reference_archive);
}


BOOST_AUTO_TEST_CASE( kondo )
{
    DmrgParameters parms;
    parms.set("max_bond_dimension", 2);
    parms.set("LATTICE_LIBRARY", custom_lattice_library);
    parms.set("LATTICE", "open Kondo lattice");
    parms.set("L",       5                  );
    parms.set("MODEL",   "Kondo lattice");
    parms.set("MEASURE_LOCAL[Local Sz]",            "Sz");
    parms.set("MEASURE_LOCAL[Local density up]",    "n_up");
    parms.set("MEASURE_LOCAL[Local exchance]",      "exchange_xy");
    parms.set("MEASURE_LOCAL[Local double occ]",    "double_occupancy");
    parms.set("MEASURE_AVERAGE[Double occ]",        "double_occupancy");
    parms.set("MEASURE_AVERAGE[Hopping]",           "fermion_hop");
    parms.set("MEASURE_AVERAGE[Density]",           "n");
    parms.set("MEASURE_CORRELATIONS[dens-dens]",    "n:n");
    parms.set("MEASURE_HALF_CORRELATIONS[cdag-c]",  "cdag_up:c_up");
    parms.set("MEASURE_LOCAL_AT[Local at meas]",    "cdag_up:Sz:c_up|(0,1,2),(1,2,4)");

    run_test<matrix, grp>(parms, result_archive, reference_archive);
}

BOOST_AUTO_TEST_SUITE_END()
