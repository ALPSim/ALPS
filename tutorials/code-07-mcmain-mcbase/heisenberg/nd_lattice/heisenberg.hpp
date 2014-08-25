/* heienberg.hpp
 * adapted from heisenberg_beispiel/heisenberg.hpp
 */

#ifndef HEISENBERG_HPP
#define HEISENBERG_HPP

#include <alps/mcbase.hpp>
#include <alps/ngs/numeric.hpp>
#include <alps/random/uniform_on_sphere_n.h>
#include <alps/lattice.h>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/array.hpp>

#include <boost/function.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/array.hpp>

#include <vector>
#include <string>

typedef boost::array<double, 3> spintype;

class ALPS_DECL heisenberg_sim : public alps::mcbase {

    public:

        heisenberg_sim(parameters_type const & parms, std::size_t seed_offset = 0);

        virtual void update();
        virtual void measure();
        virtual double fraction_completed() const;

        const spintype random_spin_accept_reject();
        const spintype random_spin();

        using alps::mcbase::save;
        virtual void save(alps::hdf5::archive & ar) const;

        using alps::mcbase::load;
        virtual void load(alps::hdf5::archive & ar);

    private:
        
        int length;
        int sweeps;
        int thermalization_sweeps;
        int total_sweeps;
        double beta;
        alps::uniform_on_sphere_n<3, double, std::vector<double> > random_spin_gen;
        alps::graph_helper<> lattice;
        std::vector<spintype > spins;
};

#endif 
