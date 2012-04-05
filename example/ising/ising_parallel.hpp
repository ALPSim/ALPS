/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * This software is part of the ALPS libraries, published under the ALPS           *
 * Library License; you can use, redistribute it and/or modify it under            *
 * the terms of the license, either version 1 or (at your option) any later        *
 * version.                                                                        *
 *                                                                                 *
 * You should have received a copy of the ALPS Library License along with          *
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       *
 * available from http://alps.comp-phys.org/.                                      *
 *                                                                                 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       *
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       *
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs.hpp>

#include <boost/lambda/lambda.hpp>

// a vector with 'tabs' on both sides: elements at i = -1 as well as
// i = size() are also accessible
template<typename T>
class tabbed_vector {
public:
  typedef T value_type;
  typedef typename std::vector<value_type>::size_type size_type;

  explicit tabbed_vector(size_type n = 0) : vector_(n + 2) {}
  explicit tabbed_vector(size_type n, value_type x) : vector_(n + 2, x) {}
  void resize(size_type n, value_type x = value_type()) { vector_.resize(n + 2, x); }

  size_type size() const { return vector_.size() - 2; }
  value_type const& operator[](int i) const { return vector_[i+1]; }
  value_type& operator[](int i) { return vector_[i+1]; }

private:
  std::vector<value_type> vector_;
};


// rename Impl to Base or similar
template<typename Impl> class parallel_ising_sim : public Impl {
    public:

        parallel_ising_sim(typename Impl::parameters_type const & params, std::size_t seed_offset) 
        {
            throw std::runtime_error("No communicator passed" + ALPS_STACKTRACE);
        }
  
        #ifdef ALPS_HAVE_MPI
// replace: template <class Arg> .... , Arg comm)
            parallel_ising_sim(typename Impl::parameters_type const & params, boost::mpi::communicator c)
                : Impl(params, c)
                , comm(c)
                , length(params["L"])
                , local_length(comm.rank() == 0 ? length - (comm.size() - 1) * (length / comm.size()) : length / comm.size())
                , thermalization_sweeps(int(params["THERMALIZATION"]))
                , total_sweeps(int(params["SWEEPS"]))
                , beta(1. / double(params["T"]))
                , spins(local_length)
            {
                create();
                std::cerr << total_sweeps << ' ' << thermalization_sweeps << std::endl;
            }

        #endif
// update()
        void do_update() {
            for (int i = 0; i < local_length; ++i) {
              double diff = 4 * (spins[i-1] ^ spins[i] + spins[i] ^ spins[i+1]) - 4;
              if (Impl::random() < 0.5 * (1 + std::tanh(-0.5 * beta * diff))) spins[i] ^= 1;
              if (i == 0) copy2left();
              if (i == local_length - 1) copy2right();
            }
        };

  // measure() ... use verbs
        void do_measurements() {
            sweeps++;
            std::cerr << comm.rank() << ' ' << sweeps << std::endl;
            if (sweeps > thermalization_sweeps) {
              double my_energy = 0;
              double my_mag = 0;
              for (int i = 0; i < local_length; ++i) {
                my_energy -= 2 * (spins[i] ^ spins[i+1]) - 1;
                my_mag += (2 * spins[i] - 1);
              }
              if (comm.rank() == 0) {
                double energy, mag;
                reduce(comm, my_energy, energy, std::plus<double>(), 0);
                reduce(comm, my_mag, mag, std::plus<double>(), 0);
                energy /= length;
                mag /= length;
                Impl::measurements["Energy"] << energy;
                Impl::measurements["Energy^2"] << energy * energy;
                Impl::measurements["Magnetization"] << mag;
                Impl::measurements["Magnetization^2"] << mag * mag;
                Impl::measurements["Magnetization^4"] << mag * mag * mag * mag;
              } else {
                boost::mpi::reduce(comm, my_energy, std::plus<double>(), 0);
                reduce(comm, my_mag, std::plus<double>(), 0);
              }
            }
        };

        double fraction_completed() const {
            return (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps));
        }

    private:
        // call it init()
        void create() {
            for(int i = 0; i < local_length; ++i)
                spins[i] = (Impl::random() < 0.5 ? 1 : 0);
            copy2right();
            copy2left();
            if (comm.rank() == 0) {
              Impl::measurements << alps::ngs::RealObservable("Energy")
                                 << alps::ngs::RealObservable("Energy^2")
                                 << alps::ngs::RealObservable("Magnetization")
                                 << alps::ngs::RealObservable("Magnetization^2")
                                 << alps::ngs::RealObservable("Magnetization^4")
              ;
            }
        }

        void copy2right() {
          if (comm.size() == 1) {
            spins[-1] = spins[local_length-1];
          } else {
            comm.send((comm.rank() + 1) % comm.size(), 0, spins[local_length-1]);
            comm.recv((comm.rank() + comm.size() - 1) % comm.size(), 0, spins[-1]);
          }
        }

        void copy2left() {
          if (comm.size() == 1) {
            spins[local_length] = spins[0];
          } else {
            comm.send((comm.rank() + comm.size() - 1) % comm.size(), 0, spins[0]);
            comm.recv((comm.rank() + 1) % comm.size(), 0, spins[local_length]);
          }
        }

        boost::mpi::communicator comm;
        int length, local_length;
        int sweeps;
        int thermalization_sweeps;
        int total_sweeps;
        double beta;
        tabbed_vector<int> spins;
};
