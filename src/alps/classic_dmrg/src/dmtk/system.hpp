// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/system.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers; the header includes <map>,
// which it uses (the engine tests include it without a prior <map>).

#ifndef DMTK_SYSTEM_H_
#define DMTK_SYSTEM_H_

//*******************************************************************
// System Class: Defines system properties and performs DMRG loops
// This version of the system class works only on symmetric systems
//*******************************************************************
 
#include "dmtk/output.hpp"
#include <array>
#include <functional>
#include <iostream>
#include <iosfwd>
#include <list>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <iomanip>
#include <cmath>
#include "dmtk/enums.hpp"
#include "dmtk/vector.hpp"
#include "dmtk/matrix.hpp"
#include "dmtk/state.hpp"
#include "dmtk/operators.hpp"
#include "dmtk/block.hpp"
#include "dmtk/lattice.hpp"
#include "dmtk/hami.hpp"
#include "dmtk/support.hpp"
#include "dmtk/util.hpp"
#include "dmtk/globals.hpp"
#include "dmtk/product.hpp"
#include "dmtk/signals.hpp"
#ifdef WITH_PTHREADS
#endif // WITH_PTHREADS


#include <filesystem>

namespace dmtk
{

using namespace std;

enum{
 SYSTEM_SIGNAL_GS,
 SYSTEM_SIGNAL_TRUNCATE,
 SYSTEM_SIGNAL_BUILD_DM,
 SYSTEM_SIGNAL_DM_READY,
 SYSTEM_SIGNAL_ROTATE,
 SYSTEM_SIGNAL_START_ITER,
 SYSTEM_SIGNAL_END_ITER,
 SYSTEM_SIGNAL_MEASURE,
};

#define COUT_PRODUCT_DEFAULT(term) \
   dmtk::out() << "PRODUCT " << term.description() << '\n';

#define COUT_PRODUCT(term) \
 if(calc_hc) \
   dmtk::out() << "PRODUCT " << term.description() << " + h.c." << '\n'; \
 else  \
   dmtk::out() << "PRODUCT " << term.description() << '\n';

template<class T> 
VectorState<T> product(const System<T>& ss, const VectorState<T>& vv);

template<class T> 
VectorState<T> product_default(const System<T>& ss, const VectorState<T>& vv, const Hami<T>* hami = nullptr, bool only_local = false);

template<class T>
VectorState<T> product_composite(System<T>& ss, const VectorState<T>& vv, int mask_hc = (MASK_PRODUCT_DEFAULT|MASK_PRODUCT_HC));

template<class T>
class AuxTerm{
  public:
    Term<T> t;
    Term<T> cterm[4]; // composite term
    BasicOp<T> top[4]; // reference operator
    const BasicOp<T> *op[4]; // original operator
    const BasicOp<T> *ref_op[4]; // complementary operator
    BasicOp<T> sum_op; // complementary operator
    size_t pos[4]; // original position in the term
    size_t b1, b2, b3;
    int mask;
    bool done;
    int nblocks;

    AuxTerm(): b1(BLOCK_NONE), b2(BLOCK_NONE), b3(BLOCK_NONE), mask(0), done(false), nblocks(0) { op[0] = op[1] = op[2] = op[3] = NULL; ref_op[0] = ref_op[1] = ref_op[2] = ref_op[3] = NULL; }
};

template<class T>
void product_term(const System<T>& ss, const AuxTerm<T> &auxt, const Term<T> &t, const VectorState<T> &vv, VectorState<T> &res);
template<class T>
void init_terms_composite(System<T>& ss, const VectorState<T>& vv);
template<class T>
void init_term_composite(System<T>& ss, const AuxTerm<T> &auxt, const Term<T> &t, const VectorState<T> &vv, VectorState<T> &res);

class FileList 
{
  private:
    std::map<std::string,std::string> _tmp_filenames;
    std::string last_filename;
    std::filesystem::path temp_dir;
  public:
    // Static construction must be silent; an invocation validates and reports
    // its selected scratch directory through set_temp_dir.
    FileList() : temp_dir(".") {}
    FileList(const char *dir) { set_temp_dir(dir); }

    void swap(FileList& other) noexcept
      {
         _tmp_filenames.swap(other._tmp_filenames);
         last_filename.swap(other.last_filename);
         temp_dir.swap(other.temp_dir);
      }

    void set_temp_dir(const char *dir) 
      { 
         std::string aux = std::string(dir);
         if(!std::filesystem::is_directory(std::filesystem::path(aux))) {
             dmtk::err() << "*** ERROR: ALPS DMRG could not open directory for temporary files. Create the directory " + aux + " or choose a different path.\n";
                         throw std::runtime_error("*** ERROR: ALPS DMRG could not open directory for temporary files. Create the directory " + aux + " or choose a different path.");
         }
         // Diagnostic; routed to stderr so callers consuming the
         // dmrg binary's stdout (or a alps binding that captures
         // stdout for machine-readable output like `alps doctor
         // --json`) aren't contaminated by it.
         dmtk::err() << "ALPS DMRG temporary files will be written to " << aux << std::endl;
         temp_dir = std::filesystem::path(aux);
      }

    const char * get_filename(const char *input) 
      {
         std::string filename(input);
         std::map<std::string,std::string>::iterator old_name = _tmp_filenames.find(filename);
         if(old_name != _tmp_filenames.end()) return (old_name->second).c_str();

         filename = filename.substr(0,filename.find_first_of('.'));
         this->last_filename = dmtk::temporary_filename((temp_dir / filename).string());
         _tmp_filenames[std::string(input)] = this->last_filename;
         dmtk::err() << "Creating temp file " << this->last_filename << std::endl;
         return this->last_filename.c_str();
      }
};

// `inline` makes this a single shared definition across every TU that
// includes system.hpp. Without it, linking the impl-lib (factory.C) and
// any binding TU (apps/dmrg/dmrg/factory.cpp → dmrg.hpp → system.hpp, plus
// the same chain in alps's dmrg_bindings.cpp) produced two
// definitions of tmp_files — both ran the FileList() ctor at
// process start (printing "ALPS DMRG temporary files will be written
// to ." twice), and both ran the destructor at exit, double-freeing
// the same temp files and SIGABRT'ing the process.
inline FileList tmp_files;

template<class T>
class System
{
  protected:
    using B = typename dmtk::Block<T>;

    Lattice _lattice;

    dmtk::Stopwatch _timer;

    bool _in_warmup;
    size_t dir;
    int _sweep;
    int iter;

    bool _use_hc;
    int _grand_canonical;

    int m;
    int m1;
    int m2;
    int m3;
    int m4;

    double _error;  
    double _truncation_error;
    double _entropy;
    bool _use_error;
    int _error_max_size;
    bool _target_subspaces; // build density matrix using all subspaces
    int _nsub; // number states per subspace

    int _verbose;
    int _calc_gap;
    bool _project;
    bool _use_k;
    bool _use_seed;
    bool _use_basic_seed;
    bool _use_composite;
    bool _grow_symmetric;
    bool _grow_outward;
    bool _custom_qns;
    bool _store_products;
    bool _measure_symmetric;
    bool _apply_hami;
    bool _apply_extern;
    bool _rotate_terms;
    bool _full_sweep; // run sweep from 1 to ls-3/ls-4
    bool _use_coef_tol; 
    double _coef_tol;
    Matrix<size_t> _sweeps;
    size_t _numsweeps;

    char _name[255];

    double _lanczos_tol;
    int _lanczos_maxiter;
    bool _lanczos_certify; // measured states: Ritz residuals down to the tolerance

    // Changed here: final_sweep's trims of what no measurement reads (see
    // there). rotate() carries the measurement operators into the new block
    // only while _rotate_corr is set; read_block refuses a block file that
    // final_sweep removed and nothing has written since.
    bool _rotate_corr = true;
    mutable std::set<std::string> _discarded_blocks;

    virtual void init_iteration(const B&b1, const B&b2, const B&b3, const B&b4, bool use_seed = false, bool create_composite = true);
    void rotate_hami(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *this_hami = nullptr);
    void rotate_terms(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *this_hami = nullptr);
    void rotate_corr(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *this_hami = nullptr);

    typename std::list<Signal<System<T> > > handlers;
    std::vector<std::string> _signal_id;

    void init_signals()
    {
        _signal_id.push_back(std::string("SYSTEM_SIGNAL_GS"));
        _signal_id.push_back(std::string("SYSTEM_SIGNAL_TRUNCATE"));
        _signal_id.push_back(std::string("SYSTEM_SIGNAL_BUILD_DM"));
        _signal_id.push_back(std::string("SYSTEM_SIGNAL_DM_READY"));
        _signal_id.push_back(std::string("SYSTEM_SIGNAL_ROTATE"));
        _signal_id.push_back(std::string("SYSTEM_SIGNAL_START_ITER"));
        _signal_id.push_back(std::string("SYSTEM_SIGNAL_END_ITER"));
        _signal_id.push_back(std::string("SYSTEM_SIGNAL_MEASURE"));
        _signal_id.push_back(std::string("SYSTEM_SIGNAL_END_SWEEP"));
    }
    int get_signal_id(std::string signal_id){
         int n = 0;
         std::vector<std::string>::iterator iter = _signal_id.begin();
         while (iter!=_signal_id.end()) {
           if (*iter==signal_id) break;
           ++iter;
           ++n;
         }
         if(iter == _signal_id.end()) {
            dmtk::err() << "ERROR: signal_connect : signal " << signal_id << " not found \n";
            return -1;
         }
         return n;
    }

  public:

    std::vector<ProductTerm<T> > product_terms;

    B rightblock;
    B leftblock;
    B site1;
    B site2;
    B newblock;

    const B* _b1;
    const B* _b2;
    const B* _b3;
    const B* _b4;

    QN qnt;
    QN qn;
    double precision;

    size_t _ntargets;
    Vector<VectorState<T> > _target;   
    Vector<double> _target_weight;
    Vector<VectorState<T> *> _project_states;

    size_t _nstates;  // excited states for diagonalization
    Vector<VectorState<T> > _state; // gs and excited states  
    Vector<VectorState<T> > _propagate_state; // states that we want to transform as we sweep (besides the ground state)

    BMatrix<T> rho;
    Basis rho_basis;
    VectorState<T> gs;
    VectorState<T> seed;
    Vector<double> energy;

    Hami<T> h;
    std::list<Hami<T> > hs;
    H<T> h12;
    H<T> h23;
    H<T> h34;
    std::list<AuxTerm<T> > aux_terms;

    std::list<BasicOp<T> > dm_ops;
    Hami<T> ops;
    Hami<T> corr;
    Hami<T> control_ops; // single operators

//  Constructor
  
    System()
    : _in_warmup(true)
    , _sweep(1)
    , iter(1)
    , _use_hc(false)
    , _grand_canonical(QN::default_mask())
    , m(10)
    , m1(2)
    , m2(2)
    , m3(2)
    , m4(2)
    , _use_error(false)
    , _error_max_size(-1)
    , _target_subspaces(false)
    , _nsub(1)
    , _verbose(0)
    , _calc_gap(0)
    , _project(false)
    , _use_k(false)
    , _use_seed(true)
    , _use_basic_seed(false)
    , _use_composite(true)
    , _grow_symmetric(true)
    , _grow_outward(false)
    , _custom_qns(true)
    , _store_products(true)
    , _measure_symmetric(false)
    , _apply_hami(true)
    , _apply_extern(true)
    , _rotate_terms(true)
    , _full_sweep(false)
    , _use_coef_tol(false)
    , _coef_tol(1.e-5) 
    , _numsweeps(0)
    , _lanczos_tol(1.e-7)
    , _lanczos_maxiter(-1)
    , _lanczos_certify(false)
    , qnt(0)
    , _ntargets(1)
    , _nstates(1)
    { _sweeps.resize(2,10); _target.resize(1); _target_weight.resize(1); set_name(""); init_signals(); };

    System(const Hami<T> &_h, const Lattice& lattice, const char *the_name)
     : _lattice(lattice)
     , _in_warmup(true)
     , _sweep(1)
     , iter(1)
     , _use_hc(false)
     , _grand_canonical(QN::default_mask())
     , _use_error(false)
     , _error_max_size(-1)
     , _target_subspaces(false)
     , _nsub(1)
     , _verbose(0)
     , _calc_gap(0)
     , _project(false)
     , _use_k(false)
     , _use_seed(true)
     , _use_basic_seed(false)
     , _use_composite(true)
     , _grow_symmetric(true)
     , _grow_outward(false)
     , _custom_qns(true)
     , _store_products(true)
     , _measure_symmetric(false)
     , _apply_hami(true)
     , _apply_extern(true)
     , _rotate_terms(true)
     , _full_sweep(false)
     , _use_coef_tol(false)
     , _coef_tol(1.e-5) 
     , _numsweeps(0)
     , _lanczos_tol(1.e-7)
     , _lanczos_maxiter(-1)
     , _lanczos_certify(false)
     , qnt(0)
     , _ntargets(1)
     , _nstates(1)
     , h(_h)
      { set_name(the_name); m1 = m2 = m3 = m4 = _h.get_site(0).dim(); m = m1 * m2; _sweeps.resize(2,10); _target.resize(1); _target_weight.resize(1); _target_weight[0] = double(1); init_signals(); }

    virtual ~System() {}
//  Methods

    void start() 
      { 
        h.reorder_terms();
        write_block(h.get_site(0),1,LEFT); 
        write_block(h.get_site(_lattice.size()-1),1,RIGHT); 
      }

    void start(size_t numsweeps, const Matrix<size_t> &nstates)
      {
        _numsweeps = numsweeps;
        _sweeps = nstates;
        start();
      }

    void resume(size_t resume_sweep, size_t resume_dir, size_t resume_iter)
      {
        for(int _sweep = resume_sweep; _sweep <= _numsweeps; _sweep++){
          dmtk::out() << "NEW SWEEP " << _sweep << " " << _sweeps(0,_sweep-1) << " " << _sweeps(1,_sweep-1) << '\n';
          sweep(_sweeps(0,_sweep-1),_sweeps(1,_sweep-1),static_cast<size_t>(resume_dir),resume_iter);
        }
      }

    void resume_default()
      {
        read_status();

        int resume_sweep = _sweep;
        size_t resume_dir = dir;
        int resume_iter = iter;

        resume(resume_sweep, resume_dir, resume_iter);
      }

    void run(int nwarmup = 20)
      {
        if(nwarmup > 0) warmup_loop(nwarmup);
        resume(1, RIGHT2LEFT, 1);
      }

    void set_error(double err, int error_max_size = -1) { _error = err; _error_max_size = error_max_size; _use_error = true; }

    void set_truncation(size_t _m) { m = _m; }
    bool in_warmup() const { return _in_warmup; }

    void set_name(const char *the_name) 
      { snprintf(_name, sizeof(_name), "%s", the_name); }
    const char* name() const { return _name; }

    void create_interactions(bool create_composite = true);
    virtual void diagonalize(bool use_seed = true);
    virtual void build_target_states(int pos = 0);
    virtual void truncate(int block, int new_size);
    virtual void rotate(int block, Block<T>& b); 
    virtual void warmup_loop(size_t t);
    virtual void warmup_loop(size_t t, const Vector<QN> &qns);
    virtual void sweep(size_t t1, size_t t2, size_t dir=RIGHT2LEFT, int _start=1);
    virtual void init_correlations();
    virtual void final_sweep(size_t t, size_t dir, int _start, bool _rotate);
    virtual void measure(const VectorState<T> *v=nullptr); 
    virtual T measure_operator(const BasicOp<T> &top, const VectorState<T> *v, bool &info);

    size_t get_iter() const { return iter; }
    size_t get_dir() const { return dir; }


    size_t block(int site) const;

    const BasicOp<T>* operator()(const BasicOp<T>& op) const;

    // Written for this tree: the block (0-3) that holds a site and the
    // site's index as that block stores it, the resolution block() and
    // operator()(op) share; -1 when no block holds the site.
    int locate(int site, int &local_site) const;
    // Written for this tree: lookup tables over the four blocks with the
    // match of operator()(op) (block by site, then name, site and internal
    // site, first match in list order), built once per call of hint(),
    // rotate_terms(), build_composite_operators(), rotate_corr() or
    // measure() instead of a list scan per lookup. Valid only while the four
    // blocks are unchanged, i.e. within that call.
    using OpKey = std::tuple<std::string, int, int>;
    struct OpIndex {
      std::array<std::map<OpKey, const BasicOp<T>*>, 4> ops;
      std::array<bool, 4> built{};
    };
    const BasicOp<T>* indexed_op(const BasicOp<T>& op, OpIndex &index) const;
    const Block<T>& operator[](size_t pos) const 
      { 
        const Block<T> *b;
        switch(pos){
          case BLOCK1: b = _b1; break;
          case BLOCK2: b = _b2; break;
          case BLOCK3: b = _b3; break;
          case BLOCK4: b = _b4; break;
        }
        return *b;
      }

    size_t size() const
      { return _b1->lattice().size()+ _b2->lattice().size()+
               _b3->lattice().size()+ _b4->lattice().size(); }
//    int site1() const 
//      { if(dir == LEFT2RIGHT) return iter; else return h.lattice().size()-iter-2; }
//    int site2() const { return site1()+1; }

    const Lattice& lattice() const { return _lattice; } 

    System<T>& set_ntargets(size_t ntargets)
      { _ntargets = ntargets >= 1 ? ntargets : 1; return *this; } 
    size_t ntargets() { return _ntargets; }

    System<T>& set_target_weights(const Vector<double>& w)
      { _target_weight = w; return *this; }
    

    System<T>& set_store_products(bool b) { _store_products = b; return *this; }
    bool store_products() const { return _store_products; }
    System<T>& set_target_subspaces(bool b, int nsub = 1) { _target_subspaces = b; _nsub = nsub; return *this; }
    bool target_subspaces() const { return _target_subspaces; }
    System<T>& set_verbose(int b) { _verbose = b; return *this; }
    System<T>& set_verbose(bool b) { _verbose = b ? 1 : 0; return *this; }
    int verbose() const { return _verbose; }
    System<T>& set_translations(bool b) { _use_k = b; return *this; }
    bool translations() const { return _use_k; }
    System<T>& set_use_hc(bool b) { _use_hc = b; return *this; }
    bool use_hc() const { return _use_hc; }
    System<T>& set_use_seed(bool b) { _use_seed = b; return *this; }
    bool use_seed() const { return _use_seed; }
    System<T>& set_use_basic_seed(bool b) { _use_basic_seed = b; return *this; }
    bool use_basic_seed() const { return _use_basic_seed; }
    System<T>& set_use_composite(bool b) { _use_composite = b; return *this; }
    bool use_composite() const { return _use_composite; }
    System<T>& set_apply_hami(bool b) { _apply_hami = b; return *this; }
    bool apply_hami() const { return _apply_hami; }
    System<T>& set_apply_extern(bool b) { _apply_extern = b; return *this; }
    bool apply_extern() const { return _apply_extern; }
    System<T>& set_use_coef_tol(bool b, double tol=1.e-5) { _use_coef_tol = b; _coef_tol = tol; return *this; }
    bool use_coef_tol() const { return _use_coef_tol; }
    double coef_tol() const { return _coef_tol; }
    System<T>& set_grow_symmetric(bool b) { _grow_symmetric = b; return *this; }
    bool grow_symmetric() const { return _grow_symmetric; }
    System<T>& set_grow_outward(bool b) { _grow_outward = b; return *this; }
    bool grow_outward() const { return _grow_outward; }
    System<T>& set_measure_symmetric(bool b) { _measure_symmetric = b; return *this; }
    bool measure_symmetric() const { return _measure_symmetric; }
    System<T>& set_grand_canonical(bool b) { _grand_canonical = b ? 0 : QN::default_mask(); return *this; }
    bool grand_canonical() const { return (_grand_canonical == 0); }
    System<T>& set_qn_mask(int m = QN::default_mask()) { _grand_canonical = m; return *this; }
    int qn_mask() const { return _grand_canonical; }
    System<T>& set_full_sweep(bool b) { _full_sweep = b; return *this; }
    bool full_sweep() const { return _full_sweep; }
    void set_rotate_terms(bool b) { _rotate_terms = b; }
    bool rotate_terms() { return _rotate_terms; }

    void hint(int _mask, bool add_local_h = true, const Hami<T> *_h = NULL);
    void build_composite_operators(const Hami<T> *_h = NULL);

    // Changed here: three hooks through which a derived system supplies its
    // own Hamiltonian (the molecular solver's normal/complementary operator
    // path, dmtk/molecular_system.hpp). init_hamiltonian ends init_iteration,
    // hamiltonian_product is the product the Lanczos solver and the one-state
    // diagonalization call, rotate_hamiltonian carries the Hamiltonian's
    // operators into a new block in rotate(). The defaults are the code that
    // stood at those places.
    virtual void init_hamiltonian(bool create_composite);
    virtual VectorState<T> hamiltonian_product(const VectorState<T>& vv);
    virtual void rotate_hamiltonian(int position, Block<T>& b, Basis& basis, Basis& rho_basis);

    System<T>& set_calc_gap(int b)
      {
        _calc_gap = b;
        if(b > 0){
          _ntargets = b+1;
          _target.resize(b+1);
          _target_weight.resize(b+1);
          _target_weight[0] = 0.875;
          for(int i = 0; i < b; i++) _target_weight[i+1] = 0.125/b;
        } else {
          _ntargets = 1;
          _target.resize(1);
          _target_weight.resize(1);
          _target_weight[0] = 1.0;
        }
        return *this;
      }
    int calc_gap() const { return _calc_gap; }
    bool project() const { return _project; }

    dmtk::Stopwatch timer() const { return _timer; }
    size_t sweep_direction() const { return dir; }
    int iteration() const { return iter; }

    void set_lanczos_tolerance(double tol) { _lanczos_tol = tol; }
    double lanczos_tolerance() const { return _lanczos_tol; }
    void set_lanczos_maxiter(int n) { _lanczos_maxiter = n; }
    int lanczos_maxiter() { return _lanczos_maxiter; }
    // Residual bound of the current solve (written for this tree). States
    // that only feed a density matrix are accepted once their Ritz energies
    // have settled to the tolerance and their residuals are below its square
    // root, the accuracy an energy settled to tol implies; the final sweep,
    // whose states are measured, certifies every residual down to the
    // tolerance itself. A negative bound disables the residual check.
    double lanczos_residual_tolerance() const
      { return _lanczos_certify ? _lanczos_tol : std::sqrt(_lanczos_tol); }
    double truncation_error() const { return _truncation_error; }
    double entropy() const { return _entropy; }

    // Streams

    void write_status() const
    {
      const char *file = tmp_files.get_filename("system.dat"); 
      ofstream outputfile(file,std::ios::out|std::ios::binary);
      if(!outputfile) {
        dmtk::err() << "*** ERROR: Could not open file " << file << '\n';
        return;
      }
      outputfile.write(reinterpret_cast<const char*>(&_numsweeps), sizeof(size_t));
      for(int i = 1; i <= _numsweeps; i++){
        int l = _sweeps(0,i);
        outputfile.write(reinterpret_cast<const char*>(&l), sizeof(int));
        l = _sweeps(1,i);
        outputfile.write(reinterpret_cast<const char*>(&l), sizeof(int));
      }
      outputfile.write(reinterpret_cast<const char*>(&_in_warmup), sizeof(bool));
      outputfile.write(reinterpret_cast<const char*>(&_sweep), sizeof(int));
      outputfile.write(reinterpret_cast<const char*>(&dir), sizeof(size_t));
      outputfile.write(reinterpret_cast<const char*>(&iter), sizeof(int));
      outputfile.write(reinterpret_cast<const char*>(&_ntargets), sizeof(int));
      qnt.write(outputfile);
      qn.write(outputfile);

      outputfile.close();
    }

    void read_status() 
    {
      const char *file = tmp_files.get_filename("system.dat"); 
      ifstream inputfile(file,std::ios::in|std::ios::binary);
      if(!inputfile) {
        dmtk::err() << "*** ERROR: Could not open file " << file << '\n';
        return;
      }
      inputfile.read(reinterpret_cast<char*>(&_numsweeps), sizeof(size_t));
      _sweeps.resize(2,_numsweeps+1);
      for(int i = 1; i <= _numsweeps; i++){
        int l;
        inputfile.read(reinterpret_cast<char*>(&l), sizeof(int));
        _sweeps(0,i) = l;
        inputfile.read(reinterpret_cast<char*>(&l), sizeof(int));
        _sweeps(1,i) = l;
      }
      inputfile.read(reinterpret_cast<char*>(&_in_warmup), sizeof(bool));
      inputfile.read(reinterpret_cast<char*>(&_sweep), sizeof(int));
      inputfile.read(reinterpret_cast<char*>(&dir), sizeof(size_t));
      inputfile.read(reinterpret_cast<char*>(&iter), sizeof(int));
      inputfile.read(reinterpret_cast<char*>(&_ntargets), sizeof(int));
      inputfile.read(reinterpret_cast<char*>(&_nstates), sizeof(int));
      qnt.read(inputfile);
      qn.read(inputfile);

      inputfile.close();
    }

    void write(std::ostream &s) const
    {
      s.write(reinterpret_cast<const char*>(&_numsweeps), sizeof(size_t));
      for(int i = 1; i <= _numsweeps; i++){
        int l = _sweeps(0,i);
        s.write(reinterpret_cast<const char*>(&l), sizeof(int));
        l = _sweeps(1,i);
        s.write(reinterpret_cast<const char*>(&l), sizeof(int));
      }
      s.write(reinterpret_cast<const char*>(&_in_warmup), sizeof(bool));
      s.write(reinterpret_cast<const char*>(&_sweep), sizeof(int));
      s.write(reinterpret_cast<const char*>(&dir), sizeof(size_t));
      s.write(reinterpret_cast<const char*>(&iter), sizeof(int));
      s.write(reinterpret_cast<const char*>(&_ntargets), sizeof(int));
      s.write(reinterpret_cast<const char*>(&_nstates), sizeof(int));
      {
        int nprop = static_cast<int>(_propagate_state.size());
        s.write(reinterpret_cast<const char*>(&nprop), sizeof(int));
      }
      qnt.write(s);
      qn.write(s);
      gs.write(s);
      for(int i = 0; i < _ntargets; i++){
        s.write(reinterpret_cast<const char*>(&_target_weight[i]), sizeof(double));
        _target[i].write(s);
      }
      for(int i = 0; i < _propagate_state.size(); i++){
        _propagate_state[i].write(s);
      }
      rho.write(s);
      rho_basis.write(s);
      rightblock.write(s);
      leftblock.write(s);
    }

    void read(std::istream &s)
    {
      s.read(reinterpret_cast<char*>(&_numsweeps), sizeof(size_t));
      _sweeps.resize(2,_numsweeps+1);
      for(int i = 1; i <= _numsweeps; i++){
        int l;
        s.read(reinterpret_cast<char*>(&l), sizeof(int));
        _sweeps(0,i) = l;
        s.read(reinterpret_cast<char*>(&l), sizeof(int));
        _sweeps(1,i) = l;
      }
      s.read(reinterpret_cast<char*>(&_in_warmup), sizeof(bool));
      s.read(reinterpret_cast<char*>(&_sweep), sizeof(int));
      s.read(reinterpret_cast<char*>(&dir), sizeof(size_t));
      s.read(reinterpret_cast<char*>(&iter), sizeof(int));
      s.read(reinterpret_cast<char*>(&_ntargets), sizeof(int));
      int nprop;
      s.read(reinterpret_cast<char*>(&nprop), sizeof(int));
      _propagate_state.resize(nprop);
      qnt.read(s);
      qn.read(s);
      gs.read(s);
      _target.resize(_ntargets);
      _target_weight.resize(_ntargets);
      for(int i = 0; i < _ntargets; i++){
        s.read(reinterpret_cast<char*>(&_target_weight[i]), sizeof(double));
        _target[i].read(s);
      }
      for(int i = 0; i < _propagate_state.size(); i++){
        _propagate_state[i].read(s);
      }
      rho.read(s);
      rho_basis.read(s);
      rightblock.read(s);
      leftblock.read(s);
    }

    void write_gs(const VectorState<T> &gs, int n, int position = RIGHT) const
    {
      char file[255];
      if(position == LEFT)
        snprintf(file, sizeof(file), "gs_%s_%i_l.dat",_name,n);
      else
        snprintf(file, sizeof(file), "gs_%s_%i_r.dat",_name,n);

      ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::binary);
      if(!outputfile) 
        dmtk::err() << "*** ERROR: Could not open file " << file << '\n';
      else{
#ifdef DMTK_DEBUG
        dmtk::out() << "WRITING GS " << file << '\n';
#endif // DMTK_DEBUG

        gs.write(outputfile);
  
        int _n = _target.size();
        outputfile.write(reinterpret_cast<const char*>(&_n), sizeof(int));
        for(int i = 0; i < _n; i++){
          _target[i].write(outputfile);
        }
  
        _n = _propagate_state.size();
        outputfile.write(reinterpret_cast<const char*>(&_n), sizeof(int));
        for(int i = 0; i < _propagate_state.size(); i++){
          _propagate_state[i].write(outputfile);
        }
        outputfile.close();

      }
    }

    void read_gs(VectorState<T> &gs, int n, int position = RIGHT)
    {
      char file[255];
      if(position == LEFT)
        snprintf(file, sizeof(file), "gs_%s_%i_l.dat",_name,n);
      else
        snprintf(file, sizeof(file), "gs_%s_%i_r.dat",_name,n);

      ifstream inputfile(tmp_files.get_filename(file),std::ios::in|std::ios::binary);
      if(!inputfile) 
        dmtk::err() << "*** ERROR: Could not open file " << file << '\n';
      else{
#ifdef DMTK_DEBUG
        dmtk::out() << "READING GS " << file << '\n';
#endif // DMTK_DEBUG

        gs.read(inputfile);
  
        inputfile.read(reinterpret_cast<char*>(&_ntargets), sizeof(int));
        _target.resize(_ntargets);
        _target_weight.resize(_ntargets);
        for(int i = 0; i < _ntargets; i++){
          _target[i].read(inputfile);
        }
  
        int nprop;
        inputfile.read(reinterpret_cast<char*>(&nprop), sizeof(int));
        _propagate_state.resize(nprop);
        for(int i = 0; i < _propagate_state.size(); i++){
          _propagate_state[i].read(inputfile);
        }

        inputfile.close();
      }
    }

    void write_rho(const BMatrix<T> &m, const Basis& basis, int n, int position) const
    {
      char file[255];
      if(position == LEFT)
        snprintf(file, sizeof(file), "rho_%s_%i_l.dat",_name,n);
      else
        snprintf(file, sizeof(file), "rho_%s_%i_r.dat",_name,n);

#ifdef DMTK_DEBUG
      dmtk::out() << "SAVING RHO " << file << '\n';
#endif // DMTK_DEBUG
      ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::binary);
      if(!outputfile) 
        dmtk::err() << "*** ERROR: Could not open file " << file << '\n';
      else{
        basis.write(outputfile);
        m.write(outputfile);
        outputfile.close();
      }
    }

    void read_rho(BMatrix<T> &m, Basis &basis, int n, int position)
    {
      char file[255];
      if(position == LEFT)
        snprintf(file, sizeof(file), "rho_%s_%i_l.dat",_name,n);
      else
        snprintf(file, sizeof(file), "rho_%s_%i_r.dat",_name,n);

#ifdef DMTK_DEBUG
      dmtk::out() << "READING RHO " << file << '\n';
#endif // DMTK_DEBUG
      ifstream inputfile(tmp_files.get_filename(file),std::ios::in|std::ios::binary);
      if(!inputfile) 
        dmtk::err() << "*** ERROR: Could not open file " << file << '\n';
      else{
        basis.read(inputfile);
        m.read(inputfile);
        inputfile.close();
      }
    }

    void write_block(const Block<T> &b, int n, int position) const
    {
      char file[255];
      if(position == LEFT)
        snprintf(file, sizeof(file), "block_%s_%i_l.dat",_name,n);
      else
        snprintf(file, sizeof(file), "block_%s_%i_r.dat",_name,n);

#ifdef DMTK_DEBUG
      dmtk::out() << "SAVING BLOCK " << file << '\n';
#endif // DMTK_DEBUG
      _discarded_blocks.erase(file); // Changed here (see final_sweep)
      ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::binary);
      if(!outputfile) 
        dmtk::err() << "*** ERROR: Could not open file " << file << '\n';
      else{
        b.write(outputfile);
        outputfile.close();
      }
    }

    void read_block(Block<T> &b, int n, int position)
    {
      if(n == 1){
        if(position == LEFT) 
          b = h.get_site(0);
        else
          b = h.get_site(lattice().size()-1);
        return; 
      }

      char file[255];
      if(position == LEFT)
        snprintf(file, sizeof(file), "block_%s_%i_l.dat",_name,n);
      else
        snprintf(file, sizeof(file), "block_%s_%i_r.dat",_name,n);

#ifdef DMTK_DEBUG
      dmtk::out() << "READING BLOCK " << file << '\n';
#endif // DMTK_DEBUG
      // Changed here (see final_sweep): a block final_sweep removed is gone,
      // not stale; refuse it rather than keep the block already in b.
      if(_discarded_blocks.count(file))
        throw std::logic_error(std::string("dmtk: ") + file + " was removed by final_sweep");
      ifstream inputfile(tmp_files.get_filename(file),std::ios::in|std::ios::binary);
      if(!inputfile) 
        dmtk::err() << "*** ERROR: Could not open file " << file << '\n';
      else{
        b.read(inputfile);
        inputfile.close();
      }
    }

    // Changed here (written for this tree, see final_sweep): removes the
    // file of block n on the given side; read_block refuses it until the
    // block is written again.
    void discard_block(int n, int position)
    {
      char file[255];
      if(position == LEFT)
        snprintf(file, sizeof(file), "block_%s_%i_l.dat",_name,n);
      else
        snprintf(file, sizeof(file), "block_%s_%i_r.dat",_name,n);
      std::error_code ec;
      std::filesystem::remove(tmp_files.get_filename(file), ec);
      _discarded_blocks.insert(file);
    }


    void write_iter(int position, bool symmetric = false) const
    {
      write_status();
      write_block(newblock,iter+1,position);
      write_rho(rho,rho_basis,iter+1,position);
      write_gs(gs,iter,position);
      if(symmetric){
        write_block(newblock,iter+1,1-position);
        write_rho(rho,rho_basis,iter+1,1-position);
        write_gs(gs,iter,1-position);
      }
    }

    System<T>& signal_connect(typename Signal<System<T> >::callback handler, size_t signal_id, void *data)
      {
         Signal<System<T> > signal(handler, this, signal_id, data);
         this->handlers.push_back(signal);
         return *this;
      }
    System<T>& signal_connect(typename Signal<System<T> >::callback handler, std::string signal_id, void *data)
      {
         int n = this->get_signal_id(signal_id);
         Signal<System<T> > signal(handler, this, n, data);
         this->handlers.push_back(signal);
         return *this;
      }

    bool signal_emit(size_t signal_id)
      {
         typename std::list<Signal<System<T> > >::reverse_iterator siter;
         for(siter = handlers.rbegin(); siter != handlers.rend(); siter++){
           Signal<System<T> > &signal = *siter;
           if(signal.signal_id() == signal_id) {
              if(!signal.emit()) return false;
           }
         }
         return true;
      }
    bool signal_emit(std::string signal_id)
      {
         int n = this->get_signal_id(signal_id);
         if(n >= 0) return this->signal_emit(n);
      }
    std::string signal_description(size_t signal_id)
      {
          return this->_signal_id[signal_id];
      }
    int signal_id(std::string signal_id)
      {
          return this->get_signal_id(signal_id);
      }
    System<T>& signal_disconnect_all()
      {
         this->handlers.clear();
         return *this;
      }
};

/////////////////////////////////////////////////////////////////////////
// warmup:
// warmup_loop: 
// First sweep through the lattice to build the blocks.
// Symmetric growing (1D): we first sweep to ls/2-1 and build right and
// left blocks at the same time. We then continue from ls/2-1 to ls-3 for
// the remaining left blocks. The system is then ready for the main loop.
// Non symmetric growing (2D): we first weep from left to right, using small
// 1-site/2-sites environment blocks on the right (to have an even total 
// number of sites all the time), until we reach the final size. We then sweep
// from right to left to build the right blocks. Finally, a sweep from left to 
// right to leave the system ready for the main loop. 
/////////////////////////////////////////////////////////////////////////
template <class T>
void
System<T>::warmup_loop(size_t t)
{
  Vector<QN> qns(lattice().ls());
  qns = qnt;
  _custom_qns = false;
  warmup_loop(t, qns);
}

template <class T>
void
System<T>::warmup_loop(size_t t, const Vector<QN> &qns)
{
  m = t;
  int ls = lattice().ls();
  _timer.restart();
  bool save_use_seed = _use_seed;
  _use_seed = false;

  char file[255];
  snprintf(file, sizeof(file), "iter_%s.dat",_name);
  ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::app);
  if(!outputfile) dmtk::err() << "*** ERROR: could not open " << file << endl;

  start();
//****************************************************************

// LEFT-TO-RIGHT sweep to get B(1)B(2)...B(L/2-1)
  dmtk::out() << "LEFT-TO-RIGHT sweep to get B(1)B(2)...B(L/2-1)\n";
  outputfile << "LEFT-TO-RIGHT sweep to get B(1)B(2)...B(L/2-1)\n";

  dir = LEFT2RIGHT;
  int sweep_max = full_sweep() ? ls-2 : ls-3;
//  int sweep_max = ls-2;
  int mid_size = (ls/2)*2 == ls ? ls/2 : ls/2+1;
  int max = _grow_symmetric ? mid_size : sweep_max;

  for(iter = 1; iter < max; iter++)
  {
    signal_emit(SYSTEM_SIGNAL_START_ITER);
//    qn = qns[_grow_symmetric?iter*2+2-1:(SGN(iter) == 1?iter+3:iter+2)];
    qn = qns[_grow_symmetric?iter*2+2-1:iter+3];

    if(_grow_symmetric)
      read_block(rightblock, 2*iter+2 <= ls ? iter : iter-1, RIGHT);
    else
      rightblock = h.get_site(iter+2);
/*
#ifdef GRANDCANONICAL
    read_block(rightblock, _grow_symmetric ? iter : 1, RIGHT);
#else // GRANDCANONICAL
    read_block(rightblock, _grow_symmetric ? iter : (SGN(iter) == 1 && _grand_canonical  == 0 ? 2 : 1), RIGHT);
#endif // GRANDCANONICAL
*/
    read_block(leftblock, iter, LEFT);
    const Block<T>& site1 = _grow_symmetric && _grow_outward ? h.get_site(ls/2) : h.get_site(iter);
    const Block<T>& site2 = _grow_symmetric ? (_grow_outward ? h.get_site(ls/2+1) : h.get_site(ls-iter-1)) : h.get_site(iter+1);

    _timer.mark();
    dmtk::out() << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;
    outputfile << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;

    m1 = leftblock.dim();
    m4 = rightblock.dim();
    m2 = site1.dim();
    m3 = site2.dim();

    init_iteration(leftblock, site1, site2, rightblock, false);
    diagonalize(false);

    m1 = m4 = std::min(m1*m2,m);

    truncate(LEFT, m1);
    if(_grow_outward){
      _b2 = &(h.get_site(ls/2-iter-2)); 
      _b1 = &leftblock; 
    }
    rotate(LEFT, newblock);
    write_iter(LEFT);


#ifndef GRANDCANONICAL 
    if((_grow_symmetric && iter < mid_size-1) || (_grand_canonical  == 0 && iter == 1)) { 
      truncate(RIGHT, m4);
      if(_grow_outward){
        _b3 = &(h.get_site(ls/2+iter+2)); 
        _b4 = &leftblock; 
      }
      rotate(RIGHT, newblock);
      write_iter(RIGHT);
    }
#endif //GRANDCANONICAL

    if(verbose() > 0) {
      dmtk::out() << "===========================================\n";
      dmtk::out() << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
      dmtk::out() << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
      outputfile << "===========================================\n";
    }
    signal_emit(SYSTEM_SIGNAL_END_ITER);
  }

  qn = qnt;
  if(_grow_symmetric){

//  LEFT-TO-RIGHT sweep to get B(L/2)B(L/2+1)...B(L-3)
    dmtk::out() << "LEFT-TO-RIGHT sweep to get B(L/2)B(L/2+1)...B(L-3)\n";
    outputfile << "LEFT-TO-RIGHT sweep to get B(L/2)B(L/2+1)...B(L-3)\n";

    for(iter = mid_size; iter < sweep_max; iter++)
    {
      signal_emit(SYSTEM_SIGNAL_START_ITER);
      read_block(rightblock, ls-iter-2, RIGHT);
      read_block(leftblock, iter, LEFT);
      const Block<T>& site1 = h.get_site(iter);
      const Block<T>& site2 = h.get_site(iter+1);
  
      _timer.mark();
      dmtk::out() << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;
      outputfile << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;
  
      m1 = leftblock.dim();
      m4 = rightblock.dim();
      m2 = site1.dim();
      m3 = site2.dim();
  
      init_iteration(leftblock, site1, site2, rightblock, false);
      diagonalize(false);
  
      m1 = std::min(m1*m2,m);
  
      truncate(LEFT, m1);
      rotate(LEFT, newblock);
      write_iter(LEFT);
  
      if(verbose() > 0) {
        dmtk::out() << "===========================================\n";
        dmtk::out() << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
        dmtk::out() << "===========================================\n";
        outputfile << "===========================================\n";
        outputfile << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
        outputfile << "===========================================\n";
      }
      signal_emit(SYSTEM_SIGNAL_END_ITER);
    }
  } else {
    dmtk::out() << "RIGHT-TO-LEFT sweep to get B(1)B(2)...B(L/2-1)\n";
    outputfile << "RIGHT-TO-LEFT sweep to get B(1)B(2)...B(L/2-1)\n";
    dir = RIGHT2LEFT;

    for(iter = 1; iter < sweep_max; iter++)
    {
      signal_emit(SYSTEM_SIGNAL_START_ITER);
      read_block(leftblock, ls-iter-2, LEFT);
      read_block(rightblock, iter, RIGHT);
      const Block<T>& site1 = h.get_site(ls-2-iter);
      const Block<T>& site2 = h.get_site(ls-2-iter+1);

      m1 = leftblock.dim();
      m4 = rightblock.dim();
      m2 = site1.dim();
      m3 = site2.dim();

      _timer.mark();
      dmtk::out() << "WARMUP RIGHT-TO-LEFT ITERATION " << iter << endl;
      outputfile << "WARMUP RIGHT-TO-LEFT ITERATION " << iter << endl;

      init_iteration(leftblock, site1, site2, rightblock, false); 
      diagonalize(false); 

      m4 = std::min(m4*m3,m);

      truncate(RIGHT, m4);
      rotate(RIGHT, newblock);
      write_iter(RIGHT);

      if(verbose() > 0) {
        dmtk::out() << "===========================================\n";
        dmtk::out() << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
        dmtk::out() << "===========================================\n";
        outputfile << "===========================================\n";
        outputfile << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
        outputfile << "===========================================\n";
      }
      signal_emit(SYSTEM_SIGNAL_END_ITER);
    }

    dir = LEFT2RIGHT;
    _use_seed = save_use_seed;

    for(iter = 1; iter < sweep_max; iter++)
    {
      signal_emit(SYSTEM_SIGNAL_START_ITER);
      read_block(rightblock, ls-iter-2, RIGHT);
      read_block(leftblock, iter, LEFT);
      const Block<T>& site1 = h.get_site(iter);
      const Block<T>& site2 = h.get_site(iter+1);
  
      m1 = leftblock.dim();
      m4 = rightblock.dim();
      m2 = site1.dim();
      m3 = site2.dim();
 
      _timer.mark();
      dmtk::out() << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;
      outputfile << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;
 
//      if(iter < sweep_max/2) {
//        init_iteration(leftblock, site1, site2, rightblock, false); 
//        diagonalize(false);
//      }else{
        init_iteration(leftblock, site1, site2, rightblock, save_use_seed); 
        diagonalize(_use_seed);
//      }
  
      m1 = std::min(m1*m2,m);
  
      truncate(LEFT, m1);
      rotate(LEFT, newblock);
      write_iter(LEFT);
  
      if(verbose() > 0) {
        dmtk::out() << "===========================================\n";
        dmtk::out() << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
        dmtk::out() << "===========================================\n";
        outputfile << "===========================================\n";
        outputfile << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
        outputfile << "===========================================\n";
      }
      signal_emit(SYSTEM_SIGNAL_END_ITER);
    }
  }
    
  qn = qnt;
  outputfile.close();
  _in_warmup = false;
  _use_seed = save_use_seed;
}

/////////////////////////////////////////////////////////////////////////
// sweep: 
// left-to-right and right-to-left sweeps, with t1(t2) number
// of states, respectively. We repeat this procedure several times
// in the main program with growing t1(t2) until we reach convergence
// with the final desired number of states.
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::sweep(size_t t1, size_t t2, size_t _dir, int start)
{
  int ls = lattice().ls();
  _in_warmup = false;
  _custom_qns = false;
  int sweep_max = full_sweep() ? ls-2 : ls-3;
  qn = qnt;

  char file[255];
  snprintf(file, sizeof(file), "iter_%s.dat",_name);
  ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::app);
  if(!outputfile) dmtk::err() << "*** ERROR: could not open " << file << endl;

  dmtk::out() << "===========================================\n";
  dmtk::out() << "NEW SWEEP\n";
  dmtk::out() << "===========================================\n";
  outputfile << "===========================================\n";
  outputfile << "NEW SWEEP\n";
  outputfile << "===========================================\n";

// RIGHT-TO-LEFT sweep
  
  dir = _dir;

  if(dir == RIGHT2LEFT){
    for(iter = start; iter < sweep_max; iter++)
    {
      signal_emit(SYSTEM_SIGNAL_START_ITER);
      read_block(leftblock, ls-iter-2, LEFT);
      read_block(rightblock, iter, RIGHT);
      const Block<T>& site1 = h.get_site(ls-2-iter);
      const Block<T>& site2 = h.get_site(ls-2-iter+1);

      m1 = leftblock.dim();
      m4 = rightblock.dim();
      m2 = site1.dim();
      m3 = site2.dim();

      _timer.mark();
      dmtk::out() << "RIGHT-TO-LEFT ITERATION " << iter << endl;
      outputfile << "RIGHT-TO-LEFT ITERATION " << iter << endl;

      init_iteration(leftblock, site1, site2, rightblock, _use_seed); 
      diagonalize(_use_seed); 

      m = t1;
      m4 = std::min(m4*m3,m);

      truncate(RIGHT, m4);
      rotate(RIGHT, newblock);
      write_iter(RIGHT);
      if(!_store_products){
        if (signal_emit(SYSTEM_SIGNAL_MEASURE)) measure();
      }

      if(verbose() > 0) {
        dmtk::out() << "===========================================\n";
        dmtk::out() << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
        dmtk::out() << "===========================================\n";
        outputfile << "===========================================\n";
        outputfile << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
        outputfile << "===========================================\n";
      }
      signal_emit(SYSTEM_SIGNAL_END_ITER);
    }

    start = 1;
  }

// LEFT-TO-RIGHT sweep

  dir = LEFT2RIGHT;

  for(iter = start; iter < sweep_max; iter++)
  {
    signal_emit(SYSTEM_SIGNAL_START_ITER);
    read_block(rightblock, ls-iter-2, RIGHT);
    read_block(leftblock, iter, LEFT);
    const Block<T>& site1 = h.get_site(iter);
    const Block<T>& site2 = h.get_site(iter+1);

    m1 = leftblock.dim();
    m4 = rightblock.dim();
    m2 = site1.dim();
    m3 = site2.dim();

    _timer.mark();
    dmtk::out() << "LEFT-TO-RIGHT ITERATION " << iter << endl;
    outputfile << "LEFT-TO-RIGHT ITERATION " << iter << endl;

    init_iteration(leftblock, site1, site2, rightblock, _use_seed); 
    diagonalize(_use_seed);

    m = t2;
    m1 = std::min(m1*m2,m);

    truncate(LEFT, m1);
    rotate(LEFT, newblock);
    write_iter(LEFT);
    if(!_store_products){
      if (signal_emit(SYSTEM_SIGNAL_MEASURE)) measure();
    }

    if(verbose() > 0) {
      dmtk::out() << "===========================================\n";
      dmtk::out() << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
      dmtk::out() << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
      outputfile << "===========================================\n";
    }
    signal_emit(SYSTEM_SIGNAL_END_ITER);
  }

  outputfile.close();
}

/////////////////////////////////////////////////////////////////////////
// init_correlations:
// We store in site blocks the operators used for measuring correlations
/////////////////////////////////////////////////////////////////////////
template <class T>
void
System<T>::init_correlations()
{

  corr.reorder_terms(true);
  return;
}

/////////////////////////////////////////////////////////////////////////
// final_sweep:
// Like a regular sweep but from 1 to ls/2-1 to build
// right and left symmetric blocks (we don't care about the others for the
// measurement). The main difference is that at each step we add the 
// measurement operators to the block (stored in the list "ops")  
/////////////////////////////////////////////////////////////////////////

template<class T>
void
System<T>::final_sweep(size_t t, size_t _dir, int _start, bool _rotate )
{
// Last iteration to get a symmetric block B(L/2-1)..B(L/2-1)

  int ls = lattice().ls();
  _in_warmup = false;
  // The states of this sweep are measured: certify their Lanczos residuals.
  const bool save_certify = _lanczos_certify;
  _lanczos_certify = true;

  char file[255];
  snprintf(file, sizeof(file), "iter_%s.dat",_name);
  ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::app);
  if(!outputfile) dmtk::err() << "*** ERROR: could not open " << file << endl;
  outputfile << "Last iteration to get a symmetric block B(L/2-1)..B(L/2-1)\n";

//****************************************************************
// We store in site blocks the operators used for measuring correlations
//****************************************************************
  init_correlations();
//****************************************************************

  dir = _dir;
  int sweep_max = full_sweep() ? ls-2 : ls-3;

  // Changed here: with stored products and a full sweep, nothing is measured
  // during the sweep, and the measure() that follows it reads the closing
  // configuration only: the left block of ls-3 sites, two sites and the last
  // site, a site block. Measurement operators rotated into the right blocks
  // of the right-to-left half, or into the left block grown at the last step,
  // never reach a measurement, so they are not built. Each left block of the
  // left-to-right half is read once, by the next step, so its file is removed
  // once read. A block's stored products of two-body terms grow as the fourth
  // power of its size (the molecular solver's RDMs), and these blocks were
  // most of the scratch. Measured values are bit for bit the same. The left
  // blocks are not on disk after this sweep; read_block refuses them.
  bool const trim = _store_products && full_sweep();

  if(dir == RIGHT2LEFT){
    for(iter = _start; iter < sweep_max; iter++) // ls/2-1; iter++);
    {
      signal_emit(SYSTEM_SIGNAL_START_ITER);
      read_block(rightblock, iter, RIGHT);
      read_block(leftblock, ls-iter-2, LEFT);
      const Block<T>& site1 = h.get_site(ls-2-iter);
      const Block<T>& site2 = h.get_site(ls-2-iter+1);

      m1 = leftblock.dim();
      m4 = rightblock.dim();
      m2 = site1.dim();
      m3 = site2.dim();
   
      _timer.mark();
      dmtk::out() << "FINAL SWEEP ITERATION " << endl;
      outputfile << "FINAL SWEEP ITERATION " << endl;
      dmtk::out() << "RIGHT-TO-LEFT ITERATION " << iter << endl;
      outputfile << "RIGHT-TO-LEFT ITERATION " << iter << endl;
  
      init_iteration(leftblock, site1, site2, rightblock, _use_seed);
  
      diagonalize(_use_seed);
  
      m = t;
      m4 = std::min(m4*m3,m);
  
      truncate(RIGHT, m4);
      _rotate_corr = !trim; // Changed here (see above)
      rotate(RIGHT, newblock);
      _rotate_corr = true;
      write_iter(RIGHT);
      if(!_store_products){
        if (signal_emit(SYSTEM_SIGNAL_MEASURE)) measure();
      }

      if(verbose() > 0) {
        dmtk::out() << "===========================================\n";
        dmtk::out() << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
        dmtk::out() << "===========================================\n";
        outputfile << "===========================================\n";
        outputfile << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
        outputfile << "===========================================\n";
      }
      signal_emit(SYSTEM_SIGNAL_END_ITER);
    }
    _start = 1;
  }

  dir = LEFT2RIGHT;
  sweep_max = full_sweep() ? ls-2 : ls/2-1;

  for(iter = _start; iter < sweep_max; iter++)
  {
    signal_emit(SYSTEM_SIGNAL_START_ITER);
    read_block(rightblock, ls-iter-2, RIGHT);
    read_block(leftblock, iter, LEFT);
    if(trim && iter > 1) discard_block(iter, LEFT); // Changed here (see above)
    const Block<T>& site1 = h.get_site(iter);
    const Block<T>& site2 = h.get_site(iter+1);

    m1 = leftblock.dim();
    m4 = rightblock.dim();
    m2 = site1.dim();
    m3 = site2.dim();

    _timer.mark();
    dmtk::out() << "FINAL SWEEP ITERATION " << endl;
    outputfile << "FINAL SWEEP ITERATION " << endl;
    dmtk::out() << "LEFT-TO-RIGHT ITERATION " << iter << endl;
    outputfile << "LEFT-TO-RIGHT ITERATION " << iter << endl;

    init_iteration(leftblock, site1, site2, rightblock, _use_seed); 

    diagonalize(_use_seed);

    m = t;
    m1 = std::min(m1*m2,m);

    truncate(LEFT, m1);
    _rotate_corr = !(trim && iter == sweep_max-1); // Changed here (see above)
    rotate(LEFT, newblock);
    _rotate_corr = true;
    write_iter(LEFT);
    if(!_store_products){
      if (signal_emit(SYSTEM_SIGNAL_MEASURE)) measure();
    }

    if(verbose() > 0) {
      dmtk::out() << "===========================================\n";
      dmtk::out() << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
      dmtk::out() << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << dmtk::format_elapsed(_timer.lap_elapsed()) << endl;
      outputfile << "===========================================\n";
    }
  }

// Calculate the ground state in the symmetric system

  if(!full_sweep()){
    dmtk::out() << "FINAL SWEEP - LAST ITERATION " << endl;
    outputfile << "FINAL SWEEP - LAST ITERATION " << endl;
    read_block(rightblock, IS_EVEN(ls) ? ls/2-1 : ls/2, RIGHT);
    read_block(leftblock, ls/2-1, LEFT);
    const Block<T>& site1 = h.get_site(ls/2-1);
    const Block<T>& site2 = h.get_site(ls/2);
    init_iteration(leftblock, site1, site2, rightblock, _use_seed);  
//  set_lanczos_tolerance(std::min(1.e-16,_lanczos_tol));
//    set_lanczos_maxiter(30);
    diagonalize(_use_seed);
    if(!_store_products){
      if (signal_emit(SYSTEM_SIGNAL_MEASURE)) measure();
    }
    if(_rotate){
      m = t;
      m1 = std::min(m1*m2,m);
                                                                               
      truncate(LEFT, m1);
      rotate(LEFT, newblock);
      write_iter(LEFT);
    } else {
      write_gs(gs, ls/2-1, LEFT);
      if(!_store_products){
        if (signal_emit(SYSTEM_SIGNAL_MEASURE)) measure();
      }
    }
  }
  signal_emit(SYSTEM_SIGNAL_END_ITER);
  _lanczos_certify = save_certify;
  
  outputfile.close();
}



/////////////////////////////////////////////////////////////////////////
// init_iteration: 
// initialize the blocks, quantum numbers, and transform the ground
// state of the previous iteration as a new seed for the new one.
/////////////////////////////////////////////////////////////////////////
template <class T>
void
System<T>::create_interactions(bool create_composite)
{
#ifndef USE_PRODUCT_DEFAULT
  if(this->verbose() > 0) dmtk::out() << "Creating interaction terms" << endl;
  hint(MASK_BLOCK1|MASK_BLOCK2);
  hint(MASK_BLOCK3|MASK_BLOCK4);
  hint(MASK_BLOCK2|MASK_BLOCK3, false);
  if(create_composite)
    build_composite_operators();
#endif
}


template <class T>
void
System<T>::init_iteration(const B&b1, const B&b2, const B&b3, const B&b4, bool use_seed, bool create_composite)
{
  _b1 = &b1;
  _b2 = &b2;
  _b3 = &b3;
  _b4 = &b4;
  int ltotal = 0, lnow = 0;
  for(int i = 0; i < h.lattice().ls(); i++) 
    ltotal += h.get_site(i).n_orbitals();
  for(int i = 0; i < size(); i++) 
    lnow += h.get_site(i).n_orbitals();

  const Block<T> *site_block = _b2;
  if(dir == RIGHT2LEFT) site_block = _b3;

  QN aux_qn(999);
  PackedBasis::const_iterator biter1, biter2, biter3, biter4;
  for(biter1 = b1.basis().subspace_begin(); biter1 != b1.basis().subspace_end(); biter1++){
    for(biter2 = b2.basis().subspace_begin(); biter2 != b2.basis().subspace_end(); biter2++){
      for(biter3 = b3.basis().subspace_begin(); biter3 != b3.basis().subspace_end(); biter3++){
        for(biter4 = b4.basis().subspace_begin(); biter4 != b4.basis().subspace_end(); biter4++){
          QN sqn = (*biter1).qn()+(*biter2).qn()+(*biter3).qn()+(*biter4).qn();
          double this_dist = 0;
          double dist = 0;
          for(int j = 0; j < QN::max_index(); j++){
            if((_grand_canonical & (1 << j)) == 0) continue;           
            double nd = (lnow * (static_cast<double>(qnt[j].get_twice())*0.5 / static_cast<double>(ltotal)));
            this_dist += (0.5*sqn[j].get_twice()-nd)*(0.5*sqn[j].get_twice()-nd);
            dist += (0.5*aux_qn[j].get_twice()-nd)*(0.5*aux_qn[j].get_twice()-nd);
          }
          if(this_dist < dist) aux_qn = sqn;
        }
      }
    }
  }
  
  if(_in_warmup && !_custom_qns) 
    qn = aux_qn;
  if(!_in_warmup)
    qn = qnt;

  if(verbose() > 0) 
    {
       dmtk::out() << "TOTAL SIZE : " << size() << " " << lnow << endl;
       dmtk::out() << "LEFT BLOCK SIZE : " << _b1->lattice().size() << " " << _b1->n_orbitals() << endl;
       dmtk::out() << "RIGHT BLOCK SIZE : " << _b4->lattice().size() << " " << _b4->n_orbitals() << endl;
       dmtk::out() << "|" << _b1->n_orbitals() << "|-|" << _b2->n_orbitals() << "|-|" << _b3->n_orbitals() << "|-|" << _b4->n_orbitals() << "|" << endl;
       dmtk::out() << "|" << _b1->lattice().size() << "|-|" << _b2->lattice().size() << "|-|" << _b3->lattice().size() << "|-|" << _b4->lattice().size() << "|" << endl;
       for(int i = 0; i < QN::max_index(); i++)
         if(_grand_canonical & (1 << i)) dmtk::out() << QN::qn_name(i) << " = " << qnt[i] << endl;
    }
  char file[255];
  snprintf(file, sizeof(file), "iter_%s.dat",_name);
  ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::app);
  if(!outputfile) dmtk::err() << "*** ERROR: could not open " << file << endl;
  outputfile << "TOTAL SIZE : " << size() << " " << lnow << endl;
  outputfile << "LEFT BLOCK SIZE : " << _b1->lattice().size() << " " << _b1->n_orbitals() << endl;
  outputfile << "RIGHT BLOCK SIZE : " << _b4->lattice().size() << " " << _b4->n_orbitals() << endl;
  outputfile << "|" << _b1->n_orbitals() << "|-|" << _b2->n_orbitals() << "|-|" << _b3->n_orbitals() << "|-|" << _b4->n_orbitals() << "|" << endl; 
  outputfile << "|" << _b1->lattice().size() << "|-|" << _b2->lattice().size() << "|-|" << _b3->lattice().size() << "|-|" << _b4->lattice().size() << "|" << endl; 

  if(verbose() > 0) 
    for(int i = 0; i < QN::max_index(); i++)
      if(_grand_canonical & (1 << i)) dmtk::out() << QN::qn_name(i) << " = " << qn[i] << endl;

  if(use_seed){
     BMatrix<T> rho1;
     Basis basis1;
     BMatrix<T> rho2;
     Basis basis2;

     dmtk::Stopwatch clock;
     clock.restart();

     if(_use_basic_seed){
       seed.set_qn_mask(qn, _grand_canonical);   
       seed.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
       seed = T(1);

     } else if(iter == 1){

       if(lattice().ls() != 4) {
         seed.set_qn_mask(qn, _grand_canonical);   
         seed.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);

         if(dir == RIGHT2LEFT){
           read_rho(rho1, basis1, lattice().ls()-3, LEFT);
           read_rho(rho2, basis2, 2, RIGHT);
           read_gs(gs, lattice().ls()-4,LEFT);
           gs.resize(_grand_canonical);

           if(verbose() > 0) dmtk::out() << "NEW SEED " << gs.size() << endl;
/*
           if(_target.size() > 1) {
             gs *= T(_target_weight[0]);
             for(int i = 1; i < _target.size(); i++){
               gs += T(_target_weight[i])*_target[i];
             }
           }
*/
           new_seed(gs, seed, rho1, rho2, basis1, basis2, LEFT);

           VectorState<T> aux;
           for(int i = 0; i < _propagate_state.size(); i++){
             aux = _propagate_state[i]; 
             aux.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
             new_seed(_propagate_state[i], aux, rho1, rho2, basis1, basis2, LEFT);
             _propagate_state[i] = aux;
           }
           if(verbose() > 0)
             dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
         } else {
           read_rho(rho1, basis1, lattice().ls()-3, RIGHT);
           read_rho(rho2, basis2, 2, LEFT);
           read_gs(gs, lattice().ls()-4,RIGHT);
           gs.resize(_grand_canonical);

           if(verbose() > 0) dmtk::out() << "NEW SEED " << gs.size() << endl;
/*
           if(_target.size() > 1){
             gs *= T(_target_weight[0]);
             for(int i = 1; i < _target.size(); i++){
               gs += T(_target_weight[i])*_target[i];
             }
           }
*/
           new_seed(gs, seed, rho1, rho2, basis1, basis2, RIGHT);

           VectorState<T> aux;
           for(int i = 0; i < _propagate_state.size(); i++){
             aux = _propagate_state[i]; 
             aux.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
             new_seed(_propagate_state[i], aux, rho1, rho2, basis1, basis2, RIGHT);
             _propagate_state[i] = aux;
           }
           if(verbose() > 0)
             dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
         }
       } else {
         seed = gs;
       }
     } else {

       seed.set_qn_mask(qn, _grand_canonical);   
       seed.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);

       if(dir == LEFT2RIGHT){
         read_rho(rho1, basis1, iter, LEFT);
         read_rho(rho2, basis2, lattice().ls()-iter-2+1, RIGHT);
         read_gs(gs, iter-1, LEFT);
         gs.resize(_grand_canonical);

         if(verbose() > 0) dmtk::out() << "NEW SEED " << gs.size() << endl;
/*
         if(_target.size() > 1){
           gs *= T(_target_weight[0]);
           for(int i = 1; i < _target.size(); i++){
             gs += T(_target_weight[i])*_target[i];
           }
         }
*/
         new_seed(gs, seed, rho1, rho2, basis1, basis2, LEFT);

         VectorState<T> aux;
         for(int i = 0; i < _propagate_state.size(); i++){
           aux = _propagate_state[i]; 
           aux.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
           new_seed(_propagate_state[i], aux, rho1, rho2, basis1, basis2, LEFT);
           _propagate_state[i] = aux;
         }
         if(verbose() > 0)
           dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
       } else {
         read_rho(rho1, basis1, iter, RIGHT);
         read_rho(rho2, basis2, lattice().ls()-iter-2+1, LEFT);
         read_gs(gs, iter-1, RIGHT);
         gs.resize(_grand_canonical);

         if(verbose() > 0) dmtk::out() << "NEW SEED " << gs.size() << endl;
/*
         if(_target.size() > 1){
           gs *= T(_target_weight[0]);
           for(int i = 1; i < _target.size(); i++){
             gs += T(_target_weight[i])*_target[i];
           }
         }
*/
         new_seed(gs, seed, rho1, rho2, basis1, basis2, RIGHT);

         VectorState<T> aux;
         for(int i = 0; i < _propagate_state.size(); i++){
           aux = _propagate_state[i]; 
           aux.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
           new_seed(_propagate_state[i], aux, rho1, rho2, basis1, basis2, RIGHT);
           _propagate_state[i] = aux;
         }
         if(verbose() > 0)
           dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
       }

     }
 
  }

  gs.set_qn_mask(qn, _grand_canonical);
  gs.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
  init_hamiltonian(create_composite);
}

template <class T>
void
System<T>::init_hamiltonian(bool create_composite)
{
  // Changed here: removed a duplicate of these two calls that rebuilt the same terms.
  create_interactions(create_composite);
  init_terms_composite(*this, this->gs);
}

#include "dmtk/lanczos.ipp"

/////////////////////////////////////////////////////////////////////////
// diagonalize:
// run lanczos to calculate the ground state
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::diagonalize(bool use_seed)
{   
  double tol = _lanczos_tol;
  Vector<double> a(100);
  Vector<double> b(100);
  Vector<double> e(_nstates);

  dmtk::Stopwatch clock;
  clock.restart();

  if(verbose() > 0)
    dmtk::out() << "Initialization time: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;

  clock.restart();

  _state.resize(_nstates);
  for(int i = 0; i < _nstates; i++){
     _state[i].set_qn_mask(qn, _grand_canonical);
     _state[i].resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
  }

//  if(verbose() > 0) 
    dmtk::out() << "NUMBER OF STATES: " << _b1->dim() << " " << _b2->dim() << " " << _b3->dim() << " " <<  _b4->dim() << " " << gs.size() << endl;

  int n = _lanczos_maxiter;
  if(_target_subspaces && _in_warmup && lattice().size() != size()) n = 0;

  char file[255];
  snprintf(file, sizeof(file), "vectors_%s.dat",_name);

  if(verbose() > 0) dmtk::out() << ">>>>>>> LANCZOS <<<<<<<" << endl;

  if(real(product(seed,seed)) < 0.1) use_seed = false;

  if(gs.size() > 1) {
//    verif_hamiltonian<T, System<T>, VectorState<T> >(*this, _state[0]);
    lanczos<T, System<T>, VectorState<T> >(*this, _state, seed, e, _nstates, a, b, n, tol, use_seed, true, tmp_files.get_filename(file), false, false, lanczos_residual_tolerance());
    gs = _state[0];
  } else {
    gs = T(1);
    VectorState<T> aux = gs;
    aux = product(*this,gs);
    T xe = product(aux,gs);
    e[0] = real(xe);
    _state[0] = gs; 
    dmtk::out() << "E0 = " << e[0] << endl;
    dmtk::out() << "ENER = " << e[0] << endl;
  }

  energy.resize(_calc_gap+1);
  energy[0] = e[0];

  if(verbose() && use_seed){
    dmtk::out() << "OVERLAP " << product(seed,gs) << endl;
  } 

  snprintf(file, sizeof(file), "iter_%s.dat",_name);
  ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::app);
  if(!outputfile) dmtk::err() << "*** ERROR: could not open " << file << endl;

  outputfile << _b1->lattice().size() << " " << _b2->lattice().size() << " " << _b3->lattice().size() << " " <<  _b4->lattice().size() << endl;
  outputfile << _b1->dim() << " " << _b2->dim() << " " << _b3->dim() << " " <<  _b4->dim() << " " << gs.size() << endl;

  outputfile.precision(10);
  outputfile << "ENER = " << e[0] << endl;
  if(verbose() > 0){
    outputfile << "Lanczos time: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
    dmtk::out() << "-------------------------------------------\n";
    dmtk::out() << "Lanczos time: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
    dmtk::out() << "-------------------------------------------\n";
  }
  outputfile.close();
  signal_emit(SYSTEM_SIGNAL_GS);
}

/////////////////////////////////////////////////////////////////////////
// build_target_states:
// build target states for the density matrix. 
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::build_target_states(int /*pos*/)
{
  if(verbose() > 0) {
    dmtk::out() << "-------------------------------------------\n";
    dmtk::out() << "Creating target states for the Density Matrix\n";
  }

  if(_calc_gap == 0){
    if(_ntargets > 1){
      _target.resize(std::min(_nstates,_ntargets));
      _target_weight.resize(std::min(_nstates,_ntargets));
      for(int i = 0; i < std::min(_nstates,_ntargets); i++){
        _target[i] = _state[i];
      }
    } else {
      _target[0] = gs;
      _target_weight[0] = double(1.);
    }
    return;
  }

  BMatrix<T> rho1;
  Basis basis1;
  BMatrix<T> rho2;
  Basis basis2;
  _target[0] = gs;

  energy.resize(_calc_gap+1);

  // Four sites keep four single-site blocks, so there is no block-basis change
  // for the seed transform. Ordinary four-site sweeps do not produce its
  // required two-site density matrices either. Use the existing random-seed
  // path instead of reading those absent matrices or starting at the penalized
  // ground state. Larger systems retain their transformed excited-state seeds.
  const bool use_seed = _use_seed && lattice().ls() != 4;
  for(int i = 0; i < _calc_gap; i++){
    int i1, i2;

    if(use_seed){
      seed = gs;
      if(iter == 1){
        if(dir == RIGHT2LEFT){
          i1 = lattice().ls() - 3;
          i2 = 2;
          read_rho(rho1, basis1, i1, LEFT);
          read_rho(rho2, basis2, i2, RIGHT);
          new_seed(_target[1+i], seed, rho1, rho2, basis1, basis2, LEFT);
        } else {
          i1 = lattice().ls() - 3;
          i2 = 2;
          read_rho(rho1, basis1, i1, RIGHT);
          read_rho(rho2, basis2, i2, LEFT);
          new_seed(_target[1+i], seed, rho1, rho2, basis1, basis2, RIGHT);
        }
      } else {
        if(dir == LEFT2RIGHT){
          i1 = iter;
          i2 = lattice().ls() - iter - 2 + 1;
          read_rho(rho1, basis1, i1, LEFT);
          read_rho(rho2, basis2, i2, RIGHT);
          new_seed(_target[1+i], seed, rho1, rho2, basis1, basis2, LEFT);
        } else {
          i1 = iter;
          i2 = lattice().ls() - iter - 2 + 1;
          read_rho(rho1, basis1, i1, RIGHT);
          read_rho(rho2, basis2, i2, LEFT);
          new_seed(_target[1+i], seed, rho1, rho2, basis1, basis2, RIGHT);
        }
      }
    }

    if(gs.size() <= i+1){
      _target[i+1] = gs;
      continue;
    }

    if(verbose() > 0) dmtk::out() << "<<<<<<< LANCZOS >>>>>>>\n";

    _project = true;
    double tol = _lanczos_tol;
    Vector<double> a(100);
    Vector<double> b(100);
    Vector<double> e(1);
    char file[255];
    snprintf(file, sizeof(file), "vectors_%s.dat",_name);
    int n = _lanczos_maxiter;

    _project_states.resize(i+1);
    for(int j = 0; j <= i; j++) _project_states[j] = &_target[j];

    lanczos<T, System<T>, VectorState<T> >(*this, _state, seed, e, 1, a, b, n, tol,
        use_seed, true, tmp_files.get_filename(file), false, true, lanczos_residual_tolerance());
    _target[i+1] = _state[0];

    _project = false;
    for(int j = 0; j <= i; j++)
      dmtk::out() << "EXCITED OVERLAP(" << j << ") = " << product(_target[j],_target[1+i]) << endl;
    dmtk::out() << "E(" << i+1 << ") = " << e[0] << endl;
    energy[i+1] = e[0];
    dmtk::out() << "GAP(" << i+1 << ") = " << energy[0] << " " << e[0] << " " << e[0] - energy[0] << endl;
  }

}
/////////////////////////////////////////////////////////////////////////
// truncate:
// build and diagonalize density matrix. 
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::truncate(int position, int new_size)
{
  int ls = lattice().ls();

  dmtk::Stopwatch clock;
  clock.restart();

//********************************************************************
//                   DENSITY MATRIX
//********************************************************************
  signal_emit(SYSTEM_SIGNAL_TRUNCATE);
//******* Calculate density matrix (block by block)

  Basis basis(_b1->basis(),_b2->basis());
  if(position == RIGHT) basis = Basis(_b3->basis(),_b4->basis());
  basis.reorder();
  rho.resize(basis);

//////////////////////////////////////////////////////////////////////////

//******* Create target states

  build_target_states(position);
  
//******* Measure control operators
  typename Hami<T>::iterator op_iter;
  for(op_iter = control_ops.begin(); op_iter != control_ops.end(); op_iter++){
    Term<T> &term = *op_iter;
    BasicOp<T> this_op = term[0].internals();
    bool info;
    if(dir == LEFT2RIGHT){
      if(iter > 1 && iter < ls-3){
        measure_operator(this_op.set_site(iter),NULL,info);
      } else if(iter == 1) {
        measure_operator(this_op.set_site(0),NULL,info);
        measure_operator(this_op.set_site(1),NULL,info);
      } else if(iter == ls-3){
        measure_operator(this_op.set_site(ls-3),NULL,info);
        measure_operator(this_op.set_site(ls-2),NULL,info);
        measure_operator(this_op.set_site(ls-1),NULL,info);
      }
    } else {
      if(iter > 1 && iter < ls-3){
        measure_operator(this_op.set_site(ls-2-iter),NULL,info);
      } else if(iter == 1) {
        measure_operator(this_op.set_site(ls-1),NULL,info);
        measure_operator(this_op.set_site(ls-2),NULL,info);
        measure_operator(this_op.set_site(ls-3),NULL,info);
      } else if(iter == ls-3){
        measure_operator(this_op.set_site(1),NULL,info);
        measure_operator(this_op.set_site(0),NULL,info);
      }
    }
  }

#ifndef USE_PRODUCT_DEFAULT
  if(_rotate_terms){
    H<T> haux(h23);
    hint(MASK_BLOCK2|MASK_BLOCK3);
    VectorState<T> v23 = gs.condense(MASK_BLOCK2|MASK_BLOCK3);
    VectorState<T> res23(v23);
    res23 = T(0);
    product(h23, v23, res23, BLOCK2, T(1), false);
    dmtk::out() << "LOCAL ENERGY " << (dir == LEFT2RIGHT ? iter : ls-2-iter) << " " << product(res23,v23)/product(v23,v23) << endl;
    h23 = haux;
  }
#endif // USE_PRODUCT_DEFAULT
//********************************************************************
//                   DENSITY MATRIX
//********************************************************************
//////////////////////////////////////////////////////////////////////////
  signal_emit(SYSTEM_SIGNAL_BUILD_DM);
  if(verbose() > 0)
    dmtk::out() << "Building Density Matrix in blocks \n";

  PackedBasis::iterator basis_iter;
  for(basis_iter = basis.subspace_begin(); basis_iter != basis.subspace_end(); basis_iter++){
    SubSpace s = *basis_iter;
    SubMatrix<T> block(s.qn(),s,s);
    rho.push_back(block);
  }

  for(int nt = 0; nt < _target.size(); nt++){
    if(_target_weight[nt] != 0.0) {
      VectorState<T> &s = _target[nt]; 
#ifdef USE_K
//    if(_in_warmup && lattice().size() != size()) s.randomize();
#endif // USE_K
      if(verbose() > 0)
        dmtk::out() << "Target vector " << nt << " " << _target_weight[nt] << endl;
      T norm = product(s,s);
      rho += s.density_matrix(position)*T(_target_weight[nt]/norm);
    }
  }

  if(verbose() > 0)
    dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
  clock.mark();

  T trace = T(0);
  for (SubMatrix<T> & block : rho) {
    trace += block.trace();
  }

  if(verbose() > 0)
    dmtk::out() << "TRACE = " << trace << endl;
  signal_emit(SYSTEM_SIGNAL_DM_READY);
//******* Diagonalize Density Matrix in blocks
  if(verbose() > 0)
    dmtk::out() << "Diagonalizing Density Matrix in blocks\n";
  Vector<double> w(basis.dim()); // eigenvalues
  rho = rho.diagonalize(w);

  if(verbose() > 0)
    dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;

//******* Calculate entropy

  _entropy = 0.0;
  for(int i = 0; i < w.size(); i++){
    if(w[i] >= 1.e-5) {
//      cout << "WARNING: w[i] <= 0 " << w[i] << endl;
      _entropy -= w[i]*log(w[i]);
    }
    if(_verbose > 0) dmtk::out() << "RHO EIGENVALUE " << i << " " << w[i] << endl;
  }
//  _entropy /= LOG2;
  dmtk::out() << "ITER = " << iter << " ENTROPY = " << _entropy << endl;
  

//******* Reorder eigenvalues and eigenvectors, and truncate
  if(verbose() > 0)
    dmtk::out() << "Reordering eigenvalues and eigenvectors, and truncating\n";

  Vector<size_t> idx(basis.dim());
  indexx<double, Vector<double> >(basis.dim(), w, idx, false);

  Vector<size_t> new_idx(basis.dim());
  for(int i = 0; i < basis.dim(); i++) new_idx(idx(i)) = i;

  if(_use_error){
    double error = 1.;
    for(int i = 0; i < basis.dim(); i++){
      error -= w(idx(i));
      if(error <= _error){
//        new_size = std::min(new_size, (int)((i+1)*1.1));
        new_size = static_cast<int>((i+1)*1.1);
        if(_error_max_size != -1)
          new_size = std::min(new_size, _error_max_size);
        dmtk::out() << "ERROR = " << _error << " " << error << endl;
        dmtk::out() << "NEW SIZE = " << i+1 << " " << new_size << endl;
        break;
      }
    }
  }

  Vector<double> vaux = w;
  rho_basis = basis;
  int new_dim = 0;
  typename BMatrix<T>::iterator biter = rho.begin();
//  bool found_first = false;
//  bool found_n = false;
//  bool use_2 = false;

//  if(_target_subspaces && _in_warmup && lattice().size() != size())
  if(_target_subspaces && _in_warmup)
//  if(_target_subspaces)
  {
    Vector<double> bweight(rho.size());
    Vector<size_t> bidx(rho.size());
    Vector<SubMatrix<T>*> bvector(rho.size());
    size_t ib = 0;
    while(biter != rho.end()){
      SubMatrix<T> &block(*biter);
      if(block.qn().n() > qnt.n()) { biter++; continue; }
//      if(block.qn().n() == qnt.n() && block.qn().kx() != qnt.kx()) { biter++; continue; }
      bvector[ib] = &block;
      Range col_range = block.col_range();
      bweight[ib] = 0;
      for(size_t col = col_range.begin(); col <= col_range.end(); col++)
        bweight[ib] += w(col);
      ib++;
      biter++;
    }
    int bsize = ib;
    indexx<double, Vector<double> >(bsize, bweight, bidx, false);
    int ngroup = std::max(_nsub,int(rho.size()/bsize));
    int nmax = 0;
    for(ib = 0; ib < bsize; ib++){
      SubMatrix<T> &block(*bvector(bidx(ib)));
      Range col_range = block.col_range();
      nmax += std::min(int(col_range.size()),ngroup);
      if(nmax > new_size) break;   
    }
    bsize = std::min(ib,rho.size());

    BMatrix<T> new_rho;
    Vector<double> waux(rho.size());
    for (SubMatrix<T> & block : rho) {
      bool found = false;
      for(ib = 0; ib < bsize; ib++){
        SubMatrix<T> *_block = bvector(bidx(ib));
        if(_block == &block) { found = true; break; }
      }
      if(!found) continue; 

      SubMatrix<T> aux_block(block); // copy of the original
      Range col_range = block.col_range();
      int new_col = 0;
      int first_col = new_dim;
      Vector<double> waux(col_range.size());
      waux(Range(0,col_range.size()-1)) = vaux(col_range);
      indexx<double, Vector<double> >(col_range.size(), waux, idx, false);

      int maxcol = std::min(size_t(ngroup), col_range.size());
      for(size_t col = 0; col < maxcol; col++){
        block.column(new_col) = aux_block.column(idx(col));
        w(new_dim) = waux(idx(col));
        rho_basis(new_dim) = basis(idx(col)+col_range.begin());
        new_col++;
        new_dim++;
      }
      if(first_col != new_dim){
        block.resize(Range(first_col,new_dim-1),block.row_range());
        new_rho.push_back(block);
      }
      if(new_dim >= new_size) break;
    }
    rho = new_rho;
  }
  else
  {
    while(biter != rho.end()){
      SubMatrix<T> &block(*biter);
      SubMatrix<T> aux_block(block); // copy of the original
      Range col_range = block.col_range();
      int new_col = 0;
      int first_col = new_dim;
      for(size_t col = col_range.begin(); col <= col_range.end(); col++){
        if(new_idx(col) < new_size){
           block.column(new_col) = aux_block.column(col-col_range.begin());
           w(new_dim) = vaux(col);
           rho_basis(new_dim) = basis(col);
           new_col++;
           new_dim++;
        }
      }
      if(first_col == new_dim){
        biter = rho.erase(biter);
      }else{
        block.resize(Range(first_col,new_dim-1),block.row_range());
        biter++;
      }
//      if(new_dim >= new_size) break;
    }
  }

  if(verbose() > 0) dmtk::out() << "NEW DIM " << new_dim << endl;

  rho_basis.resize(new_dim); // truncate basis
  rho.repack(rho_basis);

//******* Check precision after truncation

  double xweight = 0;
  for(int i = 0; i < new_dim; i++){
//    cout << i << " " << w(i) << " " << vaux(idx(i)) << endl;
    xweight += w(i);
  }

  if(verbose() > 0) dmtk::out() << "TOTAL WEIGHT = " << xweight << endl;

//  xweight = 0;
//  for(int i = new_size; i < w.size(); i++){
////    cout << i << " " << w(i) << " " << vaux(idx(i)) << endl;
//    xweight += w(i);
//  }
//
//  cout << "DISCARDED WEIGHT = " << xweight << endl;

  precision = 1.;
  for(int i = 0; i < new_dim; i++){
//    cout << i << " " << w(i) << " " << vaux(idx(i)) << endl;
    precision -= w(i);
  }
  _truncation_error = precision;

  dmtk::out() << "-------------------------------------------\n";
  dmtk::out() << "Truncation error = " << precision << endl;
  if(verbose() > 0)
    dmtk::out() << "Truncation time: " << dmtk::format_elapsed(clock.elapsed()) << endl;
  dmtk::out() << "-------------------------------------------\n";

  char file[255];
  snprintf(file, sizeof(file), "iter_%s.dat",_name);
  ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::app);
  if(!outputfile) dmtk::err() << "*** ERROR: could not open " << file << endl;

  outputfile.precision(10);

  outputfile << "Truncation error = " << precision << endl;
  outputfile.close();
}

/////////////////////////////////////////////////////////////////////////
// hint:
// build new interaction Hamiltonian between two blocks.
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::hint(int _mask, bool add_local_h, const Hami<T> *_h)
{
  if(!_rotate_terms) return;
#ifdef USE_PRODUCT_DEFAULT
  return;
#endif //USE_PRODUCT_DEFAULT
  if(!_h) _h = &this->h;
  const Hami<T> &this_hami = *_h;

  const Block<T> *pb1;
  const Block<T> *pb2;
  BasicOp<T> *_hij;
  size_t pos1, pos2;
  switch(_mask){
    case (MASK_BLOCK1|MASK_BLOCK2):
      pb1 = _b1;
      pb2 = _b2;
      _hij = &h12;
      pos1 = BLOCK1;
      pos2 = BLOCK2;
      break;
    case (MASK_BLOCK2|MASK_BLOCK3):
      pb1 = _b2;
      pb2 = _b3;
      _hij = &h23;
      pos1 = BLOCK2;
      pos2 = BLOCK3;
      break;
    case (MASK_BLOCK3|MASK_BLOCK4):
      pb1 = _b3;
      pb2 = _b4;
      _hij = &h34;
      pos1 = BLOCK3;
      pos2 = BLOCK4;
      break;
  }

  const Block<T> &b1(*pb1);
  const Block<T> &b2(*pb2);
  BasicOp<T> &hij(*_hij);
//  hij.dqn.kx2() = qnt.kx2();
//  hij.dqn.ky2() = qnt.ky2();

  Basis aux_basis(b1.basis(),b2.basis());
  aux_basis.reorder();
  hij.resize(aux_basis);
  hij = T(0);

// We use indentity matrix for the new operators since
// we don't want to rotate or truncate them

  BMatrix<T> aux_rho(aux_basis);
  aux_rho.clear();
  typename PackedBasis::const_iterator biter;
  for(biter = aux_basis.subspaces().begin(); biter != aux_basis.subspaces().end(); biter++){
    SubMatrix<T> block(biter->qn(),*biter,*biter);
    // identity block (the expression-template I<T>() helper was retired
    // with meta.h; fill the diagonal explicitly)
    block = T(0);
    for(size_t k = 0; k < std::min(block.rows(), block.cols()); ++k)
      block(k, k) = T(1);
    aux_rho.push_back(block);
  }

  dmtk::Stopwatch clock;
  clock.restart();
  if(add_local_h && apply_hami()){
    clock.mark();
    if(verbose() > 0)
      dmtk::out() << "OPERATOR H2 " << this_hami.name() << endl;
    const BasicOp<T> *hop = b2(this_hami);
    if(hop) new_operator(hij, *hop, aux_rho, aux_basis, RIGHT, T(1), false);
//    cout << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;

    clock.mark();
    if(verbose() > 0)
      dmtk::out() << "OPERATOR H1 " << this_hami.name() << endl;
    hop = b1(this_hami);
    if(hop) new_operator(hij, *hop, aux_rho, aux_basis, LEFT, T(1), false);
//   cout << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
  }

  // Changed here: operators are found through a per-call index (see
  // OpIndex) rather than a scan of the block lists per lookup, and a term's
  // factors are sorted into blocks before anything is copied, so the terms
  // that do not couple these two blocks, most of them, cost no copies.
  OpIndex index;
  typename Hami<T>::const_iterator hiter;
  for(hiter = this_hami.begin(); hiter != this_hami.end(); hiter++){
    const Term<T>& t = (*hiter);
//    if(T(t.coef()) == T(0)) continue;
    if(t.type() == TERM_EXTERN) continue;
    if(t.size() == 1 && t.type() != TERM_LOCAL && apply_hami() && t[0].name() != this_hami.name()) {
//      continue;
      const BasicOp<T>&top = t[0];
      size_t ib = block(top.site());
      if((ib == BLOCK2 && _mask == (MASK_BLOCK1|MASK_BLOCK2)) || (ib == BLOCK3 && _mask == (MASK_BLOCK3|MASK_BLOCK4)) || (ib == BLOCK1 && _mask == (MASK_BLOCK1|MASK_BLOCK2) && b1.lattice().size() == 1) || (ib == BLOCK4 && _mask == (MASK_BLOCK3|MASK_BLOCK4) && b2.lattice().size() == 1)){

        bool calc_hc = this_hami.use_hc();

        const BasicOp<T>* op= indexed_op(top, index);
        if(op){
          if (op->is_diagonal()) calc_hc = false;

          if(this->verbose() > 0){
            if(!calc_hc){
              dmtk::out() << "H12 TERM " << t.name(true) << endl;
            }else{
              dmtk::out() << "H12 TERM " << t.name(true) << " + h.c." << endl;
            }
          }

          clock.mark();

          size_t position = _mask;
          if(_mask == (MASK_BLOCK1|MASK_BLOCK2))
            position = ib == BLOCK2 ? RIGHT : LEFT;
          else
            position = ib == BLOCK4 ? RIGHT : LEFT;

          new_operator(hij, *op, aux_rho, aux_basis, position, T(t.coef()));
//          cout << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
        }
      }
    }
    if(t.size() >= 2 && ((t.type() == TERM_PRODUCT && apply_hami()) || (t.type() == TERM_EXTERN && apply_extern()))){
      bool found1 = false;
      bool found2 = false;
      bool found = true;
      Term<T> aux_term1, aux_term2;
      int bmask = 0;
      for(int i = 0; i < t.size(); i++){
        size_t const ib = block(t[i].site());
        bmask |= mask(ib);
//        cout << top.description() << " " << pos1 << " " << pos2 << " " << block(top.site()) << endl;
        if(ib == pos1) found1 = true;
        if(ib == pos2) found2 = true;
        if(ib == BLOCK_NONE){
          found = false;
          break;
        }
      }
      if(found && (bmask & _mask) == bmask && found1 && found2){ // we have a new composite operator
        for(int i = 0; i < t.size(); i++){
          size_t const ib = block(t[i].site());
          if(ib == pos1) aux_term1 *= t[i].internals();
          if(ib == pos2) aux_term2 *= t[i].internals();
        }

        BasicOp<T> top1, top2;
        if(aux_term1.size() == 1) // the piece contains a single operator
          top1 = aux_term1[0].internals();
        else
          top1 = aux_term1; // the piece contains more than one operator
        if(aux_term2.size() == 1) // the piece contains a single operator
          top2 = aux_term2[0].internals();
        else
          top2 = aux_term2; // the piece contains more than one operator
        const BasicOp<T>* op1 = indexed_op(top1, index);
        const BasicOp<T>* op2 = indexed_op(top2, index);

        if(!op1) {
          dmtk::out() << "ERROR : Operator " << top1.description() << " in term " << t.description() << " not found\n";
          continue;
        }
        if(!op2) {
          dmtk::out() << "ERROR : Operator " << top2.description() << " in term " << t.description() << " not found\n";
          continue;
        }

        bool calc_hc = this_hami.use_hc();
//        if (op1->is_diagonal() && op2->is_diagonal()) calc_hc = false;
        if (t.is_diagonal()) calc_hc = false;
//        if (t.size() > 2) {
////          calc_hc = false; // you need to include h.c. terms explicitly
//          bool is_diag = true;
//          for(int i = 0; i < t.size(); i++)
//            if(t[i].is_diagonal()) { is_diag = false; break; }
//          if(is_diag) calc_hc = false;
//        }

        if(op1 && op2){
          clock.mark();

          if(this->verbose() > 0){
            if(!calc_hc){
              dmtk::out() << "H12 TERM " << t.description() << endl;
            } else {
              dmtk::out() << "H12 TERM " << t.description() << " + h.c." << endl;
            }
          }
          new_operator(hij, *op1, *op2, pos1, pos2, aux_rho, aux_basis, T(t.coef()), calc_hc);

//          cout << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
        }
      }
    }

  }
}

//////////////////////////////////////////////////////////////////////////
// product_term:
//////////////////////////////////////////////////////////////////////////
template<class T>
void
product_term(const System<T> &ss, const AuxTerm<T> &auxt, const Term<T> &t, const VectorState<T> &vv, VectorState<T> &res)
{
  dmtk::Stopwatch aux_timer;
#ifdef WITH_PTHREADS
  DMTKglobals<T> *globals = NULL;
#else
  DMTKglobals<T> *globals = get_globals(T(0));
#endif

  aux_timer.restart();

  bool calc_hc = ss.h.use_hc();
// TODO WARNING: This may cause problems ?
//#warning If you are using use_hc and you have problems, make sure you are not including the h.c. operator in the Hamiltonian
  int ndiag = 0;
  if(t.is_diagonal()) calc_hc = false;
/*
  for(int i = 0; i < t.size(); i++){
    if(t[i].is_diagonal()) ndiag++;
  }
  if(ndiag == t.size()) calc_hc = false;
*/
///////////////////////////////////////////////
  if(ss.verbose() > 0){
    if(!calc_hc){
      dmtk::out() << "APPLYING " << t.description() << endl;
    }else{
      dmtk::out() << "APPLYING " << t.description() << " + h.c." << endl;
    }
  }

  if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2) || auxt.mask == (MASK_BLOCK2|MASK_BLOCK3) || auxt.mask == (MASK_BLOCK3|MASK_BLOCK4) || auxt.mask == (MASK_BLOCK1|MASK_BLOCK4) || auxt.mask == (MASK_BLOCK1|MASK_BLOCK3) || auxt.mask == (MASK_BLOCK2|MASK_BLOCK4)){

    aux_timer.mark();

    if(ss.use_composite()){
//      if(auxt.ref_op[0]->is_diagonal() && auxt.sum_op.is_diagonal()) calc_hc = false;
      if(auxt.b1 < auxt.b2)
        product(*(auxt.ref_op[0]), auxt.sum_op, vv, res, auxt.b1, auxt.b2, T(1), calc_hc, globals, false);
      else
        product(auxt.sum_op, *(auxt.ref_op[0]), vv, res, auxt.b2, auxt.b1, T(1), calc_hc, globals, false);
    } else {
//      if(auxt.op[size_t(auxt.b1)-1]->is_diagonal() && auxt.op[size_t(auxt.b2)-1]->is_diagonal()) calc_hc = false;
      if(auxt.b1 < auxt.b2)
        product(*auxt.op[size_t(auxt.b1)-1], *auxt.op[size_t(auxt.b2)-1], vv, res, auxt.b1, auxt.b2, t.coef(), calc_hc, globals, false);
      else
        product(*auxt.op[size_t(auxt.b2)-1], *auxt.op[size_t(auxt.b1)-1], vv, res, auxt.b2, auxt.b1, t.coef(), calc_hc, globals, false);
    }

    if(ss.verbose() > 0)
      dmtk::out() << "PRODUCT TIME= " << dmtk::format_elapsed(aux_timer.lap_elapsed()) << endl;
  }
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
    aux_timer.mark();

//    if(auxt.op[0]->is_diagonal() && auxt.op[1]->is_diagonal() && auxt.op[2]->is_diagonal() && auxt.op[3]->is_diagonal()) calc_hc = false;
    product(*(auxt.op[0]), *(auxt.op[1]), *(auxt.op[2]), *(auxt.op[3]), vv, res, BLOCK1, BLOCK2, BLOCK3, BLOCK4, T(t.coef()), calc_hc, globals);
    if(ss.verbose() > 0)
      dmtk::out() << "PRODUCT TIME= " << dmtk::format_elapsed(aux_timer.lap_elapsed()) << endl;
  }
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3)){
    aux_timer.mark();
    if(ss.use_composite()){
//      if(auxt.sum_op.is_diagonal() && auxt.ref_op[0]->is_diagonal() && auxt.ref_op[1]->is_diagonal()) calc_hc = false;
      product(auxt.sum_op, *(auxt.ref_op[0]), *(auxt.ref_op[1]), vv, res, BLOCK1, BLOCK2, BLOCK3, T(1), calc_hc, globals);
    } else {
//      if(auxt.op[0]->is_diagonal() && auxt.op[1]->is_diagonal() && auxt.op[2]->is_diagonal()) calc_hc = false;
      product(*(auxt.op[0]), *(auxt.op[1]), *(auxt.op[2]), vv, res, BLOCK1, BLOCK2, BLOCK3, T(t.coef()), calc_hc, globals);
    }

    if(ss.verbose() > 0)
      dmtk::out() << "PRODUCT TIME= " << dmtk::format_elapsed(aux_timer.lap_elapsed()) << endl;
  }
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4)){
    aux_timer.mark();

    if(ss.use_composite()){
//      if(auxt.sum_op.is_diagonal() && auxt.ref_op[0]->is_diagonal() && auxt.ref_op[1]->is_diagonal()) calc_hc = false;
      if(auxt.b2 < auxt.b3){
        product(auxt.sum_op, *(auxt.ref_op[0]), *(auxt.ref_op[1]), vv, res, BLOCK1, BLOCK2, BLOCK4, T(1), calc_hc, globals);
      } else {
        product(*(auxt.ref_op[1]), *(auxt.ref_op[0]), auxt.sum_op, vv, res, BLOCK1, BLOCK2, BLOCK4, T(1), calc_hc, globals);
      }
    } else {
//      if(auxt.op[0]->is_diagonal() && auxt.op[1]->is_diagonal() && auxt.op[3]->is_diagonal()) calc_hc = false;
      product(*(auxt.op[0]), *(auxt.op[1]), *(auxt.op[3]), vv, res, BLOCK1, BLOCK2, BLOCK4, T(t.coef()), calc_hc, globals);
    }

    if(ss.verbose() > 0)
      dmtk::out() << "PRODUCT TIME= " << dmtk::format_elapsed(aux_timer.lap_elapsed()) << endl;
  }
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4)){
    aux_timer.mark();

    if(ss.use_composite()){
//      if(auxt.sum_op.is_diagonal() && auxt.ref_op[0]->is_diagonal() && auxt.ref_op[1]->is_diagonal()) calc_hc = false;
      if(auxt.b2 < auxt.b3){
        product(auxt.sum_op, *(auxt.ref_op[0]), *(auxt.ref_op[1]), vv, res, BLOCK1, BLOCK3, BLOCK4, T(1), calc_hc, globals);
      } else {
        product(*(auxt.ref_op[1]), *(auxt.ref_op[0]), auxt.sum_op, vv, res, BLOCK1, BLOCK3, BLOCK4, T(1), calc_hc, globals);
      }
    } else {
//      if(auxt.op[0]->is_diagonal() && auxt.op[2]->is_diagonal() && auxt.op[3]->is_diagonal()) calc_hc = false;
      product(*(auxt.op[0]), *(auxt.op[2]), *(auxt.op[3]), vv, res, BLOCK1, BLOCK3, BLOCK4, T(t.coef()), calc_hc, globals);
    }

    if(ss.verbose() > 0)
      dmtk::out() << "PRODUCT TIME= " << dmtk::format_elapsed(aux_timer.lap_elapsed()) << endl;
  }
  else if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
    aux_timer.mark();

    if(ss.use_composite()) {
//      if(auxt.sum_op.is_diagonal() && auxt.ref_op[0]->is_diagonal() && auxt.ref_op[1]->is_diagonal()) calc_hc = false;
      product(*(auxt.ref_op[0]), *(auxt.ref_op[1]), auxt.sum_op, vv, res, BLOCK2, BLOCK3, BLOCK4, T(1), calc_hc, globals);
    } else {
//      if(auxt.op[1]->is_diagonal() && auxt.op[2]->is_diagonal() && auxt.op[3]->is_diagonal()) calc_hc = false;
      product(*(auxt.op[1]), *(auxt.op[2]), *(auxt.op[3]), vv, res, BLOCK2, BLOCK3, BLOCK4, T(t.coef()), calc_hc, globals);
    }

    if(ss.verbose() > 0)
      dmtk::out() << "PRODUCT TIME= " << dmtk::format_elapsed(aux_timer.lap_elapsed()) << endl;
  }
}
/////////////////////////////////////////////////////////////////////////
// rotate:
// build new block and operators, truncating with the density matrix.
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::rotate(int position, Block<T>& b)
{
  signal_emit(SYSTEM_SIGNAL_ROTATE);

  const Block<T> *pb1 = _b1;
  const Block<T> *pb2 = _b2;

  if(position == RIGHT){
    pb1 = _b3;
    pb2 = _b4;
  }

  dmtk::Stopwatch clock;
  clock.restart();

  const Block<T> &b1(*pb1);
  const Block<T> &b2(*pb2);

// We build new lattice for the new block 

  b = Block<T>(rho_basis, b1.n_orbitals()+b2.n_orbitals());
  b.set_single_site(false);
  Lattice l, l2;

  if(position == LEFT){
    l = _b1->lattice();
    l2 = _b2->lattice();
  } else {
    l = _b4->lattice();
    l2 = _b3->lattice();
  }
  for(int is = 0; is < l2.size(); is++){
    l.push_back(Site());
    l.nx() += 1;
  }
  b.set_lattice(l);

  b = T(0);
  Basis basis(b1.basis(),b2.basis());
  basis.reorder();

/////////////////////////////////////////////////////////////
  if(verbose() == 3)
  {
    Vector<const Block<T>*> _b(4);
    _b[0] = _b1; 
    _b[1] = _b2; 
    _b[2] = _b3; 
    _b[3] = _b4; 
    for(int i = 0; i < 4; i++){
      dmtk::out() << "BLOCK " << i+1 << endl;
      typename Block<T>::const_iterator biter;
      for(biter = _b[i]->begin(); biter != _b[i]->end(); biter++){
        dmtk::out() << "OPERATOR " << (*biter).name() << "(" << (*biter).site() << "," << (*biter).internal_site() << ")" << endl;
      }
    }
    dmtk::out() << "-------------\n";
  }
/////////////////////////////////////////////////////////////
//********************************************************************
// Test using identity
//********************************************************************
// TODO: comment this lines 
//
//  typename BMatrix<T>::iterator biter = rho.begin();
//
//  for(biter = rho.begin(); biter != rho.end(); biter++){
//    SubMatrix<T> &block(*biter);
//    block = I<T>(); 
//  }
//
  rotate_hamiltonian(position, b, basis, rho_basis);
// ***************************************************************************
// Correlations
// ***************************************************************************
  // Changed here: skipped where final_sweep knows no measurement reads them.
  if(_rotate_corr) rotate_corr(position, b, basis, rho_basis, &corr);
// ***************************************************************************
  rho_basis = basis;

  if(verbose() > 0) {
    dmtk::out() << "-------------------------------------------\n";
    dmtk::out() << "Rotation time: " << dmtk::format_elapsed(clock.elapsed()) << endl;
    dmtk::out() << "-------------------------------------------\n";
  }
  char file[255];
  snprintf(file, sizeof(file), "iter_%s.dat",_name);
  ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::app);
  if(!outputfile) dmtk::err() << "*** ERROR: could not open " << file << endl;

  outputfile << "Rotation time: " << dmtk::format_elapsed(clock.elapsed()) << endl;
  outputfile.close();
}

template<class T>
void
System<T>::rotate_hamiltonian(int position, Block<T>& b, Basis& basis, Basis& rho_basis)
{
// ***************************************************************************
// Hamiltonian terms for interactions between blocks
// ***************************************************************************
  rotate_terms(position, b, basis, rho_basis, nullptr);
// ***************************************************************************
// Local Hamiltonian
// ***************************************************************************
  if(verbose() > 0)
    dmtk::out() << "NEW BLOCK: HAMILTONIAN TERMS: " << this->h.name() << endl;
  rotate_hami(position, b, basis, rho_basis, nullptr); 
// ***************************************************************************
// Other Hamiltonians 
// ***************************************************************************
  for (Hami<T> & h : hs) {
    rotate_terms(position, b, basis, rho_basis, &h);
    rotate_hami(position, b, basis, rho_basis, &h);
  }
}

/////////////////////////////////////////////////////////////////////////
// rotate_terms:
// rotate hamiltonian terms 
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::rotate_terms(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *_this_hami)
{
  if(!_this_hami) _this_hami = &this->h;
  const Hami<T> &this_hami = *_this_hami;

  int _mask = MASK_BLOCK1 | MASK_BLOCK2;
  size_t pos1 = BLOCK1;
  size_t pos2 = BLOCK2;
  size_t the_block = BLOCK1;
  size_t site_block = BLOCK2;

  if(position == RIGHT){
    _mask = MASK_BLOCK3 | MASK_BLOCK4;
    pos1 = BLOCK3;
    pos2 = BLOCK4;
    site_block = BLOCK3;
    the_block = BLOCK4;
  }

  dmtk::Stopwatch clock;
  clock.restart();

  // Changed here: operators are found through a per-call index over the
  // four blocks (see OpIndex), whether the new block b already holds an
  // operator is answered from a key set kept next to it rather than a scan
  // of its list, and a term's factors are sorted into blocks before any of
  // them is copied. The nested calls below add to b, so they drop the set.
  OpIndex index;
  std::set<OpKey> in_b;
  bool in_b_valid = false;
  auto b_contains = [&](const BasicOp<T> &op) {
    if(!in_b_valid){
      in_b.clear();
      for(const BasicOp<T> &o : b) in_b.emplace(o.name(), o.site(), o.internal_site());
      in_b_valid = true;
    }
    return in_b.count(OpKey(op.name(), op.site(), op.internal_site())) != 0;
  };
  auto b_push = [&](const BasicOp<T> &op) {
    b.push_back(op);
    if(in_b_valid) in_b.emplace(op.name(), op.site(), op.internal_site());
  };

//********************************************************************
//                   NEW BLOCK
//********************************************************************
//******* Change of basis

  if(verbose() > 0) dmtk::out() << "NEW BLOCK: SINGLE OPERATORS\n";

// ***************************************************************************
// We add to the outer block the operators necessary to build the interaction
// between blocks. We take them from the single-site blocks and we add them 
// to the outer block. 
// We also add local terms of the hamiltonian on site block
// ***************************************************************************

  for (Term<T> const & t : this_hami) {
//    if(T(t.coef()) == T(0)) continue;
    if(t.size() == 1 && t[0].name() == "I") continue;
  
    bool found = false;
    const BasicOp<T> *_op = nullptr;

    int bmask = 0;

    for (BasicOp<T> const & op : t) {
      if (op.is_hami()) {
        rotate_hami(position, b, basis, rho_basis, &op.hami());
        rotate_terms(position, b, basis, rho_basis, &op.hami());
        in_b_valid = false;
      } else {
        bmask |= mask(block(op.site()));
      }
    }

    if(position == LEFT && (bmask & MASK_BLOCK2) && !(bmask & MASK_BLOCK1) && bmask != MASK_BLOCK2)
       found = true;
    if(position == RIGHT && (bmask & MASK_BLOCK3) && !(bmask & MASK_BLOCK4) && bmask != MASK_BLOCK3)
       found = true;

    if(found){
      BasicOp<T> real_op;
      for (BasicOp<T> const & op : t) {
        if (block(op.site()) == site_block && !op.is_hami()) {
          _op = indexed_op(op, index);
          if (!_op) dmtk::out() << "WARNING : Operator " << op.description() << " not found\n";
          real_op = op;
          real_op.dqn = _op->dqn;
          break;
        }
      }
      bool add = b_contains(real_op);

      if(!add && _op){ // If we haven't added it yet, we do it now
        real_op.resize(rho_basis);

        if(verbose() > 0)
          dmtk::out() << "ADDING NEW OPERATOR " << real_op.description() << endl;
        clock.mark();

        new_operator(real_op, *_op, rho, basis, 1-position);
        b_push(real_op);
        if(verbose() > 0)
          dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
      }
    }
  }

// ***************************************************************************
// Composite operators for interactions 
// ***************************************************************************
// 1- We look for terms that involve pieces in one block and pieces
//    in another block (not a site block) and we rotate them 
//    We also add TermType::TERM_LOCAL and TermType::TERM_EXTERN that are in the block

  for (Term<T> const & t : this_hami) {
    bool found = false;
    if((t.type() == TermType::TERM_PRODUCT && t.size() >= 2) || t.type() == TermType::TERM_LOCAL || t.type() == TermType::TERM_EXTERN){
      Term<T> aux_term;
      int bmask = 0;
      for(int i = 0; i < t.size(); i++){
        size_t const ib = block(t[i].site());
        if(ib == the_block) found = true;
        bmask |= mask(ib);
      }

      bool doit = false;
      if(found && position == LEFT && (bmask & MASK_BLOCK1) && ((!(bmask & MASK_BLOCK2) && bmask != MASK_BLOCK1) || (bmask == MASK_BLOCK1 && (t.type() == TermType::TERM_LOCAL || t.type() == TermType::TERM_EXTERN))))
         doit = true;
      if(found && position == RIGHT && (bmask & MASK_BLOCK4) && ((!(bmask & MASK_BLOCK3) && bmask != MASK_BLOCK4) || (bmask == MASK_BLOCK4 && (t.type() == TermType::TERM_LOCAL || t.type() == TermType::TERM_EXTERN))))
         doit = true;

      if(doit){ // we found a piece of a composite operator
        for(int i = 0; i < t.size(); i++)
          if(block(t[i].site()) == the_block) aux_term *= t[i];
        BasicOp<T> top1;
        if(aux_term.size() == 1) // the piece contains a single operator
          top1 = aux_term[0];
        else
          top1 = aux_term; // the piece contains more than one operator

        const BasicOp<T>* op1 = indexed_op(top1, index);

        if(!op1) {
          dmtk::out() << "ERROR 1: Operator " << top1.description() << " not found\n";
          continue;
        }

        BasicOp<T> new_op(top1); 
        new_op.dqn = op1->dqn;
        bool add = b_contains(new_op);

        if(!add && op1){
          new_op.resize(rho_basis);
          clock.mark();
  
          if(verbose() > 0)
            dmtk::out() << "NEW OPERATOR " << new_op.description() << endl;
          new_op = T(0);
          new_operator(new_op, *op1, rho, basis, position);
          b_push(new_op);
          if(verbose() > 0)
            dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
        }
      }
    }
  }

// 2- We look for terms that involve pieces in one block and pieces
//    in a site block and we put the pieces together

  for (Term<T> const & t : this_hami) {
    bool found1 = false;
    bool found2 = false;
    if(t.size() >= 2){
      Term<T> aux_term1, aux_term2;
      int bmask = 0;
      for(int i = 0; i < t.size(); i++){
        size_t const ib = block(t[i].site());
        bmask |= mask(ib);
        
        if(ib == pos1) found1 = true;
        if(ib == pos2) found2 = true;
      }
//    if bmask == _mask, this term will go into the Hamiltonian
//    unless it is TermType::TERM_EXTERN
      if(((bmask != _mask) || (bmask == _mask && t.type() == TermType::TERM_EXTERN)) && found1 && found2){ // we have a new composite operator
        for(int i = 0; i < t.size(); i++){
          size_t const ib = block(t[i].site());
          if(ib == pos1) aux_term1 *= t[i];
          if(ib == pos2) aux_term2 *= t[i];
        }
        BasicOp<T> top1(aux_term1);
        if(aux_term1.size() == 1)
          top1 = aux_term1[0];
        BasicOp<T> top2(aux_term2);
        if(aux_term2.size() == 1)
          top2 = aux_term2[0];
        const BasicOp<T>* op1 = indexed_op(top1, index);
        const BasicOp<T>* op2 = indexed_op(top2, index);

        if(!op1) {
          dmtk::out() << "ERROR 2: Operator " << top1.description() << " in term " << t.description() << " not found\n";
          continue;
        }
        if(!op2) {
          dmtk::out() << "ERROR 2: Operator " << top2.description() << " in term " << t.description() << " not found\n";
          continue;
        }

        Term<T> new_term; 
        for(int i = 0; i < t.size(); i++){
          size_t const ib = block(t[i].site());
          if(ib == pos1 || ib == pos2){
            new_term *= t[i];
          }
        }
        BasicOp<T> new_op(new_term); 
        new_op.dqn = op1->dqn + op2->dqn;
        bool add = b_contains(new_op);

        if(!add && op1 && op2){
          new_op.resize(rho_basis);
          clock.mark();
  
          if(verbose() > 0)
            dmtk::out() << "ADDING NEW OPERATOR " << new_op.description() << endl;
          new_op = T(0);
          new_operator(new_op, *op1, *op2, pos1, pos2, rho, basis, T(1));
          b_push(new_op);
          if(verbose() > 0)
            dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
        }
      }
    }
  }
}

// ***************************************************************************
// Measurement operators 
// ***************************************************************************
/////////////////////////////////////////////////////////////////////////
// rotate_corr:
// rotate measurement operators 
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::rotate_corr(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *this_hami)
{
  size_t pos1 = BLOCK1;
  size_t pos2 = BLOCK2;
  size_t the_block = BLOCK1;
  size_t site_block = BLOCK2;

  if(position == RIGHT){
    pos1 = BLOCK3;
    pos2 = BLOCK4;
    site_block = BLOCK3;
    the_block = BLOCK4;
  }

  dmtk::Stopwatch clock;
  clock.restart();


  const Hami<T> *ph = &corr;
  if(this_hami) ph = this_hami;
  const Hami<T> &_h = *ph;
  typename Hami<T>::const_iterator titer;

  // Changed here: as in rotate_terms, operators are found through a per-call
  // index over the four blocks (see OpIndex) and whether the new block b
  // already holds an operator is answered from a key set kept next to it,
  // rather than a scan of the block lists per lookup; with the O(n^4)
  // measurement terms of the molecular solver those scans made the final
  // sweep quadratic in the number of operators. The nested calls below add
  // to b, so they drop the set.
  OpIndex index;
  std::set<OpKey> in_b;
  bool in_b_valid = false;
  auto b_contains = [&](const BasicOp<T> &op) {
    if(!in_b_valid){
      in_b.clear();
      for(const BasicOp<T> &o : b) in_b.emplace(o.name(), o.site(), o.internal_site());
      in_b_valid = true;
    }
    return in_b.count(OpKey(op.name(), op.site(), op.internal_site())) != 0;
  };
  auto b_push = [&](const BasicOp<T> &op) {
    b.push_back(op);
    if(in_b_valid) in_b.emplace(op.name(), op.site(), op.internal_site());
  };

// 1a-Look for operators in S.ops involving the other blocks and site block 
// We also rotate individual operators

  if(!this_hami){
    for(titer = ops.begin(); titer != ops.end(); titer++){
      const Term<T>& t = (*titer);

      if(t.size() == 1 && (!_store_products && !t[0].is_hami())) continue;
 
      bool found = false; 
      typename Term<T>::const_iterator oiter;
      int bmask = 0;
      for (BasicOp<T> const & op : t) {
        if (!op.is_hami()) {
          bmask |= mask(block(op.site()));
        } else {
          continue;
        }
      }
      if(position == LEFT && (bmask & MASK_BLOCK2) && !(bmask & MASK_BLOCK1)) found = true;
      if(position == RIGHT && (bmask & MASK_BLOCK3) && !(bmask & MASK_BLOCK4)) found = true;
  
      if(found){
        Term<T> aux_term;
        BasicOp<T> top2;
        for(int i = 0; i < t.size(); i++){
          BasicOp<T> top = t[i];
          if(block(top.site()) == site_block){
            aux_term *= top;
          }
        }
        BasicOp<T> top1;
        if(aux_term.size() == 1) // the piece contains a single operator
          top1 = aux_term[0];
        else
          top1 = aux_term; // the piece contains more than one operator
  
        const BasicOp<T>* _op = indexed_op(top1, index);
  
        if(!_op) {
          dmtk::out() << "ERROR 1a: Operator " << top1.description() << " not found\n";
          continue;
        }
  
        BasicOp<T> real_op(top1); 
        real_op.dqn = _op->dqn;
        bool add = b_contains(real_op);
  
        if(!add && _op){ // If we haven't added it yet, we do it now
          real_op.resize(rho_basis);
  
          if(verbose() > 0)
            dmtk::out() << "ADDING NEW MEAS. OPERATOR 1 " << top1.name() << " " << real_op.name() << " " << real_op.description() << endl;
          clock.mark();
  
          new_operator(real_op, *_op, rho, basis, 1-position);
//          real_op.set_name(t.name().c_str());
          b_push(real_op);
          if(verbose() > 0)
            dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
        }
      }
    }
  } 

// 1b-Look for operators in S.corr involving the other blocks and site block 
// We also rotate individual operators

  for (Term<T> const & t : _h) {
//    if(T(t.coef()) == T(0)) continue;

/*
    if(t.type() == TermType::TERM_PRODUCT && t.size() == 2 && t[1].is_hami()){
      rotate_additive(position, b, basis, rho_basis, t);
    }
*/

    bool found = false;

    int bmask = 0;
    for (BasicOp<T> const & op : t) {
      if (!op.is_hami()) {
        bmask |= mask(block(op.site()));
      } else {
        rotate_terms(position, b, basis, rho_basis, &op.hami());
        rotate_hami(position, b, basis, rho_basis, &op.hami());
        in_b_valid = false;
        continue;
      }
    }
    if(t.size() == 1 && (!_store_products && !t[0].is_hami())) continue;
    if(position == LEFT && (bmask & MASK_BLOCK2) && !(bmask & MASK_BLOCK1)) found = true;
    if(position == RIGHT && (bmask & MASK_BLOCK3) && !(bmask & MASK_BLOCK4)) found = true;

    if(found){
      Term<T> aux_term;
      BasicOp<T> top2;
      for(int i = 0; i < t.size(); i++){
        BasicOp<T> top = t[i];
        if(block(top.site()) == site_block){
          aux_term *= top;
        }
      }
      if(aux_term.size() == t.size()) aux_term.coef() = t.coef();

      BasicOp<T> top1;
      if(aux_term.size() == 1) // the piece contains a single operator
        top1 = aux_term[0];
      else
        top1 = aux_term; // the piece contains more than one operator

      const BasicOp<T>* _op = indexed_op(top1, index);

      if(!_op) {
        dmtk::out() << "ERROR 1b: Operator " << top1.description() << " not found\n";
        continue;
      }

      BasicOp<T> real_op(top1); 
      real_op.dqn = _op->dqn;
      bool add = b_contains(real_op);

      if(!add && _op){ // If we haven't added it yet, we do it now
        real_op.resize(rho_basis);

        if(verbose() > 0)
          dmtk::out() << "ADDING NEW MEAS. OPERATOR 2 " << t.name() << " " << top1.name() << " " << real_op.name() << " " << real_op.description() << endl;
        clock.mark();

        new_operator(real_op, *_op, rho, basis, 1-position);

        b_push(real_op);
        if(verbose() > 0)
          dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
      }
    }
  }

// 2a- We look for terms that involve pieces in old block and pieces
// in another block (not site blocks) and we rotate them 
// We also rotate individual operators

  if(!this_hami){
    for (Term<T> const & t : ops) {
      bool found = false;

      if(t.size() == 1 && (!_store_products && !t[0].is_hami())) continue;
  
      typename Term<T>::const_iterator oiter;
      int bmask = 0;
      for (BasicOp<T> const & op : t) {
        if (!op.is_hami()) {
          bmask |= mask(block(op.site()));
        } else {
          continue;
        }
      }
      if(position == LEFT && (bmask & MASK_BLOCK1) && !(bmask & MASK_BLOCK2)) found = true;
      if(position == RIGHT && (bmask & MASK_BLOCK4) && !(bmask & MASK_BLOCK3)) found = true;
  
      if(found){ // we found a piece of a composite operator
        Term<T> aux_term;
        BasicOp<T> top2;
        for(int i = 0; i < t.size(); i++){
          BasicOp<T> top = t[i];
          if(block(top.site()) == the_block){
            aux_term *= top;
          }
        }
        if(aux_term.size() == t.size()) aux_term.coef() = t.coef();
  
        BasicOp<T> top1;
        if(aux_term.size() == 1) // the piece contains a single operator
          top1 = aux_term[0];
        else
          top1 = aux_term; // the piece contains more than one operator
  
        const BasicOp<T>* _op = indexed_op(top1, index);
  
        if(!_op) {
          dmtk::out() << "ERROR 2a: Operator " << top1.description() << " not found\n";
          continue;
        }
  
        BasicOp<T> new_op(top1); 
        new_op.dqn = _op->dqn;
        bool add = b_contains(new_op);
  
        if(!add && _op){
          new_op.resize(rho_basis);
          clock.mark();
  
          if(verbose() > 0)
            dmtk::out() << "NEW MEAS. OPERATOR " << new_op.description() << endl;
          new_op = T(0);
          new_operator(new_op, *_op, rho, basis, position);
//          new_op.set_name(t.name().c_str());
          if(verbose() > 0)
            dmtk::out() << "NEW MEAS. OPERATOR " << t.name() << " " << new_op.name() << endl;
          b_push(new_op);
          if(verbose() > 0)
            dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
        }
      }
    }
  }

// 2b- We look for terms that involve pieces in old block and pieces
// in another block (not site blocks) and we rotate them 
// We also rotate individual operators

  for(titer = _h.begin(); titer != _h.end(); titer++){
    const Term<T>& t = (*titer);
    if(t.size() == 1 && (!_store_products && !t[0].is_hami())) continue;
    bool found = false;

    int bmask = 0;
    for (BasicOp<T> const & op : t) {
      if (!op.is_hami()) {
        bmask |= mask(block(op.site()));
      } else {
        continue;
      }
    }
    if(position == LEFT && (bmask & MASK_BLOCK1) && !(bmask & MASK_BLOCK2)) found = true;
    if(position == RIGHT && (bmask & MASK_BLOCK4) && !(bmask & MASK_BLOCK3)) found = true;

    if(found){ // we found a piece of a composite operator
      Term<T> aux_term;
      BasicOp<T> top2;
      for(int i = 0; i < t.size(); i++){
        BasicOp<T> top = t[i];
        if(block(top.site()) == the_block){
          aux_term *= top;
        }
      }

      BasicOp<T> top1;
      if(aux_term.size() == 1) // the piece contains a single operator
        top1 = aux_term[0];
      else
        top1 = aux_term; // the piece contains more than one operator

      const BasicOp<T>* _op = indexed_op(top1, index);

      if(!_op) {
        if(_store_products) dmtk::out() << "ERROR 2b: Operator " << top1.description() << " not found\n";
        continue;
      }

      BasicOp<T> new_op(top1); 
      new_op.dqn = _op->dqn;
      bool add = b_contains(new_op);

      if(!add && _op){
        new_op.resize(rho_basis);
        clock.mark();

        if(verbose() > 0)
          dmtk::out() << "NEW MEAS. OPERATOR " << new_op.description() << endl;
        new_op = T(0);
        new_operator(new_op, *_op, rho, basis, position);
//        new_op.set_name(t.name().c_str());
        b_push(new_op);
        if(verbose() > 0)
          dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
      }
    }
  }

// 3- Composite operators for correlations
  if(_store_products){
    for (Term<T> const & t : _h) {
      bool found1 = false;
      bool found2 = false;
      if(t.type() == TermType::TERM_PRODUCT && t.size() >= 2){
        Term<T> aux_term1, aux_term2;
        for(int i = 0; i < t.size(); i++){
          BasicOp<T> top = t[i];
          if(block(top.site()) == pos1){
            aux_term1 *= top;
            found1 = true;
          }
          if(block(top.site()) == pos2){
            aux_term2 *= top;
            found2 = true;
          }
        } 
        if(found1 && found2){ // we have a new composite operator
          BasicOp<T> top1, top2;
          if(aux_term1.size() == 1) // the piece contains a single operator
            top1 = aux_term1[0];
          else
            top1 = aux_term1; // the piece contains more than one operator
          if(aux_term2.size() == 1) // the piece contains a single operator
            top2 = aux_term2[0];
          else
            top2 = aux_term2; // the piece contains more than one operator
          const BasicOp<T>* op1 = indexed_op(top1, index);
          const BasicOp<T>* op2 = indexed_op(top2, index);
  
          if(!op1) {
            dmtk::out() << "ERROR 3: Operator " << top1.description() << " not found\n";
            continue;
          }
          if(!op2) {
            dmtk::out() << "ERROR 3: Operator " << top2.description() << " not found\n";
            continue;
          }
 
          Term<T> new_term; 
          for(int i = 0; i < t.size(); i++){
            BasicOp<T> top = t[i];
            if(block(top.site()) == pos1 || block(top.site()) == pos2){
              new_term *= top;
            }
          }

          BasicOp<T> new_op(new_term); 
          new_op.dqn = op1->dqn + op2->dqn;
          bool add = b_contains(new_op);
  
          if(!add && op1 && op2){
            new_op.resize(rho_basis);
            clock.mark();
    
            if(verbose() > 0)
              dmtk::out() << "ADDING NEW CORR. TERM " << new_op.description() << endl;
            new_op = T(0);
            new_operator(new_op, *op1, *op2, pos1, pos2, rho, basis, T(1));
//            new_op.set_name(new_term.name().c_str());
  
            b_push(new_op);
            if(verbose() > 0)
              dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
          }
        }
      }
    }

  }

}

/////////////////////////////////////////////////////////////////////////
// rotate_hami:
// rotate hamiltonian terms 
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::rotate_hami(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *_this_hami)
{
  if(!_this_hami) _this_hami = &this->h;
  const Hami<T> &this_hami = *_this_hami;

#ifndef USE_PRODUCT_DEFAULT

  if(_this_hami == &this->h){
    BasicOp<T> _h = BasicOp<T>(this->h);
    _h.resize(rho_basis);
    if(position == LEFT){
      new_operator(_h, this->h12, rho, basis);
      b.push_back(_h);
    } else { // position == RIGHT
      new_operator(_h, this->h34, rho, basis);
      b.push_back(_h);
    }
    return;
  }

#endif // USE_PRODUCT_DEFAULT

  const Block<T> *pb1 = this->_b1;
  const Block<T> *pb2 = this->_b2;
  int _mask = MASK_BLOCK1 | MASK_BLOCK2;
  size_t pos1 = BLOCK1;
  size_t pos2 = BLOCK2;

  if(position == RIGHT){
    pb1 = this->_b3;
    pb2 = this->_b4;
    _mask = MASK_BLOCK3 | MASK_BLOCK4;
    pos1 = BLOCK3;
    pos2 = BLOCK4;
  }

  const Block<T> &b1(*pb1);
  const Block<T> &b2(*pb2);

  dmtk::Stopwatch clock;
  clock.restart();

  if(verbose() > 0)
    dmtk::out() << "HAMILTONIAN TERMS: " << this_hami.name() << endl;

  BasicOp<T> _h = BasicOp<T>(this_hami);
  _h.resize(rho_basis);

  bool found2 = false;
  bool found1 = false;
  const BasicOp<T> *local_op2 = b2(BasicOp<T>(this_hami));
  if(local_op2 && local_op2->name() == this_hami.name()){ 
    if(verbose() > 0)
      dmtk::out() << "OPERATOR " << this_hami.name() << "2" << endl;
    clock.mark();
    found2 = true;
    new_operator(_h, *local_op2, rho, basis, RIGHT);
    if(verbose() > 0)
      dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
  } else {
//    if(this_hami.name() == "H")
//      cout << "WARNING: Add local (site) Hamiltonian term\n";
//    else
//      cout << "You can ignore the warning\n";
  }

  const BasicOp<T> *local_op1 = b1(BasicOp<T>(this_hami));
  if(local_op1 && local_op1->name() == this_hami.name()){ 
    if(verbose() > 0)
      dmtk::out() << "OPERATOR " << this_hami.name() << "1" << endl;
    clock.mark();
    found1 = true;
    new_operator(_h, *local_op1, rho, basis, LEFT);
    if(verbose() > 0)
      dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
  } else {
//    if(this_hami.name() == "H")
//      cout << "WARNING: Add local (site) Hamiltonian term\n";
//    else
//      cout << "You can ignore the warning\n";
  }

  for (Term<T> const & t : this_hami) {
//    if(T(t.coef()) == T(0)) continue;
    if(t.type() == TermType::TERM_EXTERN) continue;
    if(t.type() == TermType::TERM_LOCAL) continue;
    if(t.size() == 1 && t[0].name() != this_hami.name()) {
      const BasicOp<T>&top = t[0];

/*
      // Rotate Hamiltonians terms in the original Hamiltonian
      // This is now done inside rotate_terms
      if(top.is_hami()){
        rotate_hami(position, b, basis, rho_basis, &(top.hami()));
        continue;
      }
*/
      if(top.is_hami()) continue;

      size_t ib = block(top.site());
      if((ib == BLOCK2 && position == LEFT) || (ib == BLOCK3 && position == RIGHT)){
        const BasicOp<T>* op= operator()(top);
        bool calc_hc = this_hami.use_hc();
        if (op->is_diagonal()) calc_hc = false;
                                                                                
      
        if(verbose() > 0){ 
          if(!calc_hc)
            dmtk::out() << "TERM " << t.name(true) << endl;
          else
            dmtk::out() << "TERM " << t.name(true) << " + h.c." << endl;
        }
                                                                                
        clock.mark();
        new_operator(_h, *op, rho, basis, 1-position, T(t.coef()));
        if(verbose() > 0)
          dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
      }
//if we don't have a site Hamiltonian we try to find the missing term
      if((ib == BLOCK1 && position == LEFT && !found1) || (ib == BLOCK4 && position == RIGHT && !found2) || (position == LEFT && ib == BLOCK1 && b1.lattice().size() == 1) || (position == RIGHT && ib == BLOCK4 && b2.lattice().size() == 1)){

        const BasicOp<T>* op= operator()(top);
        bool calc_hc = this_hami.use_hc();
        if (op->is_diagonal()) calc_hc = false;
                      
        if(verbose() > 0) {
          if(!calc_hc)
            dmtk::out() << "MISSING TERM " << t.name(true) << endl;
          else
            dmtk::out() << "MISSING TERM " << t.name(true) << " + h.c." << endl;
        }
                                                                                
        clock.mark();
        new_operator(_h, *op, rho, basis, position, T(t.coef()));
        if(verbose() > 0)
          dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
      }
    }
// Multiple products
    if(t.size() >= 2 && t.type() == TermType::TERM_PRODUCT){
      bool found1 = false;
      bool found2 = false;
      Term<T> aux_term1, aux_term2;
      int bmask = 0;
      for(int i = 0; i < t.size(); i++){
        BasicOp<T> top = t[i];
        bmask |= mask(block(top.site()));
        if(block(top.site()) == pos1){
          aux_term1 *= top;
          found1 = true;
        }
        if(block(top.site()) == pos2){
          aux_term2 *= top;
          found2 = true;
        }
      } 
      if((bmask & _mask) == bmask && found1 && found2){ // we have a new composite operator

        BasicOp<T> top1, top2;
        if(aux_term1.size() == 1) // the piece contains a single operator
          top1 = aux_term1[0];
        else
          top1 = aux_term1; // the piece contains more than one operator
        if(aux_term2.size() == 1) // the piece contains a single operator
          top2 = aux_term2[0];
        else
          top2 = aux_term2; // the piece contains more than one operator
        const BasicOp<T>* op1 = operator()(top1);
        const BasicOp<T>* op2 = operator()(top2);

        if(!op1) {
          dmtk::out() << "ERROR : Operator " << top1.description() << " not found\n";
          continue;
        }
        if(!op2) {
          dmtk::out() << "ERROR : Operator " << top2.description() << " not found\n";
          continue;
        }

        bool calc_hc = this_hami.use_hc();
        if (op1->is_diagonal() && op2->is_diagonal()) calc_hc = false;
        if (t.size() > 2) calc_hc = false; // you have to add the h.c. terms explicitly

        if(op1 && op2){
          clock.mark();
 
          if(verbose() > 0){
            if(!calc_hc)
              dmtk::out() << "TERM " << t.description() << endl;
            else
              dmtk::out() << "TERM " << t.description() << " + h.c." << endl;
          }
          new_operator(_h, *op1, *op2, pos1, pos2, rho, basis, T(t.coef()), calc_hc);
         
          if(verbose() > 0) 
            dmtk::out() << "Lap: " << dmtk::format_elapsed(clock.lap_elapsed()) << endl;
        }
      }
    }
  }

  b.push_back(_h);
}

/////////////////////////////////////////////////////////////////////////
// measure:
// measure physical quantities (mean value of the observables)
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::measure(const VectorState<T> *v)
{

  using namespace std;
  DMTKglobals<T> &globals = get_globals(T(0));

  const VectorState<T> *pv = nullptr;
  if(!v)
   pv = &gs;
  else
   pv = v;

  VectorState<T> res(*pv);
  VectorState<T> aux(*pv);
  VectorState<T> aux_res(*pv);
  T x = product(*pv,*pv);
  T x0 = product(*pv,*pv);

  // Changed here: the measured terms' operators are found through a per-call
  // index over the four blocks (see OpIndex), which do not change during the
  // measurement, rather than a scan of the block lists per lookup.
  OpIndex index;


//////////////////////////////////////////////////////////////////////////
  for (Term<T> & real_t : corr) {
    Term<T> t = real_t;

    if(t.size() == 1 && t[0].is_hami()) {
      res = T(0);
      res = product_default(*this, *pv, &t[0].hami());
      T val = product(*pv,res)/x;  
      real_t.set_value(val);
      continue;
    }
    if(t.size() == 1 && !t[0].is_hami()) {
      bool info;
      // apply the Term's coefficient, exactly as the 2-op and n-op
      // branches below pass t.coef() into product(): without it a
      // weighted single-operator Term (e.g. a correlation's on-site
      // A_i*B_i entry carrying its 1/multiplicity weight) is measured
      // with weight 1.
      T val = measure_operator(t[0], v, info);
      if(info) real_t.set_value(T(t.coef())*val);
      continue;
    }

// We measure terms with 2 operators first 
    if(t.size() == 2){
      BasicOp<T> top2 = t[0];
      BasicOp<T> top1 = t[1];
  
      size_t ib1 = block(top1.site());
      size_t ib2 = block(top2.site());
 
      if(ib1 != BLOCK_NONE && ib2 != BLOCK_NONE){
        if(ib1 != ib2){
            
          const BasicOp<T>* op1 = indexed_op(top1, index);
          const BasicOp<T>* op2 = indexed_op(top2, index);
 
          if(op1 && op2){
            if(op1->name() != top1.name() || op2->name() != top2.name()) 
                 continue; 
 
            res = T(0);
            product(*op1, *op2, *pv, res, ib1, ib2, T(t.coef()));
            T val = product(*pv,res)/x;
 
/* 
            if(h.lattice().type() == LATTICE_1D)
              cout << setprecision(10) << t.coef() << "*" << top1.name().c_str() << "(" << h.lattice().x(top1.site()) << "," << top1.internal_site() << ")" << top2.name().c_str() << "(" << h.lattice().x(top2.site()) << "," << top2.internal_site() << ") = " << val << endl;
            if(h.lattice().type() == LATTICE_2D)
              cout << setprecision(10) << t.coef() << "*" << top1.name().c_str() << "(" << h.lattice().x(top1.site()) << "," << h.lattice().y(top1.site()) << "," << top1.internal_site() << ")" << top2.name().c_str() << "(" << h.lattice().x(top2.site()) << "," << h.lattice().y(top2.site()) << "," << top2.internal_site() << ") = " << val << endl;
*/
  
//            cout << setprecision(10) << t.name(true) << " = " << val << endl;
  
            real_t.set_value(val);
          }
  
        } else {

          BasicOp<T> aux_op(t);
          const BasicOp<T>* op = indexed_op(aux_op, index);

 
          if(op){
            res = T(0);
            product(*op, *pv, res, ib1, T(t.coef()));
            T val = product(*pv,res)/x;

/*  
            if(h.lattice().type() == LATTICE_1D)
              cout << setprecision(10) << t.coef() << "*" << top1.name().c_str() << "(" << h.lattice().x(top1.site()) << "," << top1.internal_site() << ")" << top2.name().c_str() << "(" << h.lattice().x(top2.site()) << "," << top2.internal_site() << ") = " << val << endl;
            if(h.lattice().type() == LATTICE_2D)
              cout << setprecision(10) << t.coef() << "*" << top1.name().c_str() << "(" << h.lattice().x(top1.site()) << "," << h.lattice().y(top1.site()) << "," << top1.internal_site() << ")" << top2.name().c_str() << "(" << h.lattice().x(top2.site()) << "," << h.lattice().y(top2.site()) << "," << top2.internal_site() << ") = " << val << endl;
*/
  
//            cout << setprecision(10) << t.name(true) << " = " << val << endl;

            real_t.set_value(val);
          }
        }
      }
  
    } else {
// We measure terms with arbitrary number of operators
      t = t.reorder();
      Term<T> aux_term[4];
      res = T(0);
      for(int ib = 0; ib < 4; ib++) aux_term[ib].clear();
      bool found = true;
      [[maybe_unused]] int bmask = 0;
      for(int i = 0; i < t.size(); i++){
        BasicOp<T> top = t[i].internals();
        if(block(top.site()) == BLOCK_NONE){
          found = false;
          break;
        }
        for(int ib = 1; ib <= 4; ib++){
          if(block(top.site()) == static_cast<size_t>(ib)){
            aux_term[ib-1] *= top;
            bmask |= mask(static_cast<size_t>(ib));
          }
        }
      }
  
      if(!found) continue;

      found = true;
      const BasicOp<T> *aux_op[4];
      aux_op[0] = aux_op[1] = aux_op[2] = aux_op[3] = nullptr;
      for(int ib = 3; ib >= 0 ; ib--){
        if(aux_term[ib].size() != 0){
          BasicOp<T> top;
          if(aux_term[ib].size() == 1)
            top = aux_term[ib][0].internals();
          else
            top = aux_term[ib];

          aux_op[ib] = indexed_op(top, index);
          if(!aux_op[ib]){ // && _verbose >= 0) {
            dmtk::out() << "ERROR : Operator " << top.description() << " not found\n";
            found = false;
            break;
          }
        }
      }

      if(!found) continue;
  
      if (bmask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
        product(*(aux_op[0]), *(aux_op[1]), *(aux_op[2]), *(aux_op[3]), *pv, res, BLOCK1, BLOCK2, BLOCK3, BLOCK4, T(t.coef()), false, &globals);
        T val = product(*pv,res)/sqrt(x)/sqrt(x0);
        real_t.set_value(val);
        if(_verbose >= 0) dmtk::out() << setprecision(10) << t.description() << " = " << val << endl;
      }
      else if(bmask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3)){
        product(*(aux_op[0]), *(aux_op[1]), *(aux_op[2]), *pv, res, BLOCK1, BLOCK2, BLOCK3, T(t.coef()), false, &globals);
        T val = product(*pv,res)/sqrt(x)/sqrt(x0);
        real_t.set_value(val);
        if(_verbose >= 0) dmtk::out() << setprecision(10) << t.description() << " = " << val << endl;
      }
      else if(bmask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4)){
        product(*(aux_op[0]), *(aux_op[1]), *(aux_op[3]), *pv, res, BLOCK1, BLOCK2, BLOCK4, T(t.coef()), false, &globals);
        T val = product(*pv,res)/sqrt(x)/sqrt(x0);
        real_t.set_value(val);
        if(_verbose >= 0) dmtk::out() << setprecision(10) << t.description() << " = " << val << endl;
      }
      else if(bmask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4)){
        product(*(aux_op[0]), *(aux_op[2]), *(aux_op[3]), *pv, res, BLOCK1, BLOCK3, BLOCK4, T(t.coef()), false, &globals);
        T val = product(*pv,res)/sqrt(x)/sqrt(x0);
        real_t.set_value(val);
        if(_verbose >= 0) dmtk::out() << setprecision(10) << t.description() << " = " << val << endl;
      }
      else if(bmask == (MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
        product(*(aux_op[1]), *(aux_op[2]), *(aux_op[3]), *pv, res, BLOCK2, BLOCK3, BLOCK4, T(t.coef()), false, &globals);
        T val = product(*pv,res)/sqrt(x)/sqrt(x0);
        real_t.set_value(val);
        if(_verbose >= 0) dmtk::out() << setprecision(10) << t.description() << " = " << val << endl;
      }
      else
      {
        aux = aux_res = *pv;
        for(int ib = 3; ib >= 0 ; ib--){
          if(aux_op[ib]){
            res.set_qn_mask(aux.qn()+aux_op[ib]->dqn, res.qn_mask());
// FIXME
            res.resize(aux.b1(),aux.b2(),aux.b3(),aux.b4());
            res = T(0);

            product(*aux_op[ib], aux, res, (size_t)(ib+1));
            aux = res;
          }
        }
        T val = product(*pv,res)/sqrt(x)/sqrt(x0);
        real_t.set_value(val*t.coef());
        if(_verbose >= 0) dmtk::out() << setprecision(10) << t.description() << " = " << real_t.value() << endl;
      }
    }
  }
}

template<class T>
T
System<T>::measure_operator(const BasicOp<T> &top, const VectorState<T> *v, bool &info)
{
  using namespace std;
  info = false;

  const VectorState<T> *pv = nullptr;
  if(!v)
   pv = &gs;
  else
   pv = v;

  VectorState<T> res(*pv);
  VectorState<T> aux(*pv);
  T x = sqrt(product(*pv,*pv));
  T x0= sqrt(product(gs,gs));

  size_t ib = block(top.site());

  const BasicOp<T>* op = operator()(top);
  if(!op) {
#ifdef WITH_WARNINGS
// Your are likely to be using _sotre_products == FALSE
    dmtk::out() << "WARNING: Operator " << top.name() << "not found\n";
#endif
    return T(0);
  }

  if(op->name() != top.name()) return T(0);

  info = true;
  res.set_qn_mask(aux.qn()+op->dqn, _grand_canonical);
  res.resize(aux.b1(),aux.b2(),aux.b3(),aux.b4());
  res = T(0);

  product(*op, aux, res, ib);

  T val = product(gs,res)/x/x0;

/*
  if(h.lattice().type() == LATTICE_2D)
    cout << top.name().c_str() << "(" << h.lattice().x(top.site()) << "," << h.lattice().y(top.site()) << "," << top.internal_site() << ")";
  if(h.lattice().type() == LATTICE_1D)
    cout << top.name().c_str() << "(" << h.lattice().x(top.site()) << "," << top.internal_site() << ")";
  cout << setprecision(10) << " = " << val << endl;
*/
  return val;
}


//////////////////////////////////////////////////////////////////////
// product_default:
// This version of product uses the four blocks
/////////////////////////////////////////////////////////////////////
// Original version of product, using 4 blocks
/////////////////////////////////////////////////////////////////////
template<class T>
VectorState<T>
product_default(const System<T>& ss, const VectorState<T>& vv, const Hami<T>* hami, bool only_local)
{
  VectorState<T> res(vv);
  res = T(0);

  const Hami<T> *phami;
  if(hami == nullptr)
    phami = &ss.h;
  else
    phami = hami;
  const Hami<T> &h = *phami;

/////////////////////////////////////////////////////////////
  VectorState<T> aux(res), aux_res(res);
  for (Term<T> const & t : h) {
//    if(T(t.coef()) == T(0)) continue;
    if(t.type() == TermType::TERM_EXTERN && !ss.apply_extern()) continue;
    if(t.size() > 2){  // you have to include h.c. terms explicitly!
      Term<T> aux_term[4];

      if(ss.verbose() > 1)
        dmtk::out() << "PRODUCT TERM " << t.description() << '\n';

      bool found = true;
      int bmask = 0;
      for(int ib = 0; ib < 4; ib++) aux_term[ib].clear();
      for(int i = 0; i < t.size(); i++){
        BasicOp<T> top = t[i];
        if(ss.block(top.site()) == BLOCK_NONE){
          found = false;
          break;
        }
        bmask |= mask(ss.block(top.site()));
        for(int ib = 1; ib <= 4; ib++){
          if(ss.block(top.site()) == static_cast<size_t>(ib)){
            aux_term[ib-1] *= top;
          }
        }
      } 
 
      if(bmask != MASK_BLOCK1 && bmask != MASK_BLOCK2 && bmask != MASK_BLOCK3 && bmask != MASK_BLOCK4 && found){ 
        aux = aux_res = vv;
        for(int ib = 3; ib >= 0 ; ib--){
          if(aux_term[ib].size() != 0){
            BasicOp<T> top;
            if(aux_term[ib].size() == 1)
              top = aux_term[ib][0];
            else
              top = aux_term[ib];
    
            const BasicOp<T>* _op = ss(top);
            if(!_op) {
              dmtk::out() << "ERROR : Operator " << top.description() << " not found\n";
            }
    
            if(_op){
              if(ss.verbose() > 1)
                dmtk::out() << "APPLYING " << top.description() << endl;
    
              aux_res.set_qn_mask(aux.qn()+_op->dqn, ss.qn_mask());
              aux_res.resize(aux.b1(),aux.b2(),aux.b3(),aux.b4());
              aux_res = T(0);
    
              product(*_op, aux, aux_res, static_cast<size_t>(ib+1));
              aux = aux_res;
            }
          }
        }
        res += aux_res * t.coef(); 
      }
    }  
    if(t.size() == 2 && t.type() == TermType::TERM_PRODUCT){
      const BasicOp<T>& top1 = t[0];
      const BasicOp<T>& top2 = t[1];

      bool calc_hc = h.use_hc();
      if (top1.is_diagonal() && top2.is_diagonal()) calc_hc = false;

      size_t b1 = ss.block(top1.site());
      size_t b2 = ss.block(top2.site());

      if(only_local && mask(b1,b2) == (MASK_BLOCK3|MASK_BLOCK4)) continue;
      if(only_local && mask(b1,b2) == (MASK_BLOCK1|MASK_BLOCK2)) continue;

      if(b1 != BLOCK_NONE && b2 != BLOCK_NONE){
        if(b1 != b2){
          const BasicOp<T>* op1 = ss(top1);
          const BasicOp<T>* op2 = ss(top2);

          if(ss.verbose() > 1)
            dmtk::out() << "PRODUCT TERM " << t.description() << '\n';
//" " << top1.description() << " " << top1.site() << " " << top2.site() << " " << b1 << " " << b2 << endl;

          if(op1 && op2) {
            if(ss.verbose() > 0) {
              COUT_PRODUCT(t);
            }
            product(*op2, *op1, vv, res, b2, b1, T(t.coef()), calc_hc);
          }

        } else if(t.type() == TermType::TERM_EXTERN && ss.apply_extern()) {

          if(only_local && (b1 == BLOCK1 || b1 == BLOCK4)) continue;

          BasicOp<T> aux_op(t);
          const BasicOp<T>* op = ss(aux_op);

          if(ss.verbose() > 1) {
            if(calc_hc) 
              dmtk::out() << "PRODUCT EXTERN " << t.coef() << "*" << top1.name() << "(" << top1.site() << ")" << top2.name() << "(" << top2.site() << ") + h.c." << '\n';
            else 
              dmtk::out() << "PRODUCT EXTERN " << t.coef() << "*" << top1.name() << "(" << top1.site() << ")" << top2.name() << "(" << top2.site() << ")" << '\n';
          }
  
          if(op){
            product(*op, vv, res, b1, T(t.coef()));
            // obviously, op1 has to preserve the quantum numbers,
            // or we are working in the grand canonical 
            if(calc_hc) 
              product(*op, vv, res, b1, T(t.coef()), calc_hc);
          }

        }
      }
    }
    else if(t.size() == 1){

      const BasicOp<T>& top = t[0];
      size_t b1 = ss.block(top.site());

      if(top.is_hami() && top.name() != h.name()){
        res += product_default(ss, vv, &top.hami());
        continue;
      }

//      if(((b1 == BLOCK1 && ss.leftblock.lattice().size() != 1) || (b1 == BLOCK4 && ss.rightblock.lattice().size() != 1)) && t.type() != TermType::TERM_LOCAL && t.type() != TermType::TERM_EXTERN) continue;
//      if(((b1 == BLOCK1 && ss.leftblock.lattice().size() != 1) || (b1 == BLOCK4 && ss.rightblock.lattice().size() != 1))) continue;
//      if((b1 == BLOCK1 || b1 == BLOCK4) && t.type() != TermType::TERM_LOCAL && t.type() != TermType::TERM_EXTERN) continue;


      if(((b1 == BLOCK1 && ss.leftblock.lattice().size() != 1) || (b1 == BLOCK4 && ss.rightblock.lattice().size() != 1))) continue;

      if(b1 != BLOCK_NONE && top.name() != h.name() && !top.is_hami()){
        const BasicOp<T>* op1 = ss(top);
        bool calc_hc = h.use_hc();
        if (op1 && op1->is_diagonal()) calc_hc = false;

        if(ss.verbose() > 1) {
          if(calc_hc)
            dmtk::out() << "PRODUCT LOCAL " << t.coef() << "*" << top.name() << "(" << top.site() << ") + h.c." << '\n';
          else
            dmtk::out() << "PRODUCT LOCAL " << t.coef() << "*" << top.name() << "(" << top.site() << ")" << '\n';
        }

        if(op1){
          product(*op1, vv, res, b1, T(t.coef()));
          // obviously, op1 has to preserve the quantum numbers,
          // or we are working in the grand canonical
          if(calc_hc) {
            product(*op1, vv, res, b1, T(t.coef()), calc_hc);
          }
        }
      }
    }
  }

  const char *name = h.name().c_str();
  const BasicOp<T>* local_op4 = ss._b4->operator()(H<T>().set_name(name));
  if(local_op4 && local_op4->name() == h.name()){ 
    if(ss.verbose() > 1) 
      dmtk::out() << "PRODUCT " << name << "4\n";
/*
typename dmtk::BMatrix<T>::const_iterator biter;
for(int i = 0; i < QN::max_index(); i++) cout << local_op4->dqn[i].get_twice();
cout << endl;
for(biter = local_op4->begin(); biter != local_op4->end(); biter++){
  const SubMatrix<T> &bm = *biter;
  for(int i = 0; i < QN::max_index(); i++) cout << bm.qn()[i].get_twice();
  cout << endl;
}
*/
/*
  for(int i = 0; i < bm.cols(); i++)
  for(int j = 0; j < bm.rows(); j++)
    cout << i << " " << j << " " << bm(i,j) << endl;
*/
    product(*local_op4, vv, res, BLOCK4);
  } else {
//    if(h.name() == "H")
//      cout << "WARNING: Add local (site) Hamiltonian term\n";
//    else
//      cout << "You can ignore the warning\n";
  }
  const BasicOp<T>* local_op3 = ss._b3->operator()(H<T>().set_name(name));
  if(local_op3 && local_op3->name() == h.name()){ 
    if(ss.verbose() > 1) 
      dmtk::out() << "PRODUCT " << name << "3\n";
    product(*local_op3, vv, res, BLOCK3);
  } else {
//    if(h.name() == "H")
//      cout << "WARNING: Add local (site) Hamiltonian term\n";
//    else
//      cout << "You can ignore the warning\n";
  }
  const BasicOp<T>* local_op2 = ss._b2->operator()(H<T>().set_name(name));
  if(local_op2 && local_op2->name() == h.name()){ 
    if(ss.verbose() > 1) 
      dmtk::out() << "PRODUCT " << name << "2\n";
    product(*local_op2, vv, res, BLOCK2);
  } else {
//    if(h.name() == "H")
//      cout << "WARNING: Add local (site) Hamiltonian term\n";
//    else
//      cout << "You can ignore the warning\n";
  }
  const BasicOp<T>* local_op1 = ss._b1->operator()(H<T>().set_name(name));
  if(local_op1 && local_op1->name() == h.name()){ 
    if(ss.verbose() > 1) 
      dmtk::out() << "PRODUCT " << name << "1\n";
    product(*local_op1, vv, res, BLOCK1);
  } else {
//    if(h.name() == "H")
//      cout << "WARNING: Add local (site) Hamiltonian term\n";
//    else
//      cout << "You can ignore the warning\n";
  }

  return res;
}

/////////////////////////////////////////////////////////////////////////
// product: 
/////////////////////////////////////////////////////////////////////////
#ifdef USE_PRODUCT_DEFAULT
template<class T>
VectorState<T>
product(const System<T>& ss, const VectorState<T>& vv)
{
  VectorState<T> res(vv);
  res = T(0);
  res = product_default(ss, vv);

  if(ss.project()){
    for(int i = 0; i < ss._project_states.size(); i++){
      VectorState<T> aux = *ss._project_states[i];
      T p = product(aux,vv);
      res = res + aux * p * T(100000);
    }
  }

  return res;
}

template<class T>
VectorState<T>
System<T>::hamiltonian_product(const VectorState<T>& vv)
{
  return product(static_cast<const System<T>&>(*this), vv);
}
#else // USE_PRODUCT_DEFAULT
////////////////////////////////////////////////////////////////////////
// build_composite:
////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::build_composite_operators(const Hami<T> *this_hami)
{
  dmtk::Stopwatch aux_timer;
  aux_timer.restart();

  if(!this_hami) aux_terms.clear();

  const Hami<T> *ph = &h;
  if(this_hami) ph = this_hami;
  const Hami<T> &_h = *ph;

  // Changed here: operators are found through a per-call index rather than
  // a scan of the block lists per lookup, terms are read by reference and
  // moved into aux_terms, and each factor's block is resolved once, before
  // the AuxTerm is built, so that the terms the mask test below drops are
  // never copied.
  OpIndex index;
  std::vector<size_t> where;

  typename Hami<T>::const_iterator hiter;
  for(hiter = _h.begin(); hiter != _h.end(); hiter++){
    const Term<T> &t = *hiter;
    if(_use_coef_tol && fabs(dmtk::real(t.coef())) <= _coef_tol) continue;

    if(t.size() == 1 && t[0].is_hami()) {
      build_composite_operators(&t[0].hami());
    }

    if(t.size() >= 2 && ((apply_hami() && t.type() == TERM_PRODUCT) || (apply_extern() && t.type() == TERM_EXTERN))){
      int bmask = 0;
      bool found = true;
      where.assign(t.size(), BLOCK_NONE);
      for(int i = 0; i < t.size(); i++){
        where[i] = block(t[i].site());
        if(where[i] == BLOCK_NONE){
          found = false;
//          cout << "WARNING: *** Operator " << top.description() << " " << top.site() << " not found\n";
          break;
        }
        bmask |= mask(where[i]);
      }
      if(!(found && bmask != MASK_BLOCK1 && bmask != MASK_BLOCK2 && bmask != MASK_BLOCK3 && bmask != MASK_BLOCK4 && ((t.type() != TERM_EXTERN && bmask != (MASK_BLOCK1|MASK_BLOCK2) && bmask != (MASK_BLOCK2|MASK_BLOCK3) && bmask != (MASK_BLOCK3|MASK_BLOCK4)) || t.type() == TERM_EXTERN || this_hami != NULL))) continue;

      AuxTerm<T> new_term;
      new_term.t = (*hiter);
      for(int i = 0; i < t.size(); i++){
        const BasicOp<T> top = t[i].internals();
        new_term.cterm[where[i]-1] *= top;
        new_term.pos[where[i]-1] = i;
      }
      new_term.mask = bmask;
      new_term.nblocks = 0;
      new_term.nblocks += ((bmask & MASK_BLOCK1) == MASK_BLOCK1 ? 1 : 0);
      new_term.nblocks += ((bmask & MASK_BLOCK2) == MASK_BLOCK2 ? 1 : 0);
      new_term.nblocks += ((bmask & MASK_BLOCK3) == MASK_BLOCK3 ? 1 : 0);
      new_term.nblocks += ((bmask & MASK_BLOCK4) == MASK_BLOCK4 ? 1 : 0);
// WARNING: FIXME
if(bmask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4) || bmask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4)) new_term.nblocks = 3;
if(bmask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3) || bmask == (MASK_BLOCK4|MASK_BLOCK3|MASK_BLOCK2)) new_term.nblocks = 3;

      bool found_op = true;
      { // the mask test that guarded this block is the continue above
        for(int ib = 3; ib >= 0 ; ib--){
          if(new_term.cterm[ib].size() != 0){
            if(new_term.cterm[ib].size() == 1)
              new_term.top[ib] = new_term.cterm[ib][0];
            else
              new_term.top[ib] = new_term.cterm[ib];

            new_term.op[ib] = indexed_op(new_term.top[ib], index);
            if(!new_term.op[ib]) {
              dmtk::out() << "ERROR : Operator " << new_term.top[ib].description() << " not found\n";
              found_op = false;
            }
          }
        }

        if(found_op) aux_terms.push_back(std::move(new_term));
      }
    }
  }

  // Changed here: a leader collects the terms it can be summed with from a
  // table keyed on what the match below compares (mask, diagonality, total
  // charge change, and the operators in the leader's blocks b1 and b3)
  // instead of comparing it with every later term, which made each
  // iteration quadratic in the number of terms (n^4 for molecular
  // Hamiltonians). A term enters the table once for each (b1, b3) choice a
  // leader with its mask can make. Every matching later term joins its first
  // leader in list order. The scan this replaces erased each term it merged
  // and then stepped over the one after it, which was left to a later leader
  // or to a group of its own; when the merged term was the last one, the
  // step went past the end of the list, which libc++ and libstdc++ turn into
  // a second pass from the front. The summed operator is unchanged; only the
  // grouping, and so the rounding of the sums, can differ. Runs far from
  // convergence can amplify that: at 8 and 10 orbitals with 32 states and
  // two sweeps, the energy moves by 1e-3 and 1e-2, no more than the old
  // code moves when the integrals change by 1e-13.
  using QNKey = std::vector<short>;
  auto qn_key = [](const QN &q) {
    QNKey key;
    for(int i = 0; i < int(QN::QN_LAST); i++)
      if(QN::get_qn_mask() & (1 << i)) key.push_back(q[i].get_twice());
    return key;
  };
  using OpId = std::tuple<QNKey, std::string, int, int, int>;
  auto op_id = [&](const BasicOp<T> *op) {
    return OpId(qn_key(op->dqn), op->name(), op->site(), op->internal_site(), op->label());
  };
  using GroupKey = std::tuple<int, bool, QNKey, size_t, size_t, OpId, OpId>;
  auto variants = [](int m) -> std::vector<std::pair<size_t, size_t> > {
    switch(m){
      case (MASK_BLOCK1|MASK_BLOCK4): return {{BLOCK4, BLOCK_NONE}, {BLOCK1, BLOCK_NONE}};
      case (MASK_BLOCK2|MASK_BLOCK4): return {{BLOCK2, BLOCK_NONE}};
      case (MASK_BLOCK1|MASK_BLOCK3): return {{BLOCK3, BLOCK_NONE}};
      case (MASK_BLOCK1|MASK_BLOCK2): return {{BLOCK2, BLOCK_NONE}};
      case (MASK_BLOCK3|MASK_BLOCK4): return {{BLOCK3, BLOCK_NONE}};
      case (MASK_BLOCK2|MASK_BLOCK3): return {{BLOCK2, BLOCK_NONE}};
      case (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3): return {{BLOCK2, BLOCK3}};
      case (MASK_BLOCK4|MASK_BLOCK2|MASK_BLOCK3): return {{BLOCK2, BLOCK3}};
      case (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4): return {{BLOCK2, BLOCK4}, {BLOCK2, BLOCK1}};
      case (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4): return {{BLOCK3, BLOCK4}, {BLOCK3, BLOCK1}};
    }
    return {};
  };
  auto group_key = [&](const AuxTerm<T> &a, size_t b1, size_t b3) {
    return GroupKey(a.mask, a.t.is_diagonal(), qn_key(a.t.dqn()), b1, b3, op_id(a.op[int(b1)-1]),
                    b3 == BLOCK_NONE ? OpId() : op_id(a.op[int(b3)-1]));
  };
  typedef typename std::list<AuxTerm<T> >::iterator aux_iterator;
  std::map<GroupKey, std::vector<aux_iterator> > groups;
  if(use_composite())
    for(aux_iterator a = aux_terms.begin(); a != aux_terms.end(); a++)
      if(!a->done)
        for(auto const &v : variants(a->mask))
          groups[group_key(*a, v.first, v.second)].push_back(a);

//  cout << "COMPOSITE TERMS: " << aux_terms.size() << endl;
  typename std::list<AuxTerm<T> >::iterator titer;
  for(titer = aux_terms.begin(); titer != aux_terms.end(); titer++){
    AuxTerm<T>& auxt = (*titer);

    if(!auxt.done){

      if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK4)){
        auxt.b1 = BLOCK4, auxt.b2 = BLOCK1;
        if(operator[](BLOCK1).lattice().size() < operator[](BLOCK4).lattice().size()){
          auxt.b1 = BLOCK1; auxt.b2 = BLOCK4;
        }
        if(auxt.cterm[0].size() == 1) { auxt.b1 = BLOCK1; auxt.b2 = BLOCK4; }
        if(auxt.cterm[3].size() == 1) { auxt.b1 = BLOCK4; auxt.b2 = BLOCK1; }
      } else if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK4)){
        auxt.b1 = BLOCK2; auxt.b2 = BLOCK4;
      } else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK3)){
        auxt.b1 = BLOCK3; auxt.b2 = BLOCK1;
//------------------------------------------------------------------
      } else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2)){
        auxt.b1 = BLOCK2; auxt.b2 = BLOCK1;
      } else if(auxt.mask == (MASK_BLOCK3|MASK_BLOCK4)){
        auxt.b1 = BLOCK3; auxt.b2 = BLOCK4;
      } else if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK3)){
        auxt.b1 = BLOCK2; auxt.b2 = BLOCK3;
      } else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3)){
        auxt.b1 = BLOCK2; auxt.b2 = BLOCK1; auxt.b3 = BLOCK3;
      } else if(auxt.mask == (MASK_BLOCK4|MASK_BLOCK2|MASK_BLOCK3)){
        auxt.b1 = BLOCK2; auxt.b2 = BLOCK4; auxt.b3 = BLOCK3;
      } else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4)){
        auxt.b1 = BLOCK2; auxt.b2 = BLOCK1; auxt.b3 = BLOCK4;
        if(operator[](BLOCK1).lattice().size() < operator[](BLOCK4).lattice().size()){
          auxt.b3 = BLOCK1; auxt.b2 = BLOCK4;
        }
        if(auxt.cterm[3].size() > 1) { auxt.b3 = BLOCK1; auxt.b2 = BLOCK4; }
        if(auxt.cterm[0].size() > 1) { auxt.b3 = BLOCK4; auxt.b2 = BLOCK1; }
      } else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4)){
        auxt.b1 = BLOCK3; auxt.b2 = BLOCK1; auxt.b3 = BLOCK4;
        if(operator[](BLOCK1).lattice().size() < operator[](BLOCK4).lattice().size()){
          auxt.b3 = BLOCK1; auxt.b2 = BLOCK4;
        }
        if(auxt.cterm[3].size() > 1) { auxt.b3 = BLOCK1; auxt.b2 = BLOCK4; }
        if(auxt.cterm[0].size() > 1) { auxt.b3 = BLOCK4; auxt.b2 = BLOCK1; }
      } else
        continue;

      if(!use_composite()) continue;

      auxt.ref_op[0] = (auxt.op[int(auxt.b1)-1]);
      auxt.sum_op = *(auxt.op[int(auxt.b2)-1]);
      if(auxt.nblocks >= 3) auxt.ref_op[1] = (auxt.op[int(auxt.b3)-1]);

      typename BMatrix<T>::iterator biter;
      for(biter = auxt.sum_op.begin(); biter != auxt.sum_op.end(); biter++){
        SubMatrix<T> &sm = (*biter);
        sm *= auxt.t.coef();
      }
      auxt.done = true;

//      if(verbose() > 0)
      {
        if(auxt.nblocks >=3 ) dmtk::out() << " " << auxt.ref_op[1]->description();
        dmtk::out() << endl;
      }

      auto const group = groups.find(group_key(auxt, auxt.b1, auxt.nblocks >= 3 ? auxt.b3 : BLOCK_NONE));
      if(group == groups.end()) continue;
      int nterm = 1;
      for(aux_iterator titer2 : group->second){
        AuxTerm<T>& auxt2 = (*titer2);
        if(auxt2.done) continue;
        const Term<T>& t2 = auxt2.t;
        auxt2.done = true;
        auxt2.mask = 0;
        const BMatrix<T>& bm2 = *auxt2.op[int(auxt.b2)-1];

//////////////////////////////////////////////
// If this is the second term in the sum,
// initialize sum_op to include all the blocks
//////////////////////////////////////////////
        if(nterm == 1){
          BasicOp<T> aux_sum;
          aux_sum = auxt.sum_op.internals();
          const Block<T> *_block;
          valarray<const Block<T>* > b(4);
          b[0] = _b1;
          b[1] = _b2;
          b[2] = _b3;
          b[3] = _b4;
          _block = b[int(auxt.b2)-1];
          aux_sum.resize(_block->basis().subspaces());
          typename BMatrix<T>::iterator biter1;
          for(biter1 = auxt.sum_op.begin(); biter1 != auxt.sum_op.end(); biter1++){
            SubMatrix<T> &_sm1 = (*biter1);
            SubMatrix<T> *_sm2 = aux_sum.block(_sm1.qn());
            *_sm2 = _sm1;
          }
          auxt.sum_op = aux_sum;
        }
//////////////////////////////////////////////

        nterm++;
        typename BMatrix<T>::const_iterator biter2;
        for(biter2 = bm2.begin(); biter2 != bm2.end(); biter2++){
          const SubMatrix<T> &_sm2 = (*biter2);
          SubMatrix<T> *_sm1 = auxt.sum_op.block(_sm2.qn());
          if(_sm1){
            Matrix<T> sm2 = _sm2;
            Matrix<T> &sm1 = *_sm1;
            sm2 *= t2.coef();
            sm1 += sm2;
          } else {
            dmtk::out() << "WARNING: Block not found " << auxt.t.description() << " " << auxt2.t.description() << endl;
          }
        }
      }
    }
  }
  // The merged terms leave the list once no table entry refers to them.
  aux_terms.remove_if([](const AuxTerm<T> &a) { return a.done && a.mask == 0; });
//  cout << "COMPOSITE TERMS: " << aux_terms.size() << endl;
//  cout << "Build composite time: " << dmtk::format_elapsed(aux_timer.elapsed()) << endl;
}


template<class T>
VectorState<T>
product_composite(System<T>& ss, const VectorState<T>& vv, int mask_hc)
{
  VectorState<T> res(vv);
  res = T(0);
  dmtk::Stopwatch aux_timer;
  DMTKglobals<T> globals = get_globals(T(0));

  aux_timer.restart();
  VectorState<T> v12 = vv.condense(MASK_BLOCK1|MASK_BLOCK2);
  VectorState<T> v1234 = v12.condense(MASK_BLOCK3|MASK_BLOCK4);
  VectorState<T> res1234(v1234);
  if(ss.verbose() > 0)
    dmtk::out() << "CONDENSE TIME= " << dmtk::format_elapsed(aux_timer.lap_elapsed()) << endl;
  aux_timer.mark();
  res1234 = T(0);
  if(ss.verbose() > 0)
    dmtk::out() << "PRODUCT H12\n";
  product(ss.h12, v1234, res1234, BLOCK1, T(1), false, &globals);
  if(ss.verbose() > 0)
    dmtk::out() << dmtk::format_elapsed(aux_timer.lap_elapsed()) << endl;
  aux_timer.mark();
  if(ss.verbose() > 0)
    dmtk::out() << "PRODUCT H34\n";
  product(ss.h34, v1234, res1234, BLOCK4, T(1), false, &globals);
  if(ss.verbose() > 0)
    dmtk::out() << dmtk::format_elapsed(aux_timer.lap_elapsed()) << endl;
  aux_timer.mark();
  VectorState<T> res12 = res1234.decondense(MASK_BLOCK3|MASK_BLOCK4, v12);
  res = res12.decondense(MASK_BLOCK1|MASK_BLOCK2, vv);
  if(ss.verbose() > 0)
    dmtk::out() << "DECONDENSE TIME= " << dmtk::format_elapsed(aux_timer.lap_elapsed()) << endl;

// Blocks 2 and 3

  if(ss.verbose() > 0)
    dmtk::out() << "PRODUCT H23\n";
  VectorState<T> v23 = vv.condense(MASK_BLOCK2|MASK_BLOCK3);
  VectorState<T> res23(v23);
  res23 = T(0);
  product(ss.h23, v23, res23, BLOCK2, T(1), false, &globals);
  res += res23.decondense(MASK_BLOCK2|MASK_BLOCK3, vv);

  typename std::vector<ProductTerm<T> >::iterator term_iter;
  for(term_iter = ss.product_terms.begin(); term_iter != ss.product_terms.end(); term_iter++){
    ProductTerm<T> &pterm = *term_iter;
    product_term(pterm, vv, res, mask_hc, &globals);
  }

 if(ss.project()){
    for(int i = 0; i < ss._project_states.size(); i++){
      VectorState<T> aux = *ss._project_states[i];
      T p = product(aux,vv);
      res = res + aux * p*T(100000);
    }
  }

  if(ss.verbose() > 0)
    dmtk::out() << "Elapsed time (product): " << dmtk::format_elapsed(aux_timer.elapsed()) << endl;

  return res;
}

template<class T>
void
init_terms_composite(System<T>& ss, const VectorState<T>& vv)
{
  if(!ss.rotate_terms()) return;

  VectorState<T> res(vv);
  res = T(0);
  std::vector<ProductTerm<T> > new_terms;
  if(ss.verbose() > 0) dmtk::out() << "Initializing product terms\n";
  ss.product_terms.clear();
// Interactions between blocks (1-3,2-4,2-3), and external and local terms
  VectorState<T> aux(res), aux_res(res);
  typename Hami<T>::const_iterator hiter;
//////////////////////////////////////////////////////////////////////////
// product-4+  products of 4 or more operators
//////////////////////////////////////////////////////////////////////////
  int nterm = 0;
  typename std::list<AuxTerm<T> >::const_iterator titer;
  for(titer = ss.aux_terms.begin(); titer != ss.aux_terms.end(); titer++){
    const AuxTerm<T> &auxt = (*titer);

    Term<T> t = (auxt.t);

    init_term_composite(ss, auxt, t, vv, res);
  }
//////////////////////////////////////////////////////////////////////////
// product-hami
//////////////////////////////////////////////////////////////////////////
  for(hiter = ss.h.begin(); hiter != ss.h.end(); hiter++){
    const Term<T>& t = *hiter;

    if(t.size() == 1 && t[0].is_hami()){
      const Hami<T>& h = t[0].hami();
      const BasicOp<T> *op = ss._b1->operator()(t[0]);
      if(op) {
        new_terms = get_product_terms(*op, vv, res, BLOCK1);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }
      op = ss._b2->operator()(t[0]);
      if(op) {
        new_terms = get_product_terms(*op, vv, res, BLOCK2);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }
      op = ss._b3->operator()(t[0]);
      if(op) {
        new_terms = get_product_terms(*op, vv, res, BLOCK3);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }
      op = ss._b4->operator()(t[0]);
      if(op) {
        new_terms = get_product_terms(*op, vv, res, BLOCK4);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }
    }
  }
//////////////////////////////////////////////////////////////////////////
// product-1 (only for TERM_EXTERN || TERM_LOCAL || IDENTITY)
//////////////////////////////////////////////////////////////////////////
  for(hiter = ss.h.begin(); hiter != ss.h.end(); hiter++){
    const Term<T>& t = *hiter;

    if(ss.use_coef_tol() && fabs(dmtk::real(t.coef())) <= ss.coef_tol()) continue;

    if(t.size() == 1 && t[0].name() == "I"){
      VectorState<T> aux = vv;
      aux *= t.coef();
      res += aux;
      continue;
    }

    if((t.type() == TERM_EXTERN && ss.apply_extern() && t.size() == 1 && t[0].name() != ss.h.name()) || (t.type() == TERM_LOCAL && t[0].name() != ss.h.name())){
      size_t bmask = 0;
      for(int i = 0; i < t.size(); i++)
        bmask |= mask(ss.block(t[i].site()));

      if(bmask == MASK_BLOCK1 || bmask == MASK_BLOCK2 || bmask == MASK_BLOCK3 || bmask == MASK_BLOCK4){
        BasicOp<T> top = t[0].internals();
        const BasicOp<T>* _op = ss(top);
        if(!_op && ss.verbose() >= 0) {
          dmtk::out() << "ERROR EXTERN: Operator " << top.description() << " not found\n";
        }
        if(_op){
          bool calc_hc = ss.h.use_hc();
          if (_op && _op->is_diagonal()) calc_hc = false;

          size_t pos = ss.block(t[0].site());
          if(_op){
            new_terms = get_product_terms(*_op, vv, res, pos, T(t.coef()), false);
            ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
            // obviously, op1 has to preserve the quantum numbers,
            // or we are working in the grand canonical
            if(calc_hc) {
              new_terms = get_product_terms(*_op, vv, res, pos, T(t.coef()), calc_hc);
              ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
            }
          }
        }
      }
    }
  }
}

template<class T>
VectorState<T>
product(System<T>& ss, const VectorState<T>& vv)
{
  if(!ss.rotate_terms()) {
    dmtk::out() << "WARNING: product not possible (rotate_terms == FALSE)\n";
    return vv;
  }
  return ss.hamiltonian_product(vv);
}

template<class T>
VectorState<T>
System<T>::hamiltonian_product(const VectorState<T>& vv)
{
  return product_composite(*this, vv);
}

#endif //USE_PRODUCT_DEFAULT

//////////////////////////////////////////////////////////////////////////
// product_term:
//////////////////////////////////////////////////////////////////////////
template<class T>
void
init_term_composite(System<T> &ss, const AuxTerm<T> &auxt, const Term<T> &t, const VectorState<T> &vv, VectorState<T> &res)
{
  if(!ss.rotate_terms()) return;

  std::vector<ProductTerm<T> > new_terms;
  bool calc_hc = ss.h.use_hc();
// TODO WARNING: This may cause problems ?
//#warning If you are using use_hc and you have problems, make sure you are not including the h.c. operator in the Hamiltonian
  int ndiag = 0;
  if(t.is_diagonal()) calc_hc = false;
///////////////////////////////////////////////
  if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2) || auxt.mask == (MASK_BLOCK2|MASK_BLOCK3) || auxt.mask == (MASK_BLOCK3|MASK_BLOCK4) || auxt.mask == (MASK_BLOCK1|MASK_BLOCK4) || auxt.mask == (MASK_BLOCK1|MASK_BLOCK3) || auxt.mask == (MASK_BLOCK2|MASK_BLOCK4)){

    if(ss.use_composite()){
      if(auxt.b1 < auxt.b2){
        new_terms = get_product_terms(*(auxt.ref_op[0]), auxt.sum_op, vv, res, auxt.b1, auxt.b2, T(1), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      } else {
        new_terms = get_product_terms(auxt.sum_op, *(auxt.ref_op[0]), vv, res, auxt.b2, auxt.b1, T(1), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }
    } else {
      if(auxt.b1 < auxt.b2){
        new_terms = get_product_terms(*auxt.op[size_t(auxt.b1)-1], *auxt.op[size_t(auxt.b2)-1], vv, res, auxt.b1, auxt.b2, t.coef(), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      } else {
        new_terms = get_product_terms(*auxt.op[size_t(auxt.b2)-1], *auxt.op[size_t(auxt.b1)-1], vv, res, auxt.b2, auxt.b1, t.coef(), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }
    }
  }
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){

    new_terms = get_product_terms(*(auxt.op[0]), *(auxt.op[1]), *(auxt.op[2]), *(auxt.op[3]), vv, res, BLOCK1, BLOCK2, BLOCK3, BLOCK4, T(t.coef()), calc_hc);
    ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
  }
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3)){
    if(ss.use_composite()){
      new_terms = get_product_terms(auxt.sum_op, *(auxt.ref_op[0]), *(auxt.ref_op[1]), vv, res, BLOCK1, BLOCK2, BLOCK3, T(1), calc_hc);
      ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
    } else {
      new_terms = get_product_terms(*(auxt.op[0]), *(auxt.op[1]), *(auxt.op[2]), vv, res, BLOCK1, BLOCK2, BLOCK3, T(t.coef()), calc_hc);
      ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
    }
  }
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4)){
    if(ss.use_composite()){
      if(auxt.b2 < auxt.b3){
        new_terms = get_product_terms(auxt.sum_op, *(auxt.ref_op[0]), *(auxt.ref_op[1]), vv, res, BLOCK1, BLOCK2, BLOCK4, T(1), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      } else {
        new_terms = get_product_terms(*(auxt.ref_op[1]), *(auxt.ref_op[0]), auxt.sum_op, vv, res, BLOCK1, BLOCK2, BLOCK4, T(1), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }

    } else {
      new_terms = get_product_terms(*(auxt.op[0]), *(auxt.op[1]), *(auxt.op[3]), vv, res, BLOCK1, BLOCK2, BLOCK4, T(t.coef()), calc_hc);
      ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
    }
  }
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4)){

    if(ss.use_composite()){
      if(auxt.b2 < auxt.b3){
        new_terms = get_product_terms(auxt.sum_op, *(auxt.ref_op[0]), *(auxt.ref_op[1]), vv, res, BLOCK1, BLOCK3, BLOCK4, T(1), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      } else {
        new_terms = get_product_terms(*(auxt.ref_op[1]), *(auxt.ref_op[0]), auxt.sum_op, vv, res, BLOCK1, BLOCK3, BLOCK4, T(1), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }
    } else {
      new_terms = get_product_terms(*(auxt.op[0]), *(auxt.op[2]), *(auxt.op[3]), vv, res, BLOCK1, BLOCK3, BLOCK4, T(t.coef()), calc_hc);
      ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
    }
  }
  else if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
    if(ss.use_composite()) {
      new_terms = get_product_terms(*(auxt.ref_op[0]), *(auxt.ref_op[1]), auxt.sum_op, vv, res, BLOCK2, BLOCK3, BLOCK4, T(1), calc_hc);
      ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
    } else {
      new_terms = get_product_terms(*(auxt.op[1]), *(auxt.op[2]), *(auxt.op[3]), vv, res, BLOCK2, BLOCK3, BLOCK4, T(t.coef()), calc_hc);
      ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
    }
  }
}


// locate() keeps the block offsets that block() and operator()(op) each
// computed in two heap-allocated vectors per call; the first block whose
// sites contain `site` wins, as in the per-site scan it replaces.
template <class T>
int
System<T>::locate(int site, int &local_site) const
{
  std::array<const Block<T>*, 4> const b = {_b1, _b2, _b3, _b4};
  std::array<int, 4> offset;
  offset[0] = _in_warmup && _grow_symmetric && _grow_outward && size() < h.lattice().size() ? h.lattice().size()/2 - _b1->lattice().size() - _b2->lattice().size() : 0;
  offset[1] = _in_warmup && _grow_symmetric && _grow_outward && size() < h.lattice().size() ? offset[0] + _b1->lattice().size() : _b1->lattice().size();
  offset[2] = offset[1] + _b2->lattice().size();
  offset[3] = offset[2] + _b3->lattice().size();

  if(_in_warmup && _grow_symmetric && size() < h.lattice().size() && !_grow_outward){
    offset[3] = h.lattice().size() - _b4->lattice().size();
    offset[2] = offset[3] - _b3->lattice().size();
  }

  for(int i = 0; i < 4; i++){
    int const sites = int(b[i]->lattice().size());
    if (site >= offset[i] && site < offset[i] + sites) {
      local_site = b[i]->single_site() ? site - offset[i] : site;
      return i;
    }
  }
  return -1;
}

template <class T>
size_t
System<T>::block(int site) const
{
  int local_site;
  int const i = locate(site, local_site);
  return i < 0 ? BLOCK_NONE : size_t(i + 1);
}

template<class T>
const BasicOp<T>*
System<T>::operator()(const BasicOp<T>& op) const
{
  int local_site;
  int const i = locate(op.site(), local_site);
  if (i < 0) return 0;
  std::array<const Block<T>*, 4> const b = {_b1, _b2, _b3, _b4};
  BasicOp<T> _op(op);
  _op.set_site(local_site);
  return (b[i]->operator()(_op));
}

template<class T>
const BasicOp<T>*
System<T>::indexed_op(const BasicOp<T>& op, OpIndex &index) const
{
  int local_site;
  int const i = locate(op.site(), local_site);
  if (i < 0) return 0;
  if (!index.built[i]) {
    std::array<const Block<T>*, 4> const b = {_b1, _b2, _b3, _b4};
    for (BasicOp<T> const & o : *b[i])
      index.ops[i].emplace(OpKey(o.name(), o.site(), o.internal_site()), &o);
    index.built[i] = true;
  }
  auto const found = index.ops[i].find(OpKey(op.name(), local_site, op.internal_site()));
  return found == index.ops[i].end() ? 0 : found->second;
}

//////////////////////////////////////////////////////////////////

} // namespace dmtk

#endif // DMTK_SYSTEM_H_
