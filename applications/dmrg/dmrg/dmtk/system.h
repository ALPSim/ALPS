/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
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

#ifndef __DMTK_SYSTEM_H__
#define __DMTK_SYSTEM_H__

//*******************************************************************
// System Class: Defines system properties and performs DMRG loops
// This version of the system class works only on symmetric systems
//*******************************************************************
 
#include <iostream>
#include <iosfwd>
#include <list>
#include <string>
#include <iomanip>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "enums.h"
#include "vector.h"
#include "matrix.h"
#include "state.h"
#include "operators.h"
#include "block.h"
#include "lattice.h"
#include "hami.h"
#include "ctimer.h"
#include "util.h"
#include "globals.h"
#ifdef WITH_PTHREADS
#include <pthread.h>
#include <unistd.h>
#endif // WITH_PTHREADS

#include <sys/stat.h>

#ifdef BOOST_MSVC
#include <io.h>
#endif

#include <boost/filesystem/path.hpp>
#include <alps/utility/temporary_filename.hpp>

using namespace std;

namespace dmtk
{

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
   cout << "PRODUCT " << term.description() << endl;

#define COUT_PRODUCT(term) \
 if(calc_hc) \
   cout << "PRODUCT " << term.description() << " + h.c." << endl; \
 else  \
   cout << "PRODUCT " << term.description() << endl;

template<class T> 
VectorState<T> product(const System<T>& ss, const VectorState<T>& vv);

template<class T> 
VectorState<T> product_default(const System<T>& ss, const VectorState<T>& vv, const Hami<T>* hami = NULL, bool only_local = false);

class FileList 
{
  private:
    std::map<std::string,std::string> _tmp_filenames;
    std::string last_filename;
    boost::filesystem::path temp_dir;
  public:
    FileList() { set_temp_dir("."); };
    FileList(const char *dir) { set_temp_dir(dir); }

    void set_temp_dir(const char *dir) 
      { 
         std::string aux = std::string(dir);
#ifndef BOOST_MSVC
         struct stat dir_ptr;
         if((!stat(aux.c_str(),&dir_ptr)) == 0) {
             std::cerr << "*** ERROR: ALPS DMRG could not open directory for temporary files. Create the directory " + aux + " or choose a different path.\n";
                         boost::throw_exception(std::runtime_error("*** ERROR: ALPS DMRG could not open directory for temporary files. Create the directory " + aux + " or choose a different path."));
         }
#endif
         std::cout << "ALPS DMRG temporary files will be written to " << aux << std::endl;
         temp_dir = boost::filesystem::path(aux); 
      }

    const char * get_filename(const char *input) 
      {
         std::string filename(input);
         std::map<std::string,std::string>::iterator old_name = _tmp_filenames.find(filename);
         if(old_name != _tmp_filenames.end()) return (old_name->second).c_str();

         filename = filename.substr(0,filename.find_first_of('.'));
         this->last_filename = alps::temporary_filename((temp_dir / filename).string());
         _tmp_filenames[std::string(input)] = this->last_filename;
         std::cout << "Creating temp file " << this->last_filename << std::endl;
         return this->last_filename.c_str();
      }
};

FileList tmp_files;

template<class T>
class System
{
  protected:
    typedef typename dmtk::Block<T> B;

    Lattice _lattice;

    CTimer _timer;

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
    bool _grow_symmetric;
    bool _grow_outward;
    bool _custom_qns;
    bool _store_products;
    bool _measure_symmetric;
    bool _apply_hami;
    bool _apply_extern;
    bool _full_sweep; // run sweep from 1 to ls-3/ls-4
    bool _use_coef_tol; 
    double _coef_tol;
    Matrix<size_t> _sweeps;
    size_t _numsweeps;

    char _name[255];

    double _lanczos_tol;
    int _lanczos_maxiter;

    virtual void init_iteration(const B&b1, const B&b2, const B&b3, const B&b4, bool use_seed = false);
    void rotate_hami(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *this_hami = NULL);
    void rotate_terms(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *this_hami = NULL);
    void rotate_corr(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *this_hami = NULL);

    bool signal_emit(size_t signal_id, void * data)
      {
        if(signal_handler)
          return signal_handler(*this, signal_id, data);
        return true;
      }
    void * _data;

  public:

    bool (*signal_handler) (System<T> &, size_t signal_id, void *data);

    B rightblock;
    B leftblock;
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
    , _grow_symmetric(true)
    , _grow_outward(false)
    , _custom_qns(true)
    , _store_products(true)
    , _measure_symmetric(false)
    , _apply_hami(true)
    , _apply_extern(true)
    , _full_sweep(false)
    , _use_coef_tol(false)
    , _coef_tol(1.e-5) 
    , _numsweeps(0)
    , _lanczos_tol(1.e-7)
    , _lanczos_maxiter(-1)
    , qnt(0)
    , _ntargets(1)
    , _nstates(1)
    { _sweeps.resize(2,10); _target.resize(1); _target_weight.resize(1); set_name(""); };

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
     , _grow_symmetric(true)
     , _grow_outward(false)
     , _custom_qns(true)
     , _store_products(true)
     , _measure_symmetric(false)
     , _apply_hami(true)
     , _apply_extern(true)
     , _full_sweep(false)
     , _use_coef_tol(false)
     , _coef_tol(1.e-5) 
     , _numsweeps(0)
     , _lanczos_tol(1.e-7)
     , _lanczos_maxiter(-1)
     , qnt(0)
     , _ntargets(1)
     , _nstates(1)
     , h(_h)
      { set_name(the_name); m1 = m2 = m3 = m4 = _h.get_site(0).dim(); m = m1 * m2; _sweeps.resize(2,10); _target.resize(1); _target_weight.resize(1); _target_weight[0] = double(1); }

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
          cout << "NEW SWEEP " << _sweep << " " << _sweeps(0,_sweep-1) << " " << _sweeps(1,_sweep-1) << endl;
          sweep(_sweeps(0,_sweep-1),_sweeps(1,_sweep-1),(size_t)resume_dir,resume_iter);
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
      { sprintf(_name, "%s", the_name); }
    const char* name() const { return _name; }

    virtual void diagonalize(bool use_seed = true);
    virtual void build_target_states(int pos = 0);
    virtual void truncate(int block, int new_size);
    virtual void rotate(int block, Block<T>& b); 
    virtual void warmup_loop(size_t t);
    virtual void warmup_loop(size_t t, const Vector<QN> &qns);
    virtual void sweep(size_t t1, size_t t2, size_t dir=RIGHT2LEFT, int _start=1);
    virtual void main_loop(size_t nsweeps, size_t t);
    virtual void init_correlations();
    virtual void final_sweep(size_t t, size_t dir, int _start, bool _rotate);
    virtual void measure(const VectorState<T> *v=NULL); 
    virtual void measure_n(size_t n, const VectorState<T> *v=NULL); 
    virtual T measure_operator(const BasicOp<T> &top, const VectorState<T> *v=NULL);

    size_t get_iter() const { return iter; }
    size_t get_dir() const { return dir; }


    size_t block(int site) const;

    const BasicOp<T>* operator()(const BasicOp<T>& op) const;
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

/*
    size_t new_site(size_t b, size_t site) const
      {
        size_t nsite;

        switch(b){
          case BLOCK1: 
          case BLOCK3:
            nsite = site; break;
          case BLOCK2: 
            nsite = _b1->lattice().size() + site; break;
          case BLOCK4: 
            nsite = _b3->lattice().size() + site; break;
        }
        if(verbose() > 0) cout << "NEW SITE " << site << " => " << nsite << endl;
        return nsite;
      }
*/

    size_t size() const 
      { return _b1->lattice().size()+ _b2->lattice().size()+
               _b3->lattice().size()+ _b4->lattice().size(); }
    int site1() const 
      { if(dir == LEFT2RIGHT) return iter; else return h.lattice().size()-iter-2; }
    int site2() const { return site1()+1; }

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

    CTimer timer() const { return _timer; }
    size_t sweep_direction() const { return dir; }
    int iteration() const { return iter; }

    void set_lanczos_tolerance(double tol) { _lanczos_tol = tol; }
    double lanczos_tolerance() { return _lanczos_tol; }
    void set_lanczos_maxiter(int n) { _lanczos_maxiter = n; }
    int lanczos_maxiter() { return _lanczos_maxiter; }
    double truncation_error() const { return _truncation_error; }
    double entropy() const { return _entropy; }
    System<T>& set_data (void *data) { _data = data; return *this; }
    void *get_data () const { return _data;}

    // Streams

    void write_status() const
    {
      const char *file = tmp_files.get_filename("system.dat"); 
      ofstream outputfile(file,std::ios::out|std::ios::binary);
      if(!outputfile) {
        cerr << "*** ERROR: Could not open file " << file << endl;
        return;
      }
      outputfile.write((const char *)&_numsweeps, sizeof(size_t));
      for(int i = 1; i <= _numsweeps; i++){
        int l = _sweeps(0,i);
        outputfile.write((const char *)&l, sizeof(int));
        l = _sweeps(1,i);
        outputfile.write((const char *)&l, sizeof(int));
      }
      outputfile.write((const char *)&_in_warmup, sizeof(bool));
      outputfile.write((const char *)&_sweep, sizeof(int));
      outputfile.write((const char *)&dir, sizeof(size_t));
      outputfile.write((const char *)&iter, sizeof(int));
      outputfile.write((const char *)&_ntargets, sizeof(int));
/*
      outputfile.write((const char *)&tmp_filenames.size(), sizeof(int));
      std::map<std::string,std::string>::const_iterator miter;
      for(miter = _tmp_filenams.begin(); miter != _tmp_filenames.end(); miter++){
        outputfile.write((const char *)*miter; (*miter).size());
      } 
*/
      qnt.write(outputfile);
      qn.write(outputfile);

      outputfile.close();
    }

    void read_status() 
    {
      const char *file = tmp_files.get_filename("system.dat"); 
      ifstream inputfile(file,std::ios::in|std::ios::binary);
      if(!inputfile) {
        cerr << "*** ERROR: Could not open file " << file << endl;
        return;
      }
      inputfile.read((char *)&_numsweeps, sizeof(size_t));
      _sweeps.resize(2,_numsweeps+1);
      for(int i = 1; i <= _numsweeps; i++){
        int l;
        inputfile.read((char *)&l, sizeof(int));
        _sweeps(0,i) = l;
        inputfile.read((char *)&l, sizeof(int));
        _sweeps(1,i) = l;
      }
      inputfile.read((char *)&_in_warmup, sizeof(bool));
      inputfile.read((char *)&_sweep, sizeof(int));
      inputfile.read((char *)&dir, sizeof(size_t));
      inputfile.read((char *)&iter, sizeof(int));
      inputfile.read((char *)&_ntargets, sizeof(int));
      inputfile.read((char *)&_nstates, sizeof(int));
      qnt.read(inputfile);
      qn.read(inputfile);

      inputfile.close();
    }

    void write(std::ostream &s) const
    {
      s.write((const char *)&_numsweeps, sizeof(size_t));
      for(int i = 1; i <= _numsweeps; i++){
        int l = _sweeps(0,i);
        s.write((const char *)&l, sizeof(int));
        l = _sweeps(1,i);
        s.write((const char *)&l, sizeof(int));
      }
      s.write((const char *)&_in_warmup, sizeof(bool));
      s.write((const char *)&_sweep, sizeof(int));
      s.write((const char *)&dir, sizeof(size_t));
      s.write((const char *)&iter, sizeof(int));
      s.write((const char *)&_ntargets, sizeof(int));
      s.write((const char *)&_nstates, sizeof(int));
      s.write((const char *)&_propagate_state.size(), sizeof(int));
      qnt.write(s);
      qn.write(s);
      gs.write(s);
      for(int i = 0; i < _ntargets; i++){
        s.write((const char *)&_target_weight[i], sizeof(double));
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
      s.read((char *)&_numsweeps, sizeof(size_t));
      _sweeps.resize(2,_numsweeps+1);
      for(int i = 1; i <= _numsweeps; i++){
        int l;
        s.read((char *)&l, sizeof(int));
        _sweeps(0,i) = l;
        s.read((char *)&l, sizeof(int));
        _sweeps(1,i) = l;
      }
      s.read((char *)&_in_warmup, sizeof(bool));
      s.read((char *)&_sweep, sizeof(int));
      s.read((char *)&dir, sizeof(size_t));
      s.read((char *)&iter, sizeof(int));
      s.read((char *)&_ntargets, sizeof(int));
      int nprop;
      s.read((char *)&nprop, sizeof(int));
      _propagate_state.resize(nprop);
      qnt.read(s);
      qn.read(s);
      gs.read(s);
      _target.resize(_ntargets);
      _target_weight.resize(_ntargets);
      for(int i = 0; i < _ntargets; i++){
        s.read((char *)&_target_weight[i], sizeof(double));
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
        sprintf(file,"gs_%s_%i_l.dat",_name,n);
      else
        sprintf(file,"gs_%s_%i_r.dat",_name,n);

      ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::binary);
      if(!outputfile) 
        cerr << "*** ERROR: Could not open file " << file << endl;
      else{
#ifdef DMTK_DEBUG
        cout << "WRITING GS " << file << endl;
#endif // DMTK_DEBUG

        gs.write(outputfile);
  
        int _n = _target.size();
        outputfile.write((const char *)&_n, sizeof(int));
        for(int i = 0; i < _n; i++){
          _target[i].write(outputfile);
        }
  
        _n = _propagate_state.size();
        outputfile.write((const char *)&_n, sizeof(int));
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
        sprintf(file,"gs_%s_%i_l.dat",_name,n);
      else
        sprintf(file,"gs_%s_%i_r.dat",_name,n);

      ifstream inputfile(tmp_files.get_filename(file),std::ios::in|std::ios::binary);
      if(!inputfile) 
        cerr << "*** ERROR: Could not open file " << file << endl;
      else{
#ifdef DMTK_DEBUG
        cout << "READING GS " << file << endl;
#endif // DMTK_DEBUG

        gs.read(inputfile);
  
        inputfile.read((char *)&_ntargets, sizeof(int));
        _target.resize(_ntargets);
        _target_weight.resize(_ntargets);
        for(int i = 0; i < _ntargets; i++){
          _target[i].read(inputfile);
        }
  
        int nprop;
        inputfile.read((char *)&nprop, sizeof(int));
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
        sprintf(file,"rho_%s_%i_l.dat",_name,n);
      else
        sprintf(file,"rho_%s_%i_r.dat",_name,n);

#ifdef DMTK_DEBUG
      cout << "SAVING RHO " << file << endl;
#endif // DMTK_DEBUG
      ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::binary);
      if(!outputfile) 
        cerr << "*** ERROR: Could not open file " << file << endl;
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
        sprintf(file,"rho_%s_%i_l.dat",_name,n);
      else
        sprintf(file,"rho_%s_%i_r.dat",_name,n);

#ifdef DMTK_DEBUG
      cout << "READING RHO " << file << endl;
#endif // DMTK_DEBUG
      ifstream inputfile(tmp_files.get_filename(file),std::ios::in|std::ios::binary);
      if(!inputfile) 
        cerr << "*** ERROR: Could not open file " << file << endl;
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
        sprintf(file,"block_%s_%i_l.dat",_name,n);
      else
        sprintf(file,"block_%s_%i_r.dat",_name,n);

#ifdef DMTK_DEBUG
      cout << "SAVING BLOCK " << file << endl;
#endif // DMTK_DEBUG
      ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::binary);
      if(!outputfile) 
        cerr << "*** ERROR: Could not open file " << file << endl;
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
        sprintf(file,"block_%s_%i_l.dat",_name,n);
      else
        sprintf(file,"block_%s_%i_r.dat",_name,n);

#ifdef DMTK_DEBUG
      cout << "READING BLOCK " << file << endl;
#endif // DMTK_DEBUG
      ifstream inputfile(tmp_files.get_filename(file),std::ios::in|std::ios::binary);
      if(!inputfile) 
        cerr << "*** ERROR: Could not open file " << file << endl;
      else{
        b.read(inputfile);
        inputfile.close();
      }
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
  _timer.Start();
  bool save_use_seed = _use_seed;
  _use_seed = false;

  char file[255];
  sprintf(file,"iter_%s.dat",_name);
  ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::app);
  if(!outputfile) cerr << "*** ERROR: could not open " << file << endl;

  start();
//****************************************************************

// LEFT-TO-RIGHT sweep to get B(1)B(2)...B(L/2-1)
  cout << "LEFT-TO-RIGHT sweep to get B(1)B(2)...B(L/2-1)\n";
  outputfile << "LEFT-TO-RIGHT sweep to get B(1)B(2)...B(L/2-1)\n";

  dir = LEFT2RIGHT;
  int sweep_max = full_sweep() ? ls-2 : ls-3;
//  int sweep_max = ls-2;
  int mid_size = (ls/2)*2 == ls ? ls/2 : ls/2+1;
  int max = _grow_symmetric ? mid_size : sweep_max;

  for(iter = 1; iter < max; iter++)
  {
    signal_emit(SYSTEM_SIGNAL_START_ITER, NULL);
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

    _timer.Lap();
    cout << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;
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
      cout << "===========================================\n";
      cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
      cout << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
      outputfile << "===========================================\n";
    }
    signal_emit(SYSTEM_SIGNAL_END_ITER, NULL);
  }

  qn = qnt;
  if(_grow_symmetric){

//  LEFT-TO-RIGHT sweep to get B(L/2)B(L/2+1)...B(L-3)
    cout << "LEFT-TO-RIGHT sweep to get B(L/2)B(L/2+1)...B(L-3)\n";
    outputfile << "LEFT-TO-RIGHT sweep to get B(L/2)B(L/2+1)...B(L-3)\n";

    for(iter = mid_size; iter < sweep_max; iter++)
    {
      signal_emit(SYSTEM_SIGNAL_START_ITER, NULL);
      read_block(rightblock, ls-iter-2, RIGHT);
      read_block(leftblock, iter, LEFT);
      const Block<T>& site1 = h.get_site(iter);
      const Block<T>& site2 = h.get_site(iter+1);
  
      _timer.Lap();
      cout << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;
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
        cout << "===========================================\n";
        cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
        cout << "===========================================\n";
        outputfile << "===========================================\n";
        outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
        outputfile << "===========================================\n";
      }
      signal_emit(SYSTEM_SIGNAL_END_ITER, NULL);
    }
  } else {
    cout << "RIGHT-TO-LEFT sweep to get B(1)B(2)...B(L/2-1)\n";
    outputfile << "RIGHT-TO-LEFT sweep to get B(1)B(2)...B(L/2-1)\n";
    dir = RIGHT2LEFT;

    for(iter = 1; iter < sweep_max; iter++)
    {
      signal_emit(SYSTEM_SIGNAL_START_ITER, NULL);
      read_block(leftblock, ls-iter-2, LEFT);
      read_block(rightblock, iter, RIGHT);
      const Block<T>& site1 = h.get_site(ls-2-iter);
      const Block<T>& site2 = h.get_site(ls-2-iter+1);

      m1 = leftblock.dim();
      m4 = rightblock.dim();
      m2 = site1.dim();
      m3 = site2.dim();

      _timer.Lap();
      cout << "WARMUP RIGHT-TO-LEFT ITERATION " << iter << endl;
      outputfile << "WARMUP RIGHT-TO-LEFT ITERATION " << iter << endl;

      init_iteration(leftblock, site1, site2, rightblock, false); 
      diagonalize(false); 

      m4 = std::min(m4*m3,m);

      truncate(RIGHT, m4);
      rotate(RIGHT, newblock);
      write_iter(RIGHT);

      if(verbose() > 0) {
        cout << "===========================================\n";
        cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
        cout << "===========================================\n";
        outputfile << "===========================================\n";
        outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
        outputfile << "===========================================\n";
      }
      signal_emit(SYSTEM_SIGNAL_END_ITER, NULL);
    }

    dir = LEFT2RIGHT;
    _use_seed = save_use_seed;

    for(iter = 1; iter < sweep_max; iter++)
    {
      signal_emit(SYSTEM_SIGNAL_START_ITER, NULL);
      read_block(rightblock, ls-iter-2, RIGHT);
      read_block(leftblock, iter, LEFT);
      const Block<T>& site1 = h.get_site(iter);
      const Block<T>& site2 = h.get_site(iter+1);
  
      m1 = leftblock.dim();
      m4 = rightblock.dim();
      m2 = site1.dim();
      m3 = site2.dim();
 
      _timer.Lap();
      cout << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;
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
        cout << "===========================================\n";
        cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
        cout << "===========================================\n";
        outputfile << "===========================================\n";
        outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
        outputfile << "===========================================\n";
      }
      signal_emit(SYSTEM_SIGNAL_END_ITER, NULL);
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
  sprintf(file,"iter_%s.dat",_name);
  ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::app);
  if(!outputfile) cerr << "*** ERROR: could not open " << file << endl;

  cout << "===========================================\n";
  cout << "NEW SWEEP\n";
  cout << "===========================================\n";
  outputfile << "===========================================\n";
  outputfile << "NEW SWEEP\n";
  outputfile << "===========================================\n";

// RIGHT-TO-LEFT sweep
  
  dir = _dir;

  if(dir == RIGHT2LEFT){
    for(iter = start; iter < sweep_max; iter++)
    {
      signal_emit(SYSTEM_SIGNAL_START_ITER, NULL);
      read_block(leftblock, ls-iter-2, LEFT);
      read_block(rightblock, iter, RIGHT);
      const Block<T>& site1 = h.get_site(ls-2-iter);
      const Block<T>& site2 = h.get_site(ls-2-iter+1);

      m1 = leftblock.dim();
      m4 = rightblock.dim();
      m2 = site1.dim();
      m3 = site2.dim();

      _timer.Lap();
      cout << "RIGHT-TO-LEFT ITERATION " << iter << endl;
      outputfile << "RIGHT-TO-LEFT ITERATION " << iter << endl;

      init_iteration(leftblock, site1, site2, rightblock, _use_seed); 
      diagonalize(_use_seed); 

      m = t1;
      m4 = std::min(m4*m3,m);

      truncate(RIGHT, m4);
      rotate(RIGHT, newblock);
      write_iter(RIGHT);

      if(verbose() > 0) {
        cout << "===========================================\n";
        cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
        cout << "===========================================\n";
        outputfile << "===========================================\n";
        outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
        outputfile << "===========================================\n";
      }
      signal_emit(SYSTEM_SIGNAL_END_ITER, NULL);
    }

    start = 1;
  }

// LEFT-TO-RIGHT sweep

  dir = LEFT2RIGHT;

  for(iter = start; iter < sweep_max; iter++)
  {
    signal_emit(SYSTEM_SIGNAL_START_ITER, NULL);
    read_block(rightblock, ls-iter-2, RIGHT);
    read_block(leftblock, iter, LEFT);
    const Block<T>& site1 = h.get_site(iter);
    const Block<T>& site2 = h.get_site(iter+1);

    m1 = leftblock.dim();
    m4 = rightblock.dim();
    m2 = site1.dim();
    m3 = site2.dim();

    _timer.Lap();
    cout << "LEFT-TO-RIGHT ITERATION " << iter << endl;
    outputfile << "LEFT-TO-RIGHT ITERATION " << iter << endl;

    init_iteration(leftblock, site1, site2, rightblock, _use_seed); 
    diagonalize(_use_seed);

    m = t2;
    m1 = std::min(m1*m2,m);

    truncate(LEFT, m1);
    rotate(LEFT, newblock);
    write_iter(LEFT);

    if(verbose() > 0) {
      cout << "===========================================\n";
      cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
      cout << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
      outputfile << "===========================================\n";
    }
    signal_emit(SYSTEM_SIGNAL_END_ITER, NULL);
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

  char file[255];
  sprintf(file,"iter_%s.dat",_name);
  ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::app);
  if(!outputfile) cerr << "*** ERROR: could not open " << file << endl;
  outputfile << "Last iteration to get a symmetric block B(L/2-1)..B(L/2-1)\n";

//****************************************************************
// We store in site blocks the operators used for measuring correlations
//****************************************************************
  init_correlations();
//****************************************************************

  dir = _dir;
  int sweep_max = full_sweep() ? ls-2 : ls-3;

  if(dir == RIGHT2LEFT){
    for(iter = _start; iter < sweep_max; iter++) // ls/2-1; iter++);
    {
      signal_emit(SYSTEM_SIGNAL_START_ITER, NULL);
      read_block(rightblock, iter, RIGHT);
      read_block(leftblock, ls-iter-2, LEFT);
      const Block<T>& site1 = h.get_site(ls-2-iter);
      const Block<T>& site2 = h.get_site(ls-2-iter+1);

      m1 = leftblock.dim();
      m4 = rightblock.dim();
      m2 = site1.dim();
      m3 = site2.dim();
   
      _timer.Lap();
      cout << "FINAL SWEEP ITERATION " << endl;
      outputfile << "FINAL SWEEP ITERATION " << endl;
      cout << "RIGHT-TO-LEFT ITERATION " << iter << endl;
      outputfile << "RIGHT-TO-LEFT ITERATION " << iter << endl;
  
      init_iteration(leftblock, site1, site2, rightblock, _use_seed);
  
      diagonalize(_use_seed);
      if(!_store_products) measure();
  
      m = t;
      m4 = std::min(m4*m3,m);
  
      truncate(RIGHT, m4);
      rotate(RIGHT, newblock);
      write_iter(RIGHT);

      if(verbose() > 0) {
        cout << "===========================================\n";
        cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
        cout << "===========================================\n";
        outputfile << "===========================================\n";
        outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
        outputfile << "===========================================\n";
      }
      signal_emit(SYSTEM_SIGNAL_END_ITER, NULL);
    }
    _start = 1;
  }

  dir = LEFT2RIGHT;
  sweep_max = full_sweep() ? ls-2 : ls/2-1;

  for(iter = _start; iter < sweep_max; iter++)
  {
    signal_emit(SYSTEM_SIGNAL_START_ITER, NULL);
    read_block(rightblock, ls-iter-2, RIGHT);
    read_block(leftblock, iter, LEFT);
    const Block<T>& site1 = h.get_site(iter);
    const Block<T>& site2 = h.get_site(iter+1);

    m1 = leftblock.dim();
    m4 = rightblock.dim();
    m2 = site1.dim();
    m3 = site2.dim();

    _timer.Lap();
    cout << "FINAL SWEEP ITERATION " << endl;
    outputfile << "FINAL SWEEP ITERATION " << endl;
    cout << "LEFT-TO-RIGHT ITERATION " << iter << endl;
    outputfile << "LEFT-TO-RIGHT ITERATION " << iter << endl;

    init_iteration(leftblock, site1, site2, rightblock, _use_seed); 

    diagonalize(_use_seed);
    if(!_store_products) measure();

    m = t;
    m1 = std::min(m1*m2,m);

    truncate(LEFT, m1);
    rotate(LEFT, newblock);
    write_iter(LEFT);

    if(verbose() > 0) {
      cout << "===========================================\n";
      cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
      cout << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
      outputfile << "===========================================\n";
    }
  }

// Calculate the ground state in the symmetric system

  if(!full_sweep()){
    cout << "FINAL SWEEP - LAST ITERATION " << endl;
    outputfile << "FINAL SWEEP - LAST ITERATION " << endl;
    read_block(rightblock, IS_EVEN(ls) ? ls/2-1 : ls/2, RIGHT);
    read_block(leftblock, ls/2-1, LEFT);
    const Block<T>& site1 = h.get_site(ls/2-1);
    const Block<T>& site2 = h.get_site(ls/2);
    init_iteration(leftblock, site1, site2, rightblock, _use_seed);  
//  set_lanczos_tolerance(std::min(1.e-16,_lanczos_tol));
//    set_lanczos_maxiter(30);
    diagonalize(_use_seed);
    if(!_store_products) measure();
    if(_rotate){
      m = t;
      m1 = std::min(m1*m2,m);
                                                                               
      truncate(LEFT, m1);
      rotate(LEFT, newblock);
      write_iter(LEFT);
    } else {
      write_gs(gs, ls/2-1, LEFT);
    }
  }
  signal_emit(SYSTEM_SIGNAL_END_ITER, NULL);
  
  outputfile.close();
}

/////////////////////////////////////////////////////////////////////////
// main_loop: 
// if we don't care about increasing the number of states 
// in the block we can just run main_loop with a given number of states.
/////////////////////////////////////////////////////////////////////////

template<class T>
void
System<T>::main_loop(size_t nsweeps, size_t t)
{
  m = t;
  for(int _sweep = 1; _sweep <= nsweeps; _sweep++) { sweep(m,m); }
}

/////////////////////////////////////////////////////////////////////////
// init_iteration: 
// initialize the blocks, quantum numbers, and transform the ground
// state of the previous iteration as a new seed for the new one.
/////////////////////////////////////////////////////////////////////////
template <class T>
void
System<T>::init_iteration(const B&b1, const B&b2, const B&b3, const B&b4, bool use_seed)
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
            double nd = (lnow * ((double)qnt[j].get_twice()*0.5 / (double)ltotal));
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
       cout << "TOTAL SIZE : " << size() << " " << lnow << endl;
       cout << "LEFT BLOCK SIZE : " << _b1->lattice().size() << " " << _b1->n_orbitals() << endl;
       cout << "RIGHT BLOCK SIZE : " << _b4->lattice().size() << " " << _b4->n_orbitals() << endl;
       cout << "|" << _b1->n_orbitals() << "|-|" << _b2->n_orbitals() << "|-|" << _b3->n_orbitals() << "|-|" << _b4->n_orbitals() << "|" << endl; 
       cout << "|" << _b1->lattice().size() << "|-|" << _b2->lattice().size() << "|-|" << _b3->lattice().size() << "|-|" << _b4->lattice().size() << "|" << endl; 
       for(int i = 0; i < QN::max_index(); i++)
         if(_grand_canonical & (1 << i)) cout << QN::qn_name(i) << " = " << qnt[i] << endl;
    }
  char file[255];
  sprintf(file,"iter_%s.dat",_name);
  ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::app);
  if(!outputfile) cerr << "*** ERROR: could not open " << file << endl;
  outputfile << "TOTAL SIZE : " << size() << " " << lnow << endl;
  outputfile << "LEFT BLOCK SIZE : " << _b1->lattice().size() << " " << _b1->n_orbitals() << endl;
  outputfile << "RIGHT BLOCK SIZE : " << _b4->lattice().size() << " " << _b4->n_orbitals() << endl;
  outputfile << "|" << _b1->n_orbitals() << "|-|" << _b2->n_orbitals() << "|-|" << _b3->n_orbitals() << "|-|" << _b4->n_orbitals() << "|" << endl; 
  outputfile << "|" << _b1->lattice().size() << "|-|" << _b2->lattice().size() << "|-|" << _b3->lattice().size() << "|-|" << _b4->lattice().size() << "|" << endl; 

  if(verbose() > 0) 
    for(int i = 0; i < QN::max_index(); i++)
      if(_grand_canonical & (1 << i)) cout << QN::qn_name(i) << " = " << qn[i] << endl;

  if(use_seed){
     BMatrix<T> rho1;
     Basis basis1;
     BMatrix<T> rho2;
     Basis basis2;

     CTimer clock;
     clock.Start(); 

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

           if(verbose() > 0) cout << "NEW SEED " << gs.size() << endl; 
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
             cout << "Lap: " << clock.LapTime().c_str() << endl;
         } else {
           read_rho(rho1, basis1, lattice().ls()-3, RIGHT);
           read_rho(rho2, basis2, 2, LEFT);
           read_gs(gs, lattice().ls()-4,RIGHT);
           gs.resize(_grand_canonical);

           if(verbose() > 0) cout << "NEW SEED " << gs.size() << endl; 
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
             cout << "Lap: " << clock.LapTime().c_str() << endl;
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

         if(verbose() > 0) cout << "NEW SEED " << gs.size() << endl; 
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
           cout << "Lap: " << clock.LapTime().c_str() << endl;
       } else {
         read_rho(rho1, basis1, iter, RIGHT);
         read_rho(rho2, basis2, lattice().ls()-iter-2+1, LEFT);
         read_gs(gs, iter-1, RIGHT);
         gs.resize(_grand_canonical);

         if(verbose() > 0) cout << "NEW SEED " << gs.size() << endl; 
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
           cout << "Lap: " << clock.LapTime().c_str() << endl;
       }

     }
 
  }

  gs.set_qn_mask(qn, _grand_canonical);
  gs.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);

}

#include "lanczos.cc"

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

  CTimer clock;
  clock.Start(); 

  if(verbose() > 0)
    cout << "Initialization time: " << clock.LapTime().c_str() << endl;

  clock.Start();

  _state.resize(_nstates);
  for(int i = 0; i < _nstates; i++){
     _state[i].set_qn_mask(qn, _grand_canonical);
     _state[i].resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
  }

//  if(verbose() > 0) 
    cout << "NUMBER OF STATES: " << _b1->dim() << " " << _b2->dim() << " " << _b3->dim() << " " <<  _b4->dim() << " " << gs.size() << endl;

  int n = _lanczos_maxiter;
  if(_target_subspaces && _in_warmup && lattice().size() != size()) n = 0;

  char file[255];
  sprintf(file,"vectors_%s.dat",_name);

  if(verbose() > 0) cout << ">>>>>>> LANCZOS <<<<<<<" << endl;

  if(real(product(seed,seed)) < 0.1) use_seed = false;

  if(gs.size() > 1) {
//    verif_hamiltonian<T, System<T>, VectorState<T> >(*this, _state[0]);
    lanczos<T, System<T>, VectorState<T> >(*this, _state, seed, e, _nstates, a, b, n, tol, use_seed, true, tmp_files.get_filename(file));
    gs = _state[0];
  } else {
    gs = T(1);
    VectorState<T> aux = gs;
    aux = product(*this,gs);
    T xe = product(aux,gs);
    e[0] = real(xe);
    _state[0] = gs; 
    cout << "E0 = " << e[0] << endl;
    cout << "ENER = " << e[0] << endl;
  }

  energy.resize(_calc_gap+1);
  energy[0] = e[0];

  if(verbose() && use_seed){
    cout << "OVERLAP " << product(seed,gs) << endl;
  } 

  sprintf(file,"iter_%s.dat",_name);
  ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::app);
  if(!outputfile) cerr << "*** ERROR: could not open " << file << endl;

  outputfile << _b1->lattice().size() << " " << _b2->lattice().size() << " " << _b3->lattice().size() << " " <<  _b4->lattice().size() << endl;
  outputfile << _b1->dim() << " " << _b2->dim() << " " << _b3->dim() << " " <<  _b4->dim() << " " << gs.size() << endl;

  outputfile.precision(10);
  outputfile << "ENER = " << e[0] << endl;
  if(verbose() > 0){
    outputfile << "Lanczos time: " << clock.LapTime().c_str() << endl;
    cout << "-------------------------------------------\n";
    cout << "Lanczos time: " << clock.LapTime().c_str() << endl;
    cout << "-------------------------------------------\n";
  }
  outputfile.close();
  signal_emit(SYSTEM_SIGNAL_GS, NULL);
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
    cout << "-------------------------------------------\n";
    cout << "Creating target states for the Density Matrix\n";
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

  for(int i = 0; i < _calc_gap; i++){
    int i1, i2;

    if(_use_seed){
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

    if(verbose() > 0) cout << "<<<<<<< LANCZOS >>>>>>>\n";

    _project = true;
    double tol = _lanczos_tol;
    Vector<double> a(100);
    Vector<double> b(100);
    Vector<double> e(1);
    char file[255];
    sprintf(file,"vectors_%s.dat",_name);
    int n = _lanczos_maxiter;

    _project_states.resize(i+1);
    for(int j = 0; j <= i; j++) _project_states[j] = &_target[j];

    lanczos<T, System<T>, VectorState<T> >(*this, _state, seed, e, 1, a, b, n, tol, _use_seed, true, tmp_files.get_filename(file));
    _target[i+1] = _state[0];

    _project = false;
    for(int j = 0; j <= i; j++)
      cout << "EXCITED OVERLAP(" << j << ") = " << product(_target[j],_target[1+i]) << endl;
    cout << "E(" << i+1 << ") = " << e[0] << endl;
    energy[i+1] = e[0];
    cout << "GAP(" << i+1 << ") = " << energy[0] << " " << e[0] << " " << e[0] - energy[0] << endl;
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

  CTimer clock;
  clock.Start();

//********************************************************************
//                   DENSITY MATRIX
//********************************************************************
  signal_emit(SYSTEM_SIGNAL_TRUNCATE, NULL);
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
    BasicOp<T> this_op = term[0];
    if(dir == LEFT2RIGHT){
      if(iter > 1 && iter < ls-3){
        measure_operator(this_op.set_site(iter));
      } else if(iter == 1) {
        measure_operator(this_op.set_site(0));
        measure_operator(this_op.set_site(1));
      } else if(iter == ls-3){
        measure_operator(this_op.set_site(ls-3));
        measure_operator(this_op.set_site(ls-2));
        measure_operator(this_op.set_site(ls-1));
      }
    } else {
      if(iter > 1 && iter < ls-3){
        measure_operator(this_op.set_site(ls-2-iter));
      } else if(iter == 1) {
        measure_operator(this_op.set_site(ls-1));
        measure_operator(this_op.set_site(ls-2));
        measure_operator(this_op.set_site(ls-3));
      } else if(iter == ls-3){
        measure_operator(this_op.set_site(1));
        measure_operator(this_op.set_site(0));
      }
    }
  }

//////////////////////////////////////////////////////////////////////////
  signal_emit(SYSTEM_SIGNAL_BUILD_DM, NULL);
  if(verbose() > 0)
    cout << "Building Density Matrix in blocks \n";

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
        cout << "Target vector " << nt << " " << _target_weight[nt] << endl;
      T norm = product(s,s);
      rho += s.density_matrix(position)*T(_target_weight[nt]/norm);
    }
  }

  if(verbose() > 0)
    cout << "Lap: " << clock.LapTime().c_str() << endl;
  clock.Lap();

  T trace = T(0);
  typename BMatrix<T>::iterator biter;
  for(biter = rho.begin(); biter != rho.end(); biter++){
    SubMatrix<T> &block(*biter);
    trace += block.trace();
  }

  if(verbose() > 0)
    cout << "TRACE = " << trace << endl; 
  signal_emit(SYSTEM_SIGNAL_DM_READY, &rho);
//******* Diagonalize Density Matrix in blocks
  if(verbose() > 0)
    cout << "Diagonalizing Density Matrix in blocks\n";
  Vector<double> w(basis.dim()); // eigenvalues
  rho = rho.diagonalize(w);

  if(verbose() > 0)
    cout << "Lap: " << clock.LapTime().c_str() << endl;

//******* Calculate entropy

  _entropy = 0.0;
  for(int i = 0; i < w.size(); i++){
    if(w[i] >= 1.e-5) {
//      cout << "WARNING: w[i] <= 0 " << w[i] << endl;
      _entropy -= w[i]*log(w[i]);
    }
    if(_verbose > 0) cout << "RHO EIGENVALUE " << i << " " << w[i] << endl;
  }
//  _entropy /= LOG2;
  cout << "ITER = " << iter << " ENTROPY = " << _entropy << endl;
  

//******* Reorder eigenvalues and eigenvectors, and truncate
  if(verbose() > 0)
    cout << "Reordering eigenvalues and eigenvectors, and truncating\n";

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
        new_size = (int)((i+1)*1.1);
        if(_error_max_size != -1)
          new_size = std::min(new_size, _error_max_size);
        cout << "ERROR = " << _error << " " << error << endl;
        cout << "NEW SIZE = " << i+1 << " " << new_size << endl;
        break;
      }
    }
  }

  Vector<double> vaux = w;
  rho_basis = basis;
  int new_dim = 0;
  biter = rho.begin();
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
    for(biter = rho.begin(); biter != rho.end(); biter++){
      SubMatrix<T> &block(*biter);
      bool found = false;
      bool first = false;
      for(ib = 0; ib < bsize; ib++){
        SubMatrix<T> *_block = bvector(bidx(ib));
        if(_block == &block) { found = true; first = (ib == 0); break; }
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

  if(verbose() > 0) cout << "NEW DIM " << new_dim << endl;

  rho_basis.resize(new_dim); // truncate basis
  rho.repack(rho_basis);

//******* Check precision after truncation

  double xweight = 0;
  for(int i = 0; i < new_dim; i++){
//    cout << i << " " << w(i) << " " << vaux(idx(i)) << endl;
    xweight += w(i);
  }

  if(verbose() > 0) cout << "TOTAL WEIGHT = " << xweight << endl;

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

  cout << "-------------------------------------------\n";
  cout << "Truncation error = " << precision << endl;
  if(verbose() > 0)
    cout << "Truncation time: " << clock.TotalTime().c_str() << endl;
  cout << "-------------------------------------------\n";

  char file[255];
  sprintf(file,"iter_%s.dat",_name);
  ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::app);
  if(!outputfile) cerr << "*** ERROR: could not open " << file << endl;

  outputfile.precision(10);

  outputfile << "Truncation error = " << precision << endl;
  outputfile.close();
}

/////////////////////////////////////////////////////////////////////////
// rotate:
// build new block and operators, truncating with the density matrix.
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::rotate(int position, Block<T>& b)
{
  signal_emit(SYSTEM_SIGNAL_ROTATE, NULL);

  const Block<T> *pb1 = _b1;
  const Block<T> *pb2 = _b2;

  if(position == RIGHT){
    pb1 = _b3;
    pb2 = _b4;
  }

  CTimer clock;
  clock.Start();

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
      cout << "BLOCK " << i+1 << endl;
      typename Block<T>::const_iterator biter;
      for(biter = _b[i]->begin(); biter != _b[i]->end(); biter++){
        cout << "OPERATOR " << (*biter).name() << "(" << (*biter).site() << "," << (*biter).internal_site() << ")" << endl;
      }
    }
    cout << "-------------\n";
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
// ***************************************************************************
// Hamiltonian terms for interactions between blocks
// ***************************************************************************
  rotate_terms(position, b, basis, rho_basis, NULL);
// ***************************************************************************
// Local Hamiltonian
// ***************************************************************************
  if(verbose() > 0)
    cout << "NEW BLOCK: HAMILTONIAN TERMS: " << this->h.name() << endl;
  rotate_hami(position, b, basis, rho_basis, NULL); 
// ***************************************************************************
// Other Hamiltonians 
// ***************************************************************************
  typename std::list<Hami<T> >::iterator hiter;
  for(hiter = hs.begin(); hiter != hs.end(); hiter++){
    rotate_terms(position, b, basis, rho_basis, &(*hiter));
    rotate_hami(position, b, basis, rho_basis, &(*hiter)); 
  }
// ***************************************************************************
// Correlations
// ***************************************************************************
  rotate_corr(position, b, basis, rho_basis, &corr); 
// ***************************************************************************
  rho_basis = basis;

  if(verbose() > 0) {
    cout << "-------------------------------------------\n";
    cout << "Rotation time: " << clock.TotalTime().c_str() << endl;
    cout << "-------------------------------------------\n";
  }
  char file[255];
  sprintf(file,"iter_%s.dat",_name);
  ofstream outputfile(tmp_files.get_filename(file),std::ios::out|std::ios::app);
  if(!outputfile) cerr << "*** ERROR: could not open " << file << endl;

  outputfile << "Rotation time: " << clock.TotalTime().c_str() << endl;
  outputfile.close();
}

/////////////////////////////////////////////////////////////////////////
// rotate_tems:
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

  CTimer clock;
  clock.Start();

//********************************************************************
//                   NEW BLOCK
//********************************************************************
//******* Change of basis

  if(verbose() > 0) cout << "NEW BLOCK: SINGLE OPERATORS\n";

// ***************************************************************************
// We add to the outer block the operators necessary to build the interaction
// between blocks. We take them from the single-site blocks and we add them 
// to the outer block. 
// We also add local terms of the hamiltonian on site block
// ***************************************************************************

  typename Hami<T>::const_iterator titer;
  for(titer = this_hami.begin(); titer != this_hami.end(); titer++){
    const Term<T>& t = (*titer);
//    if(T(t.coef()) == T(0)) continue;
  
    bool found = false; 
    const BasicOp<T> *_op = NULL;

    typename Term<T>::const_iterator oiter;
    int bmask = 0;

    for(oiter = t.begin(); oiter != t.end(); oiter++){
      if(oiter->is_hami()) {
        rotate_hami(position, b, basis, rho_basis, &(oiter->hami()));
        rotate_terms(position, b, basis, rho_basis, &(oiter->hami()));
      } else {
        bmask |= mask(block(oiter->site()));
      }
    }

    if(position == LEFT && (bmask & MASK_BLOCK2) && !(bmask & MASK_BLOCK1) && bmask != MASK_BLOCK2)
       found = true;
    if(position == RIGHT && (bmask & MASK_BLOCK3) && !(bmask & MASK_BLOCK4) && bmask != MASK_BLOCK3)
       found = true;

    if(found){
      BasicOp<T> real_op;
      for(oiter = t.begin(); oiter != t.end(); oiter++){
        if(block(oiter->site()) == site_block && !oiter->is_hami()){
          _op = operator()(*oiter);
          if(!_op) cout << "WARNING : Operator " << (*oiter).description() << " not found\n";
          real_op = *oiter;
          real_op.dqn = _op->dqn;
          break;
        }
      }
      bool add = b.contains(real_op);

      if(!add && _op){ // If we haven't added it yet, we do it now
        real_op.resize(rho_basis);

        if(verbose() > 0)
          cout << "ADDING NEW OPERATOR " << real_op.description() << endl;
        clock.Lap();

        new_operator(real_op, *_op, rho, basis, 1-position);
        b.push_back(real_op);
        if(verbose() > 0)
          cout << "Lap: " << clock.LapTime().c_str() << endl;
      }
    }
  }

// ***************************************************************************
// Composite operators for interactions 
// ***************************************************************************
// 1- We look for terms that involve pieces in one block and pieces
//    in another block (not a site block) and we rotate them 
//    We also add TERM_LOCAL and TERM_EXTERN that are in the block

  for(titer = this_hami.begin(); titer != this_hami.end(); titer++){
    const Term<T>& t = (*titer);
    bool found = false;
    if((t.type() == TERM_PRODUCT && t.size() >= 2) || t.type() == TERM_LOCAL || t.type() == TERM_EXTERN){
      Term<T> aux_term;
      BasicOp<T> top2;
      int bmask = 0;
      for(int i = 0; i < t.size(); i++){
        BasicOp<T> top = t[i];
        if(block(top.site()) == the_block){
          aux_term *= top;
          found = true;
        }
        bmask |= mask(block(top.site()));
      }

      bool doit = false;
      if(found && position == LEFT && (bmask & MASK_BLOCK1) && ((!(bmask & MASK_BLOCK2) && bmask != MASK_BLOCK1) || (bmask == MASK_BLOCK1 && (t.type() == TERM_LOCAL || t.type() == TERM_EXTERN))))
         doit = true;
      if(found && position == RIGHT && (bmask & MASK_BLOCK4) && ((!(bmask & MASK_BLOCK3) && bmask != MASK_BLOCK4) || (bmask == MASK_BLOCK4 && (t.type() == TERM_LOCAL || t.type() == TERM_EXTERN))))
         doit = true;

      if(doit){ // we found a piece of a composite operator
        BasicOp<T> top1;
        if(aux_term.size() == 1) // the piece contains a single operator
          top1 = aux_term[0];
        else
          top1 = aux_term; // the piece contains more than one operator

        const BasicOp<T>* op1 = operator()(top1);

        if(!op1) {
          cout << "ERROR 1: Operator " << top1.description() << " not found\n";
          continue;
        }

        BasicOp<T> new_op(top1); 
        new_op.dqn = op1->dqn;
        bool add = b.contains(new_op);

        if(!add && op1){
          new_op.resize(rho_basis);
          clock.Lap();
  
          if(verbose() > 0)
            cout << "NEW OPERATOR " << new_op.description() << endl;
          new_op = T(0);
          new_operator(new_op, *op1, rho, basis, position);
          b.push_back(new_op);
          if(verbose() > 0)
            cout << "Lap: " << clock.LapTime().c_str() << endl;
        }
      }
    }
  }

// 2- We look for terms that involve pieces in one block and pieces
//    in a site block and we put the pieces together

  for(titer = this_hami.begin(); titer != this_hami.end(); titer++){
    const Term<T>& t = (*titer);
    bool found1 = false;
    bool found2 = false;
    if(t.size() >= 2){
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
//    if bmask == _mask, this term will go into the Hamiltonian
//    unless it is TERM_EXTERN
      if(((bmask != _mask) || (bmask == _mask && t.type() == TERM_EXTERN)) && found1 && found2){ // we have a new composite operator
        BasicOp<T> top1(aux_term1);
        if(aux_term1.size() == 1)
          top1 = aux_term1[0];
        BasicOp<T> top2(aux_term2);
        if(aux_term2.size() == 1)
          top2 = aux_term2[0];
        const BasicOp<T>* op1 = operator()(top1);
        const BasicOp<T>* op2 = operator()(top2);

        if(!op1) {
          cout << "ERROR 2: Operator " << top1.description() << " in term " << t.description() << " not found\n";
          continue;
        }
        if(!op2) {
          cout << "ERROR 2: Operator " << top2.description() << " in term " << t.description() << " not found\n";
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
        bool add = b.contains(new_op);

        if(!add && op1 && op2){
          new_op.resize(rho_basis);
          clock.Lap();
  
          if(verbose() > 0)
            cout << "ADDING NEW OPERATOR " << new_op.description() << endl;
          new_op = T(0);
          new_operator(new_op, *op1, *op2, pos1, pos2, rho, basis, T(1));
          b.push_back(new_op);
          if(verbose() > 0)
            cout << "Lap: " << clock.LapTime().c_str() << endl;
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

  CTimer clock;
  clock.Start();


  const Hami<T> *ph = &corr;
  if(this_hami) ph = this_hami;
  const Hami<T> &_h = *ph;

// 1a-Look for operators in S.ops involving the other blocks and site block 
// We also rotate individual operators

  typename Hami<T>::const_iterator titer;
  if(!this_hami){
    for(titer = ops.begin(); titer != ops.end(); titer++){
      const Term<T>& t = (*titer);
 
      bool found = false; 

      typename Term<T>::const_iterator oiter;
      int bmask = 0;
      for(oiter = t.begin(); oiter != t.end(); oiter++){
        if(!oiter->is_hami()) {
          bmask |= mask(block(oiter->site()));
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
  
        const BasicOp<T>* _op = operator()(top1);
  
        if(!_op) {
          cout << "ERROR 1a: Operator " << top1.description() << " not found\n";
          continue;
        }
  
        BasicOp<T> real_op(top1); 
        real_op.dqn = _op->dqn;
        bool add = b.contains(real_op);
  
        if(!add && _op){ // If we haven't added it yet, we do it now
          real_op.resize(rho_basis);
  
          if(verbose() > 0)
            cout << "ADDING NEW MEAS. OPERATOR " << real_op.description() << endl;
          clock.Lap();
  
          new_operator(real_op, *_op, rho, basis, 1-position);
          b.push_back(real_op);
          if(verbose() > 0)
            cout << "Lap: " << clock.LapTime().c_str() << endl;
        }
      }
    }
  } 

// 1b-Look for operators in S.corr involving the other blocks and site block 
// We also rotate individual operators

  for(titer = _h.begin(); titer != _h.end(); titer++){
    const Term<T>& t = (*titer);
//    if(T(t.coef()) == T(0)) continue;
 
/* 
    if(t.type() == TERM_PRODUCT && t.size() == 2 && t[1].is_hami()){
      rotate_additive(position, b, basis, rho_basis, t); 
    }
*/
 
    bool found = false; 

    typename Term<T>::const_iterator oiter;
    int bmask = 0;
    for(oiter = t.begin(); oiter != t.end(); oiter++){
      if(!oiter->is_hami()) {
        bmask |= mask(block(oiter->site()));
      } else {
        rotate_terms(position, b, basis, rho_basis, &(oiter->hami()));
        rotate_hami(position, b, basis, rho_basis, &(oiter->hami()));
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
      if(aux_term.size() == t.size()) aux_term.coef() = t.coef();

      BasicOp<T> top1;
      if(aux_term.size() == 1) // the piece contains a single operator
        top1 = aux_term[0];
      else
        top1 = aux_term; // the piece contains more than one operator

      const BasicOp<T>* _op = operator()(top1);

      if(!_op) {
        cout << "ERROR 1b: Operator " << top1.description() << " not found\n";
        continue;
      }

      BasicOp<T> real_op(top1); 
      real_op.dqn = _op->dqn;
      bool add = b.contains(real_op);

      if(!add && _op){ // If we haven't added it yet, we do it now
        real_op.resize(rho_basis);

        if(verbose() > 0)
          cout << "ADDING NEW MEAS. OPERATOR " << real_op.description() << endl;
        clock.Lap();

        new_operator(real_op, *_op, rho, basis, 1-position);
        b.push_back(real_op);
        if(verbose() > 0)
          cout << "Lap: " << clock.LapTime().c_str() << endl;
      }
    }
  }

// 2a- We look for terms that involve pieces in old block and pieces
// in another block (not site blocks) and we rotate them 
// We also rotate individual operators

  if(!this_hami){
    for(titer = ops.begin(); titer != ops.end(); titer++){
      const Term<T>& t = (*titer);
      bool found = false;
  
      typename Term<T>::const_iterator oiter;
      int bmask = 0;
      for(oiter = t.begin(); oiter != t.end(); oiter++){
        if(!oiter->is_hami()) {
          bmask |= mask(block(oiter->site()));
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
  
        const BasicOp<T>* _op = operator()(top1);
  
        if(!_op) {
          cout << "ERROR 2a: Operator " << top1.description() << " not found\n";
          continue;
        }
  
        BasicOp<T> new_op(top1); 
        new_op.dqn = _op->dqn;
        bool add = b.contains(new_op);
  
        if(!add && _op){
          new_op.resize(rho_basis);
          clock.Lap();
  
          if(verbose() > 0)
            cout << "NEW MEAS. OPERATOR " << new_op.description() << endl;
          new_op = T(0);
          new_operator(new_op, *_op, rho, basis, position);
          b.push_back(new_op);
          if(verbose() > 0)
            cout << "Lap: " << clock.LapTime().c_str() << endl;
        }
      }
    }
  }

// 2b- We look for terms that involve pieces in old block and pieces
// in another block (not site blocks) and we rotate them 
// We also rotate individual operators

  for(titer = _h.begin(); titer != _h.end(); titer++){
    const Term<T>& t = (*titer);
    bool found = false;

    typename Term<T>::const_iterator oiter;
    int bmask = 0;
    for(oiter = t.begin(); oiter != t.end(); oiter++){
      if(!oiter->is_hami()) {
        bmask |= mask(block(oiter->site()));
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

      const BasicOp<T>* _op = operator()(top1);

      if(!_op) {
        cout << "ERROR 2b: Operator " << top1.description() << " not found\n";
        continue;
      }

      BasicOp<T> new_op(top1); 
      new_op.dqn = _op->dqn;
      bool add = b.contains(new_op);

      if(!add && _op){
        new_op.resize(rho_basis);
        clock.Lap();

        if(verbose() > 0)
          cout << "NEW MEAS. OPERATOR " << new_op.description() << endl;
        new_op = T(0);
        new_operator(new_op, *_op, rho, basis, position);
        b.push_back(new_op);
        if(verbose() > 0)
          cout << "Lap: " << clock.LapTime().c_str() << endl;
      }
    }
  }

// 3- Composite operators for correlations
  if(_store_products){
    for(titer = _h.begin(); titer != _h.end(); titer++){
      const Term<T>& t = (*titer);
      bool found1 = false;
      bool found2 = false;
      if(t.type() == TERM_PRODUCT && t.size() >= 2){
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
          const BasicOp<T>* op1 = operator()(top1);
          const BasicOp<T>* op2 = operator()(top2);
  
          if(!op1) {
            cout << "ERROR 3: Operator " << top1.description() << " not found\n";
            continue;
          }
          if(!op2) {
            cout << "ERROR 3: Operator " << top2.description() << " not found\n";
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
          bool add = b.contains(new_op);
  
          if(!add && op1 && op2){
            new_op.resize(rho_basis);
            clock.Lap();
    
            if(verbose() > 0)
              cout << "ADDING NEW CORR. TERM " << new_op.description() << endl;
            new_op = T(0);
            new_operator(new_op, *op1, *op2, pos1, pos2, rho, basis, T(1));
  
            b.push_back(new_op);
            if(verbose() > 0)
              cout << "Lap: " << clock.LapTime().c_str() << endl;
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

  const Block<T> *pb1 = _b1;
  const Block<T> *pb2 = _b2;
  int _mask = MASK_BLOCK1 | MASK_BLOCK2;
  size_t pos1 = BLOCK1;
  size_t pos2 = BLOCK2;

  if(position == RIGHT){
    pb1 = _b3;
    pb2 = _b4;
    _mask = MASK_BLOCK3 | MASK_BLOCK4;
    pos1 = BLOCK3;
    pos2 = BLOCK4;
  }

  const Block<T> &b1(*pb1);
  const Block<T> &b2(*pb2);

  CTimer clock;
  clock.Start();

  if(verbose() > 0)
    cout << "HAMILTONIAN TERMS: " << this_hami.name() << endl;

  BasicOp<T> _h = BasicOp<T>(this_hami);
  _h.resize(rho_basis);

  bool found2 = false;
  bool found1 = false;
  const BasicOp<T> *local_op2 = b2(BasicOp<T>(this_hami));
  if(local_op2 && local_op2->name() == this_hami.name()){ 
    if(verbose() > 0)
      cout << "OPERATOR " << this_hami.name() << "2" << endl;
    clock.Lap();
    found2 = true;
    new_operator(_h, *local_op2, rho, basis, RIGHT);
    if(verbose() > 0)
      cout << "Lap: " << clock.LapTime().c_str() << endl;
  } else {
//    if(this_hami.name() == "H")
//      cout << "WARNING: Add local (site) Hamiltonian term\n";
//    else
//      cout << "You can ignore the warning\n";
  }

  const BasicOp<T> *local_op1 = b1(BasicOp<T>(this_hami));
  if(local_op1 && local_op1->name() == this_hami.name()){ 
    if(verbose() > 0)
      cout << "OPERATOR " << this_hami.name() << "1" << endl;
    clock.Lap();
    found1 = true;
    new_operator(_h, *local_op1, rho, basis, LEFT);
    if(verbose() > 0)
      cout << "Lap: " << clock.LapTime().c_str() << endl;
  } else {
//    if(this_hami.name() == "H")
//      cout << "WARNING: Add local (site) Hamiltonian term\n";
//    else
//      cout << "You can ignore the warning\n";
  }

  typename Hami<T>::const_iterator hiter;
  for(hiter = this_hami.begin(); hiter != this_hami.end(); hiter++){
    const Term<T>& t = (*hiter);
//    if(T(t.coef()) == T(0)) continue;
    if(t.type() == TERM_EXTERN) continue;
    if(t.type() == TERM_LOCAL) continue;
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
            cout << "TERM " << t.name(true) << endl;
          else
            cout << "TERM " << t.name(true) << " + h.c." << endl;
        }
                                                                                
        clock.Lap();
        new_operator(_h, *op, rho, basis, 1-position, T(t.coef()));
        if(verbose() > 0)
          cout << "Lap: " << clock.LapTime().c_str() << endl;
      }
//if we don't have a site Hamiltonian we try to find the missing term
      if((ib == BLOCK1 && position == LEFT && !found1) || (ib == BLOCK4 && position == RIGHT && !found2) || (position == LEFT && ib == BLOCK1 && b1.lattice().size() == 1) || (position == RIGHT && ib == BLOCK4 && b2.lattice().size() == 1)){

        const BasicOp<T>* op= operator()(top);
        bool calc_hc = this_hami.use_hc();
        if (op->is_diagonal()) calc_hc = false;
                      
        if(verbose() > 0) {
          if(!calc_hc)
            cout << "MISSING TERM " << t.name(true) << endl;
          else
            cout << "MISSING TERM " << t.name(true) << " + h.c." << endl;
        }
                                                                                
        clock.Lap();
        new_operator(_h, *op, rho, basis, position, T(t.coef()));
        if(verbose() > 0)
          cout << "Lap: " << clock.LapTime().c_str() << endl;
      }
    }
// Multiple products
    if(t.size() >= 2 && t.type() == TERM_PRODUCT){
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
          cout << "ERROR : Operator " << top1.description() << " not found\n";
          continue;
        }
        if(!op2) {
          cout << "ERROR : Operator " << top2.description() << " not found\n";
          continue;
        }

        bool calc_hc = this_hami.use_hc();
        if (op1->is_diagonal() && op2->is_diagonal()) calc_hc = false;
        if (t.size() > 2) calc_hc = false; // you have to add the h.c. terms explicitly

        if(op1 && op2){
          clock.Lap();
 
          if(verbose() > 0){
            if(!calc_hc)
              cout << "TERM " << t.description() << endl;
            else
              cout << "TERM " << t.description() << " + h.c." << endl;
          }
          new_operator(_h, *op1, *op2, pos1, pos2, rho, basis, T(t.coef()), calc_hc);
         
          if(verbose() > 0) 
            cout << "Lap: " << clock.LapTime().c_str() << endl;
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
System<T>::measure_n(size_t n , const VectorState<T> *v )
{

  using namespace std;

  const VectorState<T> *pv = NULL;
  if(!v)
   pv = &gs;
  else
   pv = v;

  VectorState<T> res(*pv);
  VectorState<T> aux(*pv);
  T x = sqrt(product(*pv,*pv));
  T x0 = sqrt(product(gs,gs));

  typename Hami<T>::iterator titer;

  for(titer = corr.begin(); titer != corr.end(); titer++){
    Term<T>& t = *titer;

    if(n != 0 && t.size() != n) continue;

    aux = *pv;

    typename Term<T>::iterator iter;
    iter = t.end();
    iter--;
    while(true){
      BasicOp<T> top = *iter;

      size_t ib = block(top.site());

      const BasicOp<T>* op = operator()(top);
      if(op){ 
        res.set_qn_mask(aux.qn()+op->dqn, _grand_canonical);
        res.resize(aux.b1(),aux.b2(),aux.b3(),aux.b4());
        res = T(0);

        product(*op, aux, res, ib);
        aux = res;
      }

      if(iter == t.begin()) break;
      iter--;
    }
    T val = T(0);
    if(gs.qn().equal(res.qn(), _grand_canonical) || _grand_canonical == 0)
      val = t.coef()*product(gs,res)/x/x0;

/*
    cout << setprecision(10) << t.coef();
    for(iter = t.begin(); iter != t.end(); iter++){
      BasicOp<T> &top = *iter;
      if(h.lattice().type() == LATTICE_2D)
        cout << "*" << top.name().c_str() << "(" << h.lattice().x(top.site()) << "," << h.lattice().y(top.site()) << "," << top.internal_site() << ")";
      if(h.lattice().type() == LATTICE_1D)
        cout << "*" << top.name().c_str() << "(" << h.lattice().x(top.site()) << "," << top.internal_site() << ")";
    }
    cout << " = " << val << endl;
*/

//// cout << setprecision(10) << t.name(true) << " = " << val << endl;
    t.set_value(val);
  }

}

template<class T>
void
System<T>::measure(const VectorState<T> *v)
{

  using namespace std;
  //DMTKglobals<T> &globals = get_globals(T(0));

  const VectorState<T> *pv = NULL;
  if(!v)
   pv = &gs;
  else
   pv = v;

  VectorState<T> res(*pv);
  VectorState<T> aux(*pv);
  VectorState<T> aux_res(*pv);
  T x = product(*pv,*pv);
  T x0 = product(*pv,*pv);


  typename Hami<T>::iterator titer;
//////////////////////////////////////////////////////////////////////////
  for(titer = corr.begin(); titer != corr.end(); titer++){
    Term<T> t = *titer;
    Term<T> &real_t = *titer;

    if(t.size() == 1 && t[0].is_hami()) {
      res = T(0);
      res = product_default(*this, *pv, &t[0].hami());
      T val = product(*pv,res)/x;  
//      cout << setprecision(10) << t[0].hami().name() << " = " << val << endl;
      real_t.set_value(val);
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
            
          const BasicOp<T>* op1 = operator()(top1);
          const BasicOp<T>* op2 = operator()(top2);
  
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
  
////          cout << setprecision(10) << t.name(true) << " = " << val << endl;
  
            real_t.set_value(val);
          }
  
        } else {

          BasicOp<T> aux_op(t);
          const BasicOp<T>* op = operator()(aux_op);
  
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
  
////          cout << setprecision(10) << t.name(true) << " = " << val << endl;

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
      int bmask = 0;
      for(int i = 0; i < t.size(); i++){
        BasicOp<T> top = t[i];
        if(block(top.site()) == BLOCK_NONE){
          found = false;
          break;
        }
        for(int ib = 1; ib <= 4; ib++){
          if(block(top.site()) == (size_t)ib){
            aux_term[ib-1] *= top;
            bmask |= mask((size_t)ib);
          }
        }
      } 
  
      if(!found) continue;
  
      const BasicOp<T> *aux_op[4];
      aux_op[0] = aux_op[1] = aux_op[2] = aux_op[3] = NULL;
      for(int ib = 3; ib >= 0 ; ib--){
        if(aux_term[ib].size() != 0){
          BasicOp<T> top;
          if(aux_term[ib].size() == 1)
            top = aux_term[ib][0];
          else
            top = aux_term[ib];
  
          aux_op[ib] = operator()(top);
          if(!aux_op[ib]) {
            cout << "ERROR : Operator " << top.description() << " not found\n";
          }
        }
      }

      aux = aux_res = *pv;
      for(int ib = 3; ib >= 0 ; ib--){
        if(aux_op[ib]){
          res.set_qn_mask(aux.qn()+aux_op[ib]->dqn, _grand_canonical);
          res.resize(aux.b1(),aux.b2(),aux.b3(),aux.b4());
          res = T(0);

          product(*aux_op[ib], aux, res, (size_t)(ib+1));
          aux = res;
        }
      }
      T val = product(*pv,res)/sqrt(x)/sqrt(x0);
      real_t.set_value(val*t.coef());
//      cout << setprecision(10) << t.description() << " = " << real_t.value() << endl;
    }
  }
}

template<class T>
T
System<T>::measure_operator(const BasicOp<T> &top, const VectorState<T> *v)
{
  using namespace std;

  const VectorState<T> *pv = NULL;
  if(!v)
   pv = &gs;
  else
   pv = v;

  VectorState<T> res(*pv);
  VectorState<T> aux(*pv);
  T x = sqrt(product(*pv,*pv));
  T x0= sqrt(product(gs,gs));

  size_t ib = block(top.site());

  cout << "MEASURE OPERATOR " << top.name() << endl;
  const BasicOp<T>* op = operator()(top);
  if(op->name() != top.name()) return T(0);

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
  if(hami == NULL)
    phami = &ss.h;
  else
    phami = hami;
  const Hami<T> &h = *phami;

/////////////////////////////////////////////////////////////
  typename Hami<T>::const_iterator hiter;
  VectorState<T> aux(res), aux_res(res);
  for(hiter = h.begin(); hiter != h.end(); hiter++){
    const Term<T>& t = (*hiter);
//    if(T(t.coef()) == T(0)) continue;   
    if(t.type() == TERM_EXTERN && !ss.apply_extern()) continue;
    if(t.size() > 2){  // you have to include h.c. terms explicitly!
      Term<T> aux_term[4];

      if(ss.verbose() > 1)
        cout << "PRODUCT TERM " << t.description() << endl;

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
          if(ss.block(top.site()) == (size_t)ib){
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
              cout << "ERROR : Operator " << top.description() << " not found\n";
            }
    
            if(_op){
              if(ss.verbose() > 1)
                cout << "APPLYING " << top.description() << endl;
    
              aux_res.set_qn_mask(aux.qn()+_op->dqn, ss.qn_mask());
              aux_res.resize(aux.b1(),aux.b2(),aux.b3(),aux.b4());
              aux_res = T(0);
    
              product(*_op, aux, aux_res, (size_t)(ib+1));
              aux = aux_res;
            }
          }
        }
        res += aux_res * t.coef(); 
      }
    }  
    if(t.size() == 2 && t.type() == TERM_PRODUCT){
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
            cout << "PRODUCT TERM " << t.description() << endl;
//" " << top1.description() << " " << top1.site() << " " << top2.site() << " " << b1 << " " << b2 << endl;

          if(op1 && op2) {
            if(ss.verbose() > 0) {
              COUT_PRODUCT(t);
            }
            product(*op2, *op1, vv, res, b2, b1, T(t.coef()), calc_hc);
          }

        } else if(t.type() == TERM_EXTERN && ss.apply_extern()) {

          if(only_local && (b1 == BLOCK1 || b1 == BLOCK4)) continue;

          BasicOp<T> aux_op(t);
          const BasicOp<T>* op = ss(aux_op);

          if(ss.verbose() > 1) {
            if(calc_hc) 
              cout << "PRODUCT EXTERN " << t.coef() << "*" << top1.name() << "(" << top1.site() << ")" << top2.name() << "(" << top2.site() << ") + h.c." << endl;
            else 
              cout << "PRODUCT EXTERN " << t.coef() << "*" << top1.name() << "(" << top1.site() << ")" << top2.name() << "(" << top2.site() << ")" << endl;
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

//      if(((b1 == BLOCK1 && ss.leftblock.lattice().size() != 1) || (b1 == BLOCK4 && ss.rightblock.lattice().size() != 1)) && t.type() != TERM_LOCAL && t.type() != TERM_EXTERN) continue;
//      if(((b1 == BLOCK1 && ss.leftblock.lattice().size() != 1) || (b1 == BLOCK4 && ss.rightblock.lattice().size() != 1))) continue;
//      if((b1 == BLOCK1 || b1 == BLOCK4) && t.type() != TERM_LOCAL && t.type() != TERM_EXTERN) continue;


      if(((b1 == BLOCK1 && ss.leftblock.lattice().size() != 1) || (b1 == BLOCK4 && ss.rightblock.lattice().size() != 1))) continue;

      if(b1 != BLOCK_NONE && top.name() != h.name() && !top.is_hami()){
        const BasicOp<T>* op1 = ss(top);
        bool calc_hc = h.use_hc();
        if (op1 && op1->is_diagonal()) calc_hc = false;

        if(ss.verbose() > 1) {
          if(calc_hc)
            cout << "PRODUCT LOCAL " << t.coef() << "*" << top.name() << "(" << top.site() << ") + h.c." << endl;
          else
            cout << "PRODUCT LOCAL " << t.coef() << "*" << top.name() << "(" << top.site() << ")" << endl;
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
      cout << "PRODUCT " << name << "4\n";
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
      cout << "PRODUCT " << name << "3\n";
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
      cout << "PRODUCT " << name << "2\n";
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
      cout << "PRODUCT " << name << "1\n";
    product(*local_op1, vv, res, BLOCK1);
  } else {
//    if(h.name() == "H")
//      cout << "WARNING: Add local (site) Hamiltonian term\n";
//    else
//      cout << "You can ignore the warning\n";
  }

  return res;
}

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

template <class T>
size_t 
System<T>::block(int site) const
{
  typename Lattice::const_iterator iter;
  std::vector<const Block<T>* > b(4);
  std::vector<int> offset(4);
  b[0] = _b1;
  b[1] = _b2;
  b[2] = _b3;
  b[3] = _b4;
  offset[0] = _in_warmup && _grow_symmetric && _grow_outward && size() < h.lattice().size() ? h.lattice().size()/2 - _b1->lattice().size() - _b2->lattice().size() : 0;
  offset[1] = _in_warmup && _grow_symmetric && _grow_outward && size() < h.lattice().size() ? offset[0] + _b1->lattice().size() : _b1->lattice().size();
  offset[2] = offset[1] + _b2->lattice().size();
  offset[3] = offset[2] + _b3->lattice().size();

  if(_in_warmup && _grow_symmetric && size() < h.lattice().size() && !_grow_outward){
    offset[3] = h.lattice().size() - _b4->lattice().size();
    offset[2] = offset[3] - _b3->lattice().size();
  }

  for(size_t i = 0; i < 4; i++){
    const Lattice &l = b[i]->lattice();
    int n = 0;
    for(iter = l.begin(); iter != l.end(); iter++){
        if((n++)+offset[i] == site){
           return (size_t)(i+1); 
        }
    }
  }
  return BLOCK_NONE;
}

template<class T>
const BasicOp<T>* 
System<T>::operator()(const BasicOp<T>& op) const
{
  typename Lattice::const_iterator iter;
  std::vector<const Block<T>* > b(4);
  std::vector<int> offset(4);
  b[0] = _b1;
  b[1] = _b2;
  b[2] = _b3;
  b[3] = _b4;
  offset[0] = _in_warmup && _grow_symmetric && _grow_outward && size() < h.lattice().size() ? h.lattice().size()/2 - _b1->lattice().size() - _b2->lattice().size() : 0;
  offset[1] = _in_warmup && _grow_symmetric && _grow_outward && size() < h.lattice().size() ? offset[0] + _b1->lattice().size() : _b1->lattice().size();
  offset[2] = offset[1] + _b2->lattice().size();
  offset[3] = offset[2] + _b3->lattice().size();

  if(_in_warmup && _grow_symmetric && size() < h.lattice().size() && !_grow_outward){
    offset[3] = h.lattice().size() - _b4->lattice().size();
    offset[2] = offset[3] - _b3->lattice().size();
  }
  for(int i = 0; i < 4; i++){
    const Lattice &l = b[i]->lattice();
    int n = 0;
    for(iter = l.begin(); iter != l.end(); iter++){
        if((n++)+offset[i] == op.site()){
           BasicOp<T> _op(op);
           if(b[i]->single_site()) _op.set_site(op.site() - offset[i]);
           
// cout << op.name() << " " << op.site() << " " << op.internal_site() << " " << i << " " <<  _op.site() << endl;
           return (b[i]->operator()(_op));
        }
    }
  }
  return 0;
}

//////////////////////////////////////////////////////////////////

} // namespace dmtk

#endif // __DMTK_SYSTEM_H__
