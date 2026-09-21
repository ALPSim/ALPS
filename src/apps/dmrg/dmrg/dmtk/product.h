#ifndef __DMTK_PRODUCT_H__
#define __DMTK_PRODUCT_H__

namespace dmtk
{

template <class T>
struct ProductTerm
{
    T coef;
    int nops;
    bool do_hc;
    int m1;
    int m2;
    int m3;
    int m4;
    const SubMatrix<T> *block1;
    const SubMatrix<T> *block2;
    const SubMatrix<T> *block3;
    const SubMatrix<T> *block4;
    StateSpace vspace;
    StateSpace res_space;
    StateSpace vspace_hc;
    StateSpace res_space_hc;
};

////////////////////////////////////////////////////////////////////
// get product terms
////////////////////////////////////////////////////////////////////
template<class T>
std::vector<ProductTerm<T> >
get_product_terms(const BasicOp<T> &op,
        const VectorState<T> &v, VectorState<T> &res,
        size_t m, T coef = T(1), bool hc = false)
{
  std::vector<ProductTerm<T> > product_terms;
  Vector<QN> dqn(5);
  if(!hc) 
    dqn(m) += op.dqn; 
  else 
    dqn(m) -= op.dqn;

  bool do_hc = hc;
  T real_coef = (!hc) ? coef : dmtk::conj(coef);

  const SubMatrix<T> *_block = NULL;

  typename vector<SubMatrix<T> *>::iterator biter;
  typename VectorState<T>::const_iterator siter;
  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++){
    StateSpace ss = *siter;

    if(!hc){
        _block = op.block(ss[m].qn());
    }else{
        _block = op.block(ss[m].qn()-op.dqn);
    }

    if(!_block) continue;

    cstate_slice<T> v_slice = v(ss);
    state_slice<T> res_slice = res(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    if(res_slice.size() == 0) continue; 

    int sign = 1;
    if(op.fermion()){
      for(int ib = m-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign();
    }

    ProductTerm<T> pterm;
    pterm.coef = T(sign)*real_coef;
    pterm.nops = 1;
    pterm.m1 = m;
    pterm.block1 = _block;
    pterm.vspace = ss;
    pterm.do_hc = do_hc;
    pterm.res_space = res.get_qn_space(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2], ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    product_terms.push_back(pterm); 
  }
  return product_terms;
}

template<class T>
std::vector<ProductTerm<T> >
get_product_terms(const BasicOp<T> &op1, const BasicOp<T>& op2,
        const VectorState<T> &v, VectorState<T> &res,
        size_t m1, size_t m2, T coef = T(1), bool hc = false)
{
  std::vector<ProductTerm<T> > product_terms;
  int _mask = mask(m1,m2);
  int _m1, _m2;
  Vector<QN> dqn(5);
  dqn(m2) += op2.dqn;
  dqn(m1) += op1.dqn;
  Vector<QN> dqn2(5);
  dqn2(m2) += op2.dqn;

  const BasicOp<T> *_op1;
  const BasicOp<T> *_op2;
  int sign0 = 1;
  if(m1 <= m2){
    _op1 = &op1;
    _op2 = &op2;
    _m1 = m1;
    _m2 = m2;
  } else {
    _op1 = &op2;
    _op2 = &op1;
    _m1 = m2;
    _m2 = m1;
    if(_op1->fermion() && _op2->fermion()) sign0 = -1;
  }

  const SubMatrix<T> *_block1 = NULL;
  const SubMatrix<T> *_block2 = NULL;

  typename VectorState<T>::const_iterator siter;
  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++){
    StateSpace ss = *siter;


    _block2 = _op2->block(ss[_m2].qn());
    if(m1 != m2){
        _block1 = _op1->block(ss[_m1].qn());
    } else {
        _block1 = _op1->block(ss[_m1].qn() + _op2->dqn);
    }
    if(!_block1 || !_block2) continue;

    cstate_slice<T> v_slice = v(ss);
    state_slice<T> res_slice = res(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    if(res_slice.size() == 0) continue; 

    cstate_slice<T> v_slice_hc = v(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    state_slice<T> res_slice_hc = res(ss[1].qn(),ss[2].qn(),ss[3].qn(),ss[4].qn());
    if(v_slice_hc.size() == 0 || res_slice_hc.size() == 0) continue; 

    int sign = sign0;

    if(op2.fermion()){
      for(int ib = m2-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign();
    }
    if(op1.fermion()){
      for(int ib = m1-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn2[ib].fermion_sign();
    }

    ProductTerm<T> pterm;
    pterm.coef = T(sign)*coef;
    pterm.nops = 2;
    pterm.block1 = _block1;
    pterm.block2 = _block2;
    pterm.m1 = _m1;
    pterm.m2 = _m2;
    pterm.do_hc = hc;
    pterm.vspace = ss;
    pterm.res_space = res.get_qn_space(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2], ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    pterm.vspace_hc = v.get_qn_space(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2], ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    pterm.res_space_hc = ss;
    product_terms.push_back(pterm); 
  }
  return product_terms;
}

template<class T>
std::vector<ProductTerm<T> >
get_product_terms(const BasicOp<T> &op1, const BasicOp<T>& op2,
        const BasicOp<T> &op3, const BasicOp<T>& op4,
        const VectorState<T> &v, VectorState<T> &res,
        size_t m1, size_t m2, size_t m3, size_t m4,
        T coef = T(1), bool hc = false)
{
  std::vector<ProductTerm<T> > product_terms;
  Vector<QN> dqn(5);
  size_t _m1, _m2, _m3, _m4;
  dqn(m4) += op4.dqn;
  dqn(m3) += op3.dqn;
  dqn(m2) += op2.dqn;
  dqn(m1) += op1.dqn;
  Vector<QN> dqn2(5);
  dqn2(m2) += op2.dqn;
  Vector<QN> dqn3(5);
  dqn3(m3) += op3.dqn;
  Vector<QN> dqn4(5);
  dqn4(m4) += op4.dqn;

  const BasicOp<T> *_op1;
  const BasicOp<T> *_op2;
  const BasicOp<T> *_op3;
  const BasicOp<T> *_op4;

  Vector<size_t> m(4), indx(4);
  Vector<const BasicOp<T>* > ops(4);
  ops(0) = &op1;
  ops(1) = &op2;
  ops(2) = &op3;
  ops(3) = &op4;
  m(0) = int(m1);
  m(1) = int(m2);
  m(2) = int(m3);
  m(3) = int(m4);
  indexx<size_t,Vector<size_t> >(4, m, indx);
  _m1 = (size_t)m(indx(0));
  _m2 = (size_t)m(indx(1));
  _m3 = (size_t)m(indx(2));
  _m4 = (size_t)m(indx(3));
  _op1 = ops(indx(0));
  _op2 = ops(indx(1));
  _op3 = ops(indx(2));
  _op4 = ops(indx(3));

  const SubMatrix<T> *_block1 = NULL;
  const SubMatrix<T> *_block2 = NULL;
  const SubMatrix<T> *_block3 = NULL;
  const SubMatrix<T> *_block4 = NULL;

  typename VectorState<T>::const_iterator siter;
  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++){
    StateSpace ss = *siter;

      _block4 = _op4->block(ss[_m4].qn());
      _block3 = _op3->block(ss[_m3].qn());
      _block2 = _op2->block(ss[_m2].qn());
      _block1 = _op1->block(ss[_m1].qn());

    if(!_block1 || !_block2 || !_block3 || !_block4) continue;

    cstate_slice<T> v_slice = v(ss);
    state_slice<T> res_slice = res(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    if(res_slice.size() == 0) continue; 

    cstate_slice<T> v_slice_hc = v(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    state_slice<T> res_slice_hc = res(ss[1].qn(),ss[2].qn(),ss[3].qn(),ss[4].qn());
    if(hc && (v_slice_hc.size() == 0 || res_slice_hc.size() == 0)) continue; 

    int sign = 1;
    if(op4.fermion()){
      for(int ib = m4-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign();
    }
    if(op3.fermion()){
      for(int ib = m3-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn4[ib].fermion_sign();
    }
    if(op2.fermion()){
      for(int ib = m2-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn3[ib].fermion_sign()*dqn4[ib].fermion_sign();
    }
    if(op1.fermion()){
      for(int ib = m1-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn2[ib].fermion_sign()*dqn3[ib].fermion_sign()*dqn4[ib].fermion_sign();
    }


    ProductTerm<T> pterm;
    pterm.coef = T(sign)*coef;
    pterm.nops = 4;
    pterm.block1 = _block1;
    pterm.block2 = _block2;
    pterm.block3 = _block3;
    pterm.block4 = _block4;
    pterm.m1 = _m1;
    pterm.m2 = _m2;
    pterm.m3 = _m3;
    pterm.m4 = _m4;
    pterm.do_hc = hc;
    pterm.vspace = ss;
    pterm.res_space = res.get_qn_space(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2], ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    pterm.vspace_hc = v.get_qn_space(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2], ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    pterm.res_space_hc = ss;
    product_terms.push_back(pterm); 
  }

  return product_terms;
}

template<class T>
std::vector<ProductTerm<T> >
get_product_terms(const BasicOp<T> &op1, const BasicOp<T>& op2,
        const BasicOp<T> &op3, 
        const VectorState<T> &v, VectorState<T> &res,
        size_t m1, size_t m2, size_t m3,
        T coef = T(1), bool hc = false)
{
  std::vector<ProductTerm<T> > product_terms;
  Vector<QN> dqn(5);
  size_t _m1, _m2, _m3;
  dqn(m3) += op3.dqn;
  dqn(m2) += op2.dqn;
  dqn(m1) += op1.dqn;
  Vector<QN> dqn2(5);
  dqn2(m2) += op2.dqn;
  Vector<QN> dqn3(5);
  dqn3(m3) += op3.dqn;

  const BasicOp<T> *_op1;
  const BasicOp<T> *_op2;
  const BasicOp<T> *_op3;

  Vector<size_t> m(3), indx(3);
  Vector<const BasicOp<T>* > ops(3);
  ops(0) = &op1;
  ops(1) = &op2;
  ops(2) = &op3;
  m(0) = int(m1);
  m(1) = int(m2);
  m(2) = int(m3);
  indexx<size_t,Vector<size_t> >(3, m, indx);
  _m1 = (size_t)m(indx(0));
  _m2 = (size_t)m(indx(1));
  _m3 = (size_t)m(indx(2));
  _op1 = ops(indx(0));
  _op2 = ops(indx(1));
  _op3 = ops(indx(2));

  const SubMatrix<T> *_block1 = NULL;
  const SubMatrix<T> *_block2 = NULL;
  const SubMatrix<T> *_block3 = NULL;

  typename VectorState<T>::const_iterator siter;
  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++){
    StateSpace ss = *siter;

    _block3 = _op3->block(ss[_m3].qn());
    _block2 = _op2->block(ss[_m2].qn());
    _block1 = _op1->block(ss[_m1].qn());

    if(!_block1 || !_block2 || !_block3) continue;

    cstate_slice<T> v_slice = v(ss);
    state_slice<T> res_slice = res(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    if(res_slice.size() == 0) continue; 

    cstate_slice<T> v_slice_hc = v(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    state_slice<T> res_slice_hc = res(ss[1].qn(),ss[2].qn(),ss[3].qn(),ss[4].qn());
    if(hc && (v_slice_hc.size() == 0 || res_slice_hc.size() == 0)) continue; 

    int sign = 1;
    if(op3.fermion()){
      for(int ib = m3-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign();
    }
    if(op2.fermion()){
      for(int ib = m2-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn3[ib].fermion_sign();
    }
    if(op1.fermion()){
      for(int ib = m1-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn2[ib].fermion_sign()*dqn3[ib].fermion_sign();
    }

    ProductTerm<T> pterm;
    pterm.coef = T(sign)*coef;
    pterm.nops = 3;
    pterm.block1 = _block1;
    pterm.block2 = _block2;
    pterm.block3 = _block3;
    pterm.m1 = _m1;
    pterm.m2 = _m2;
    pterm.m3 = _m3;
    pterm.do_hc = hc;
    pterm.vspace = ss;
    pterm.res_space = res.get_qn_space(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2], ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    pterm.vspace_hc = v.get_qn_space(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2], ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    pterm.res_space_hc = ss;
    product_terms.push_back(pterm); 
  }
  return product_terms;
}


template<class T>
void
product_term(const ProductTerm<T> &pterm,
        const VectorState<T> &v, VectorState<T> &res, int mask_hc, 
        DMTKglobals<T> *globals = NULL)
{
  switch(pterm.nops){
    case 1:
      product_term1(pterm, v, res, mask_hc, globals);
      break;
    case 2:
      product_term2(pterm, v, res, mask_hc, globals);
      break;
    case 3:
      product_term3(pterm, v, res, mask_hc, globals);
      break;
    case 4:
      product_term4(pterm, v, res, mask_hc, globals);
      break;
  }
}

template<class T>
void
product_term1(const ProductTerm<T> &pterm,
        const VectorState<T> &v, VectorState<T> &res, int mask_hc, 
        DMTKglobals<T> *globals = NULL)
{
  Matrix<T> *_maux1, *_maux2;
  if(globals){
    _maux1 = &globals->m1;
    _maux2 = &globals->m2;
  } else {
    _maux1 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
    _maux2 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
  }
  Matrix<T> &maux1 = *_maux1;
  Matrix<T> &maux2 = *_maux2;

  StateSpace ss = pterm.vspace;
  cstate_slice<T> v_slice = v(pterm.vspace);
  state_slice<T> res_slice = res(pterm.res_space);

  const SubMatrix<T> &block1 = *pterm.block1;
  bool real_do_hc = pterm.do_hc;
  if(mask_hc & MASK_PRODUCT_HC) real_do_hc = true;
  char do_hc = (!real_do_hc) ? 'N' : 'C';
  T coef = pterm.coef;

  switch(pterm.m1){
      case BLOCK1:
        for(int i2 = 0; i2 < ss[2].dim(); i2++){
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,i3,slice(0,v_slice.size4(),1));
            gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,i3,slice(0,res_slice.size4(),1));
            maux2.reshape(subres.size2(),subres.size1());
            maux1 = subv;
            matrix_matrix_product(do_hc,'T',static_cast<const Matrix<T>&>(block1),maux1,maux2,coef);
            subres.transpose() += maux2.array();
          }
        }
        break;
      case BLOCK2:
        for(int i1 = 0; i1 < ss[1].dim(); i1++){
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),i3,slice(0,v_slice.size4(),1));
            gslice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),i3,slice(0,res_slice.size4(),1));
            maux1 = subv;
            maux2.reshape(subres.size2(),subres.size1());
            matrix_matrix_product(do_hc,'T',static_cast<const Matrix<T>&>(block1),maux1,maux2,coef);
            subres.transpose() += maux2.array();
          }
        }
        break;
      case BLOCK3:
        for(int i1 = 0; i1 < ss[1].dim(); i1++){
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),slice(0,v_slice.size3(),1),i4);
            gslice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),slice(0,res_slice.size3(),1),i4);
            maux1 = subv;
            maux2.reshape(subres.size1(),subres.size2());
            matrix_matrix_product(do_hc,'N',static_cast<const Matrix<T>&>(block1),maux1,maux2,coef);
            subres += maux2.array();
          }
        }
        break;
      case BLOCK4:
        for(int i2 = 0; i2 < ss[2].dim(); i2++){
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,i3,slice(0,v_slice.size4(),1));
            gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,i3,slice(0,res_slice.size4(),1));
            maux1 = subv;
            maux2.reshape(subres.size1(),subres.size2());
            matrix_matrix_product(do_hc,'N',static_cast<const Matrix<T>&>(block1),maux1,maux2,coef);
            subres += maux2.array();
          }
        }
        break;
  }
  if(!globals){
    delete(_maux1);
    delete(_maux2);
  }
}

template<class T>
void
product_term2(const ProductTerm<T> &pterm,
             const VectorState<T> &v, VectorState<T> &res, int mask_hc, 
             DMTKglobals<T> *globals = NULL, bool use_condensed = false)
{
  int _mask = mask(pterm.m1,pterm.m2);
  Vector<T> *_vaux1, *_vaux2;
  Matrix<T> *_maux1, *_maux2, *_maux3;
  if(globals){
    _vaux1 = &globals->v1;
    _vaux2 = &globals->v2;
    _maux1 = &globals->m1;
    _maux2 = &globals->m2;
    _maux3 = &globals->m3;
  } else {
    _vaux1 = new Vector<T>(MIN_VECTOR_SIZE);
    _vaux2 = new Vector<T>(MIN_VECTOR_SIZE);
    _maux1 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
    _maux2 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
    _maux3 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
  }
  Vector<T> &vaux1 = *_vaux1;
  Vector<T> &vaux2 = *_vaux2;
  Matrix<T> &maux1 = *_maux1;
  Matrix<T> &maux2 = *_maux2;
  Matrix<T> &maux3 = *_maux3;

  StateSpace ss = pterm.vspace;
  cstate_slice<T> v_slice = v(pterm.vspace);
  state_slice<T> res_slice = res(pterm.res_space);
  cstate_slice<T> v_slice_hc = v(pterm.vspace_hc);
  state_slice<T> res_slice_hc = res(pterm.res_space_hc);

  const SubMatrix<T> &block1 = *pterm.block1;
  const SubMatrix<T> &block2 = *pterm.block2;
  bool do_hc = pterm.do_hc;
  T coef = pterm.coef;

  switch(_mask){
      case (MASK_BLOCK1|MASK_BLOCK2):
        if(mask_hc & MASK_PRODUCT_DEFAULT){
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),
                                           slice(0,v_slice.size2(),1),i3,i4);
            gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),
                                              slice(0,res_slice.size2(),1),i3,i4);
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef,T(1));
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef,T(1));
              subres = maux2.array();
            }
          }
        }
        if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice_hc(slice(0,v_slice_hc.size1(),1),
                                              slice(0,v_slice_hc.size2(),1),i3,i4);
            gslice_iter<T> subres = res_slice_hc(slice(0,res_slice_hc.size1(),1),
                                                 slice(0,res_slice_hc.size2(),1),i3,i4);
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef,T(1),true);
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef,T(1),true);
              subres = maux2.array();
            }
          }
        }
        break;
      case (MASK_BLOCK1|MASK_BLOCK3):
        if(mask_hc & MASK_PRODUCT_DEFAULT){
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,
                                           slice(0,v_slice.size3(),1),i4);
            gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,
                                              slice(0,res_slice.size3(),1),i4);
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef,T(1));
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef,T(1));
              subres = maux2.array();
            }
          }
        }
        if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice_hc(slice(0,v_slice_hc.size1(),1),i2,
                                              slice(0,v_slice_hc.size3(),1),i4);
            gslice_iter<T> subres = res_slice_hc(slice(0,res_slice_hc.size1(),1),i2,
                                                 slice(0,res_slice_hc.size3(),1),i4);
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef,T(1),true);
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef,T(1),true);
              subres = maux2.array();
            }
          }
        }
        break;
     case (MASK_BLOCK1|MASK_BLOCK4):
        if(mask_hc & MASK_PRODUCT_DEFAULT){
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,i3,
                                           slice(0,v_slice.size4(),1));
            gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,i3,
                                              slice(0,res_slice.size4(),1));
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef,T(1));
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef,T(1));
              subres = maux2.array();
            }
          }
        }
        if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice_hc(slice(0,v_slice_hc.size1(),1),i2,i3,
                                              slice(0,v_slice_hc.size4(),1));
            gslice_iter<T> subres = res_slice_hc(slice(0,res_slice_hc.size1(),1),i2,i3,
                                                 slice(0,res_slice_hc.size4(),1));
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef,T(1),true);
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef,T(1),true);
              subres = maux2.array();
            }
          }
        }
        break;
     case (MASK_BLOCK2|MASK_BLOCK3):
        if(mask_hc & MASK_PRODUCT_DEFAULT){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),
                                           slice(0,v_slice.size3(),1),i4);
            gslice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),
                                              slice(0,res_slice.size3(),1),i4);
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef,T(1));
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef,T(1));
              subres = maux2.array();
            }
          }
        }
        if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice_hc(i1,slice(0,v_slice_hc.size2(),1),
                                           slice(0,v_slice_hc.size3(),1),i4);
            gslice_iter<T> subres = res_slice_hc(i1,slice(0,res_slice_hc.size2(),1),
                                              slice(0,res_slice_hc.size3(),1),i4);
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef,T(1),true);
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef,T(1),true);
              subres = maux2.array();
            }
          }
        }
        break;
     case (MASK_BLOCK2|MASK_BLOCK4):
        if(mask_hc & MASK_PRODUCT_DEFAULT){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),i3,
                                           slice(0,v_slice.size4(),1));
            gslice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),i3,
                                              slice(0,res_slice.size4(),1));
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef,T(1));
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef,T(1));
              subres = maux2.array();
            }
          }
        }
        if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice_hc(i1,slice(0,v_slice_hc.size2(),1),i3,
                                              slice(0,v_slice_hc.size4(),1));
            gslice_iter<T> subres = res_slice_hc(i1,slice(0,res_slice_hc.size2(),1),i3,
                                                 slice(0,res_slice_hc.size4(),1));
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef,T(1),true);
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef,T(1),true);
              subres = maux2.array();
            }
          }
        }
        break;
     case (MASK_BLOCK3|MASK_BLOCK4):
        if(mask_hc & MASK_PRODUCT_DEFAULT){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++){
            cgslice_iter<T> subv = v_slice(i1,i2,slice(0,v_slice.size3(),1),
                                           slice(0,v_slice.size4(),1));
            gslice_iter<T> subres = res_slice(i1,i2,slice(0,res_slice.size3(),1),
                                              slice(0,res_slice.size4(),1));
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef,T(1));
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef,T(1));
              subres = maux2.array();
            }
          }
        }
        if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++){
            cgslice_iter<T> subv = v_slice_hc(i1,i2,slice(0,v_slice_hc.size3(),1),
                                              slice(0,v_slice_hc.size4(),1));
            gslice_iter<T> subres = res_slice_hc(i1,i2,slice(0,res_slice_hc.size3(),1),
                                                 slice(0,res_slice_hc.size4(),1));
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef,T(1),true);
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef,T(1),true);
              subres = maux2.array();
            }
          }
        }
        break;
      case MASK_BLOCK1:
        if(mask_hc & MASK_PRODUCT_DEFAULT){
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,i3,i4);
            slice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,i3,i4);
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1));
            subres = vaux2.array();
          }
        }
        if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cslice_iter<T> subv = v_slice_hc(slice(0,v_slice_hc.size1(),1),i2,i3,i4);
            slice_iter<T> subres = res_slice_hc(slice(0,res_slice_hc.size1(),1),i2,i3,i4);
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1),true);
            subres = vaux2.array();
          }
        }
        break;
      case MASK_BLOCK2:
        if(mask_hc & MASK_PRODUCT_DEFAULT){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),i3,i4);
            slice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),i3,i4);
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1));
            subres = vaux2.array();
          }
        }
        if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cslice_iter<T> subv = v_slice_hc(i1,slice(0,v_slice_hc.size2(),1),i3,i4);
            slice_iter<T> subres = res_slice_hc(i1,slice(0,res_slice_hc.size2(),1),i3,i4);
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1),true);
            subres = vaux2.array();
          }
        }
        break;
      case MASK_BLOCK3:
        if(mask_hc & MASK_PRODUCT_DEFAULT){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cslice_iter<T> subv = v_slice(i1,i2,slice(0,v_slice.size3(),1),i4);
            slice_iter<T> subres = res_slice(i1,i2,slice(0,res_slice.size3(),1),i4);
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1));
            subres = vaux2.array();
          }
        }
        if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cslice_iter<T> subv = v_slice_hc(i1,i2,slice(0,v_slice_hc.size3(),1),i4);
            slice_iter<T> subres = res_slice_hc(i1,i2,slice(0,res_slice_hc.size3(),1),i4);
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1),true);
            subres = vaux2.array();
          }
        }
        break;
     case MASK_BLOCK4:
        if(mask_hc & MASK_PRODUCT_DEFAULT){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cslice_iter<T> subv = v_slice(i1,i2,i3,slice(0,v_slice.size4(),1));
            slice_iter<T> subres = res_slice(i1,i2,i3,slice(0,res_slice.size4(),1));
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1));
            subres = vaux2.array();
          }
        }
        if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cslice_iter<T> subv = v_slice_hc(i1,i2,i3,slice(0,v_slice_hc.size4(),1));
            slice_iter<T> subres = res_slice_hc(i1,i2,i3,slice(0,res_slice_hc.size4(),1));
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1),true);
            subres = vaux2.array();
          }
        }
        break;
  }
  if(!globals){
    delete(_vaux1);
    delete(_vaux2);
    delete(_maux1);
    delete(_maux2);
    delete(_maux3);
  }
}

template<class T>
void
product_term4(const ProductTerm<T> &pterm,
             const VectorState<T> &v, VectorState<T> &res, int mask_hc, 
             DMTKglobals<T> *globals = NULL)
{
  StateSpace ss = pterm.vspace;
  cstate_slice<T> v_slice = v(pterm.vspace);
  state_slice<T> res_slice = res(pterm.res_space);
  cstate_slice<T> v_slice_hc = v(pterm.vspace_hc);
  state_slice<T> res_slice_hc = res(pterm.res_space_hc);

  const SubMatrix<T> &block1 = *pterm.block1;
  const SubMatrix<T> &block2 = *pterm.block2;
  const SubMatrix<T> &block3 = *pterm.block3;
  const SubMatrix<T> &block4 = *pterm.block4;
  bool do_hc = pterm.do_hc;
  T coef = pterm.coef;

  if(mask_hc & MASK_PRODUCT_DEFAULT){
    for(int i1 = 0; i1 < ss[1].dim(); i1++)
      for(int i2 = 0; i2 < ss[2].dim(); i2++)
        for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++)
            for(int j1 = 0; j1 < res_slice.size1(); j1++)
              for(int j2 = 0; j2 < res_slice.size2(); j2++)
                for(int j3 = 0; j3 < res_slice.size3(); j3++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice(j1,j2,j3,j4) += coef*v_slice(i1,i2,i3,i4)*block1(i1,j1)*block2(i2,j2)*block3(i3,j3)*block4(i4,j4);
  }
  if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
    for(int i1 = 0; i1 < ss[1].dim(); i1++)
      for(int i2 = 0; i2 < ss[2].dim(); i2++)
        for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++)
            for(int j1 = 0; j1 < res_slice.size1(); j1++)
              for(int j2 = 0; j2 < res_slice.size2(); j2++)
                for(int j3 = 0; j3 < res_slice.size3(); j3++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice_hc(i1,i2,i3,i4) += dmtk::conj(coef)*v_slice_hc(j1,j2,j3,j4)*dmtk::conj(block1(i1,j1)*block2(i2,j2)*block3(i3,j3)*block4(i4,j4));
  }
}

template<class T>
void
product_term3(const ProductTerm<T> &pterm,
             const VectorState<T> &v, VectorState<T> &res, int mask_hc, 
             DMTKglobals<T> *globals = NULL)
{
  StateSpace ss = pterm.vspace;
  cstate_slice<T> v_slice = v(pterm.vspace);
  state_slice<T> res_slice = res(pterm.res_space);
  cstate_slice<T> v_slice_hc = v(pterm.vspace_hc);
  state_slice<T> res_slice_hc = res(pterm.res_space_hc);

  const SubMatrix<T> &block1 = *pterm.block1;
  const SubMatrix<T> &block2 = *pterm.block2;
  const SubMatrix<T> &block3 = *pterm.block3;
  bool do_hc = pterm.do_hc;
  T coef = pterm.coef;
  int _m1 = pterm.m1;
  int _m2 = pterm.m2;
  int _m3 = pterm.m3;

  if(_m1 == BLOCK1 && _m2 == BLOCK2 && _m3 == BLOCK3){
    if(mask_hc & MASK_PRODUCT_DEFAULT){
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j1 = 0; j1 < res_slice.size1(); j1++)
                for(int j2 = 0; j2 < res_slice.size2(); j2++)
                  for(int j3 = 0; j3 < res_slice.size3(); j3++)
                    res_slice(j1,j2,j3,i4) += coef*v_slice(i1,i2,i3,i4)*block1(i1,j1)*block2(i2,j2)*block3(i3,j3);
    }
    if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j1 = 0; j1 < res_slice.size1(); j1++)
                for(int j2 = 0; j2 < res_slice.size2(); j2++)
                  for(int j3 = 0; j3 < res_slice.size3(); j3++)
                      res_slice_hc(i1,i2,i3,i4) += dmtk::conj(coef)*v_slice_hc(j1,j2,j3,i4)*dmtk::conj(block1(i1,j1)*block2(i2,j2)*block3(i3,j3));
    }
  } else if(_m1 == BLOCK1 && _m2 == BLOCK2 && _m3 == BLOCK4){
    if(mask_hc & MASK_PRODUCT_DEFAULT){
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j1 = 0; j1 < res_slice.size1(); j1++)
                for(int j2 = 0; j2 < res_slice.size2(); j2++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice(j1,j2,i3,j4) += coef*v_slice(i1,i2,i3,i4)*block1(i1,j1)*block2(i2,j2)*block3(i4,j4);
    }
    if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j1 = 0; j1 < res_slice.size1(); j1++)
                for(int j2 = 0; j2 < res_slice.size2(); j2++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice_hc(i1,i2,i3,i4) += dmtk::conj(coef)*v_slice_hc(j1,j2,i3,j4)*dmtk::conj(block1(i1,j1)*block2(i2,j2)*block3(i4,j4));
    }
 } else if(_m1 == BLOCK1 && _m2 == BLOCK3 && _m3 == BLOCK4){
    if(mask_hc & MASK_PRODUCT_DEFAULT){
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j1 = 0; j1 < res_slice.size1(); j1++)
                for(int j3 = 0; j3 < res_slice.size3(); j3++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice(j1,i2,j3,j4) += coef*v_slice(i1,i2,i3,i4)*block1(i1,j1)*block2(i3,j3)*block3(i4,j4);
    }
    if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j1 = 0; j1 < res_slice.size1(); j1++)
                for(int j3 = 0; j3 < res_slice.size3(); j3++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice_hc(i1,i2,i3,i4) += dmtk::conj(coef)*v_slice_hc(j1,i2,j3,j4)*dmtk::conj(block1(i1,j1)*block2(i3,j3)*block3(i4,j4));
    }
  } else if(_m1 == BLOCK2 && _m2 == BLOCK3 && _m3 == BLOCK4){
    if(mask_hc & MASK_PRODUCT_DEFAULT){
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j2 = 0; j2 < res_slice.size2(); j2++)
                for(int j3 = 0; j3 < res_slice.size3(); j3++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice(i1,j2,j3,j4) += coef*v_slice(i1,i2,i3,i4)*block1(i2,j2)*block2(i3,j3)*block3(i4,j4);
    }
    if(do_hc && (mask_hc & MASK_PRODUCT_HC)){
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j2 = 0; j2 < res_slice.size2(); j2++)
                for(int j3 = 0; j3 < res_slice.size3(); j3++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice_hc(i1,i2,i3,i4) += dmtk::conj(coef)*v_slice_hc(i1,j2,j3,j4)*dmtk::conj(block1(i2,j2)*block2(i3,j3)*block3(i4,j4));
    }
  }
}

} // namespace dmtk
#endif // __DMTK_PRODUCT_H__
