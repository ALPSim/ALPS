/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2009 by Matthias Troyer <troyer@comp-phys.org>,
*                       Synge Todo <wistaria@comp-phys.org>
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

/* $Id$ */

#include "factory.h"
#include "abstractspinsim.h"

#include "spinsim.h"

#include "ising.h"
#include "on.h"
#include "potts.h"

void SpinFactory::print_copyright(std::ostream& out) const
{
  out << "Generic classical Monte Carlo program using local or cluster updates\n"
      << "  available from http://alps.comp-phys.org/\n"
      << "  copyright(c) 1999-2007 by Matthias Troyer <troyer@comp-phys.org>\n"
      << "                            Mathias Koerner <mkoerner@comp-phys.org>\n"
      << " for details see the publication:\n"
      << " A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).\n\n";
}


alps::scheduler::MCSimulation* SpinFactory::make_task(
        const alps::ProcessList& w,
        const boost::filesystem::path& fn) const
{
  return new alps::scheduler::MCSimulation(w,fn);
}

/** 
 * Returns the number of white-space separated elements in the string.
 *
 * \param str the string to be analyzed.
 */
int SpinFactory::countElements(const std::string& str) const
{
  int counter = 0;
  int pos = 0;
  if (str[0] == '\0') return 0;
  while (true) {
    while (str[pos] == ' ') pos++;
    if (str[pos] == '\0') return counter;
    counter++;
    while ((str[pos] !='\0') && (str[pos] != ' ')) pos++;
  }
}

/**
 * Returns the number of elements in the string with the largest number of elements.
 */
int SpinFactory::findDominantMatrixString(const alps::Parameters& parms) const
{
  int count = 0;
  /* should be done more elegantly ... */
  if (parms.defined("J")) 
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, countElements(parms["J"]));
  if (parms.defined("J'"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, countElements(parms["J'"]));
  if (parms.defined("J0"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, countElements(parms["J0"]));
  if (parms.defined("J1"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, countElements(parms["J1"]));
  if (parms.defined("J2"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, countElements(parms["J2"]));
  if (parms.defined("J3"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, countElements(parms["J3"]));
  if (parms.defined("D"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, countElements(parms["D"]));
  if (parms.defined("D'"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, countElements(parms["D'"]));
  if (parms.defined("D0"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, countElements(parms["D0"]));
  if (parms.defined("D1"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, countElements(parms["D1"]));
  if (parms.defined("D2"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, countElements(parms["D2"]));
  if (parms.defined("D3"))
    count = std::max BOOST_PREVENT_MACRO_SUBSTITUTION (count, countElements(parms["D3"]));

  return count;
}

/**
 * Outputs an error message and aborts the computation if an invalid number of
 * elements is provided in the input for a matrix parameter.
 */
void SpinFactory::produceError(const alps::Parameters& parms) const
{
  std::cerr << "Invalid parameters given:\n"
            << "Cannot produce model " << parms["MODEL"] << " with strings "
            << parms["D"] << " for D and " << parms["J"] 
            << " for J.\nAborting\n";
  boost::throw_exception(std::runtime_error("invalid parameters"));
  std::exit(-5);
}

/**
 * Creates a worker, according to the parameters. The type of the produced 
 * worker depends on the model given in the parameter file and on the number of
 * input elements for the matrix parameters. The chosen matrix type is the most
 * general one (over all matrix parameters) of all most specialized ones that
 *  are compatible with the user input.
 */
alps::scheduler::Worker* SpinFactory::make_worker(
        const alps::ProcessList& where,
        const alps::Parameters& parms,int node) const
{
  int maxElemCount = findDominantMatrixString(parms);
  if (parms["MODEL"]=="Ising")
    if (maxElemCount == 1) 
      return new SpinSim<IsingMoment,MIdMatrix<double,1> >(where,parms,node);
    else {
      produceError(parms);
      return 0;
    }
  else if (parms["MODEL"]=="O(4)")
    return new SpinSim<ONMoment<4>,MIdMatrix<double,4> >(where,parms,node);
  else if (parms["MODEL"]=="Heisenberg")
    switch (maxElemCount) {
      case 1:
        return new SpinSim< ONMoment<3>, MIdMatrix<double,3> >(where, parms, node);
        break;
      case 3:
        return new SpinSim< ONMoment<3>, DiagMatrix<double,3> >(where, parms, node);
        break;
      case 6:
        return new SpinSim< ONMoment<3>, SquareMatrix<double,3> >(where, parms, node);
        break;
      case 9:
        return new SpinSim< ONMoment<3>, SquareMatrix<double,3> >(where, parms, node);
        break;
      default:
        produceError(parms);
        return 0;
    }
  else if (parms["MODEL"]=="XY")
    switch (maxElemCount) {
      case 1:
        return new SpinSim< XYMoment, MIdMatrix<double,2> >(where, parms, node);
        break;
      case 2:
        return new SpinSim< XYMoment, DiagMatrix<double,2> >(where, parms, node);
        break;
      case 3:
        return new SpinSim< XYMoment, SquareMatrix<double,2> >(where, parms, node);
        break;
      case 4:
        return new SpinSim< XYMoment, SquareMatrix<double,2> >(where,parms,node);
        break;
      default:
        produceError(parms);
        return 0;
    }
  
  else if (parms["MODEL"]=="Potts")
    switch (int(parms["q"])) {
    case 3:
      return new SpinSim< PottsMoment<3>, MIdMatrix<double,2> >(where,parms,node);
    case 4:
      return new SpinSim< PottsMoment<4>, MIdMatrix<double,2> >(where,parms,node);
    case 10:
      return new SpinSim< PottsMoment<10>, MIdMatrix<double,2> >(where,parms,node);
    default:
      boost::throw_exception(std::runtime_error(std::string(parms["q"]) 
            + "-state Potts model not implemented in file factory.C"));
    }
  else 
    boost::throw_exception(std::runtime_error(std::string(parms["MODEL"])
              + " model not implemented in file factory.C"));
  return 0;
}
