/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2008 by Salvatore R. Manmana <Salva@theo3.physik.uni-stuttgart.de>,
*                            Reinhard M. Noack <Reinhard.Noack@physik.uni-marburg.de>,
*                            Ian McCulloch <ianmcc@physics.uq.edu.au>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#include "dmrg3.h"

// Single-particle DMRG main function

// Output information is controlled by the OUTPUT_LEVEL parameter.
// Output levels are cumulative
// level 0 : no messages
// level 1 : show energy at end of each sweep
// level 2 : show energy at each iteration
// level 3 : show wavefunction at end of each sweep
// level 4 : show detailed debug information

int main() 
{
  try
  {
    // Construct a System from the input parameters
    std::cout << "ALPS noninteractiong DMRG application version 1.0\n"
      << "  available from http://alps.comp-phys.org/\n"
      << "  copyright (c) 2003-2004 by Salvatore R. Manmana <Salva@theo3.physik.uni-stuttgart.de>\n"
      << "                             Reinhard M. Noack <Reinhard.Noack@physik.uni-marburg.de>\n"
      << "                             Ian McCulloch <ianmcc@physik.rwth-aachen.de>\n"
      << "  for details see Chapter 3 of \"Density Matrix Renormalization - A New Numerical\n"
      << "  Method in Physics\", edited by I. Peschel, X. Wang, M. Kaulke, and K. Hallberg\n"
      << "  (Springer, Berlin, 1999)\n"
      << "\n"
      << "  WARNING: This program exists for demonstration purposes only, and may not converge\n"
      << "  well for complicated potentials.  If the results matter, use one-body exact\n"
      << "  diagonalization instead.\n";

    alps::Parameters parameters;
    std::cin >> parameters;

    System DMRGSystem(parameters);

    int Prec = parameters.value_or_default("PRECISION", 10);
    std::cout.precision(Prec);

    int nsweeps = parameters["SWEEPS"];

    DMRGSystem.show_info(std::cout);
    std::cout << "Number of sweeps = " << nsweeps << std::endl;

    int output_level = parameters.value_or_default("OUTPUT_LEVEL", 1);

    std::string OutFile = parameters.value_or_default("WAVEFUNCTION_FILE", std::string());

    // The DMRG linearizes the lattice onto a chain; this mapping is given by
    // System::site_at(i) for the i'th site of the DMRG lattice.  Use this to construct
    // an array of single-site blocks.
    std::vector<Block> site_block;
    for (int i = 0; i < DMRGSystem.num_sites(); ++i)
    {
      site_block.push_back(Block(DMRGSystem, DMRGSystem.site_at(i)));
    }

    // ground state wave function
    Wavefunction psi;
    // GS energy
    double energy = 0.0;                    

    std::stack<Block> LeftBlocks, RightBlocks;
  
    // Warmup sweep, using the infinite-system algorithm

    LeftBlocks.push(site_block[0]);

    int Lmax = DMRGSystem.num_sites() - 3;
    int Lmin = 1;
  
    if (output_level >= 2) std::cout << "Starting warmup sweep:" << std::endl;

    // start iterations, end when system size reached Lmax (s. above):
    for(int i = 1; i < Lmax; i++)
    {
      // compose the total system 'superblock'
      // in the first step, this is the simple Hamiltonian
      // for a four site system:
      Superblock S(LeftBlocks.top(), site_block[i], site_block[i+1], site_block[i+2], output_level >= 4);
      
      // ground state vector and energy for this iteration step
      boost::tie(energy, psi) = S.GetGroundState(); 

      // output of wave function and energy 
      if (output_level >= 2) std::cout << "i = " << i << ", energy = " << energy << std::endl;
      if (output_level >= 4) psi.print_debug(std::cout);

      // reduced density matrix of left system block:
      DensityMatrix rho(psi, DensityMatrix::Left);
      
      // new block Hamiltonian in 'DMRG basis'
      // for given values of t, t_prime
      LeftBlocks.push(Block(LeftBlocks.top(), site_block[i], rho, output_level >= 4));
    }
  
    // Finite System Sweeps

    TotalWavefunction Psi;

    // initial right block  
    RightBlocks.push(site_block[DMRGSystem.num_sites()-1]);

    if (output_level >= 2) std::cout << "Starting finite-size sweeping:" << std::endl;

    for (int swp = 1; swp <= nsweeps; swp++)
    {
      // Right to Left:
      if (output_level >= 2) std::cout << "Right to left sweep:" << std::endl;
      for(int i = Lmax; i > Lmin; i--)
      {
        Superblock S(LeftBlocks.top(), site_block[i], site_block[i+1], RightBlocks.top());

        // ground state vector and energy for this iteration step
        boost::tie(energy, psi) = S.GetGroundState(); 

        // output of wave function and energy 
        if (output_level >= 2) std::cout << "i = " << i << ", energy = " << energy << std::endl;
        if (output_level >= 4) psi.print_debug(std::cout);

        DensityMatrix rho(psi,DensityMatrix::Right);

        LeftBlocks.pop();
        RightBlocks.push(Block(RightBlocks.top(), site_block[i+1], rho, output_level >= 4));
      } // right to left sweep

      // left to right      
      if (output_level >= 2) std::cout << "Left to right sweep:" << std::endl;
      for (int i = 1; i < Lmax; i++)
      {
        Superblock S(LeftBlocks.top(), site_block[i], site_block[i+1], RightBlocks.top(), output_level >= 4);

        // ground state vector and energy for this iteration step
        boost::tie(energy, psi) = S.GetGroundState(); 

        // output of wave function and energy 
        if (output_level >= 2) std::cout << "i = " << i << ", energy = " << energy << std::endl;
        if (output_level >= 4) psi.print_debug(std::cout);

        // Determine the total wavefunction, if necessary
        if (i == Lmax-1 && (swp == nsweeps || output_level >= 3))
          Psi = TotalWavefunction(LeftBlocks.top(), site_block[i], site_block[i+1], RightBlocks.top(), psi);

        DensityMatrix rho(psi, DensityMatrix::Left);
      
        RightBlocks.pop();
        LeftBlocks.push(Block(LeftBlocks.top(), site_block[i], rho, output_level >= 4));
      } // left to right sweep

      if (output_level >= 2) std::cout << '\n';
      if (output_level >= 1)
        std::cout << "***End of sweep " << swp << ", energy = " << energy << std::endl;

      if (output_level >= 3)
        std::cout << "***End of sweep " << swp << ", wavefuntion = \n" << Psi << std::endl;
    }

    if (output_level >= 1) std::cout << '\n';
    std::cout << "***final result:  energy = " << energy << std::endl;
    std::cout << "***final wavefunction:\n" << Psi;

    // Write the wavefunction to a file if requested
    if (!OutFile.empty())
    {
      std::ofstream Out(OutFile.c_str());
      Out.precision(Prec);
      Out << Psi;
      Out.close();
    }

    return 0;

  }
  catch (std::exception& e)
  {
    std::cout << "EXCEPTION: " << e.what() << std::endl;
  }
  catch (...)
  {
    std::cout << "UNKNOWN EXCEPTION!\n";
  }

  return 1;
  
} // main end


