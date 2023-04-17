// clang-format off
//
// ==============================================================================
//       MM\ MMMMMMMM\ MMMMMMMM\ MM\  MMMMMM\  MM\ MMMMMMMM\ MMMMMMMM\ MM\      =
//      MM  |MM  _____|\__MM  __|MM |MM  __MM\ MM |\__MM  __|MM  _____|\MM\     =
//     MM  / MM |         MM |   MM |MM /  MM |MM |   MM |   MM |       \MM\    =
//    MM  /  MMMMM\       MM |   MM |MM |  MM |MM |   MM |   MMMMM\      \MM\   =
//    \MM<   MM  __|      MM |   MM |MM |  MM |MM |   MM |   MM  __|     MM  |  =
//     \MM\  MM |         MM |   MM |MM MM\MM |MM |   MM |   MM |       MM  /   =
//      \MM\ MMMMMMMM\    MM |   MM |\MMMMMM / MM |   MM |   MMMMMMMM\ MM  /    =
//       \__|\________|   \__|   \__| \___MMM\ \__|   \__|   \________|\__/     =
//                                        \___|                                 =
// ==============================================================================
// ==============================================================================
// ========== Energy    Transfer       Quantum        Time   Evolution ==========
// ================================ Version 1.0.0 ===============================
// ================================== April 2023 ================================
// ==============================================================================
// ===================== J.-R. Li, Y. Zhai, Z. Qu, and H. Li ====================
// ==============================================================================
// This program performs the time evolution based on quantum mechanics.
// The detailed theory can be found in the manual.
// The advantage of this program is that it is free from complex computations,
// and utilize Chebyshev propagator to get the population evolution.
// 
// This program is released under MIT licence.
// 
// Copyright 2023 the authors of ETQTE
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the “Software”), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// clang-format on

#define TOML_ENABLE_FLOAT16 0
#define TOML_ENABLE_UNRELEASED_FEATURES 0
#define TOML_ENABLE_WINDOWS_COMPAT 0
#define TOML_ENABLE_FORMATTERS 0
#define FMT_HEADER_ONLY 1
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <fmt/compile.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <toml.hh>
#include <utility>
#include <vector>

namespace etqte
{
class isosys
{
public:
  // General settings
  size_t nstate;
  std::vector<std::string> statenames;
  std::string matrixfname;

  // Workspace
  Eigen::MatrixXd H;
  Eigen::VectorXd E;
  Eigen::MatrixXd psi;
  size_t istate;
  size_t nstep;
  double tstep;
  std::vector<Eigen::VectorXd> rho;
  Eigen::VectorXd rho0;

public:
  // Constructor
  isosys(size_t nstate,
         std::vector<std::string> statenames,
         std::string matrixfname)
    : nstate(nstate)
    , statenames(statenames)
    , matrixfname(matrixfname)
    , H()
    , E()
    , psi()
    , istate(0)
    , nstep(0)
    , tstep(0.0)
    , rho()
    , rho0(nstate)
  {
    readhamiltonian();
  }

  // Constructor
  isosys(size_t nstate,
         std::vector<size_t> states_selected,
         size_t nstatep,
         std::vector<std::string> statenamesp,
         std::string matrixfname)
    : nstate(nstate)
    , statenames()
    , matrixfname(matrixfname)
    , H()
    , E()
    , psi()
    , istate(0)
    , nstep(0)
    , tstep(0.0)
    , rho()
    , rho0(nstate)
  {
    for (size_t i = 0; i != nstate; ++i)
    {
      statenames.push_back(statenamesp[i]);
    }
    readhamiltonian(nstatep, states_selected);
  }

  // This function read in the elements of effective Hamiltonian
  // Here I follow Jia-Rui's convention.
  int readhamiltonian()
  {
    std::ifstream io(matrixfname);
    int tmp;
    for (size_t i = 0; i != nstate; ++i)
    {
      for (size_t j = 0; j != nstate; ++j)
      {
        io >> tmp >> tmp >> H(i, j);
      }
    }
    io.close();

    // Apparently the matrices are in eV, Errrrr...
    // Turn hamiltonian from eV to Hartree. Idiot.
    H /= 27.21138602;
    return 0;
  }

  // This function read in the elements of effective Hamiltonian
  // Here I follow Jia-Rui's convention.
  int readhamiltonian(size_t nstatep, std::vector<size_t> states_selected)
  {
    std::ifstream io(matrixfname);
    int tmp;
    Eigen::MatrixXd H(nstatep, nstatep);
    for (size_t i = 0; i != nstatep; ++i)
    {
      for (size_t j = 0; j != nstatep; ++j)
      {
        io >> tmp >> tmp >> H(i, j);
      }
    }
    io.close();

    this->H = H(states_selected, states_selected);
    // Apparently the matrices are in eV, Errrrr...
    // Turn hamiltonian from eV to Hartree. Idiot.
    this->H /= 27.21138602;
    return 0;
  }

  static int writetimeseries(std::string fname,
                             size_t nstate,
                             const std::vector<std::string>& statenames,
                             size_t nstep,
                             double tstep,
                             const std::vector<Eigen::VectorXd>& instval,
                             const Eigen::VectorXd& infval)
  {
    FILE* output = fopen(fname.c_str(), "w");
    fmt::print(output, "#             time");
    for (size_t state = 0; state != nstate; ++state)
    {
      fmt::print(output, FMT_COMPILE("{:>13s}"), statenames[state]);
    }
    fmt::print(output, "\n");
    for (size_t step = 0; step != nstep; ++step)
    {
      fmt::print(output, FMT_COMPILE("{:18.5f}"), step * tstep);
      for (size_t state = 0; state != nstate; ++state)
      {
        fmt::print(output, FMT_COMPILE("{:13.5f}"), instval[step][state]);
      }
      fmt::print(output, "\n");
    }
    // We use a long enough time to present an infinite long time
    fmt::print(output, "\n");
    fmt::print(output, FMT_COMPILE("{:18.5f}"), nstep * tstep * 1.2);
    for (size_t state = 0; state != nstate; ++state)
    {
      fmt::print(output, FMT_COMPILE("{:13.5f}"), infval[state]);
    }
    fmt::print(output, "\n");
    fmt::print(output, FMT_COMPILE("{:18.5f}"), nstep * tstep * 2);
    for (size_t state = 0; state != nstate; ++state)
    {
      fmt::print(output, FMT_COMPILE("{:13.5f}"), infval[state]);
    }
    fmt::print(output, "\n");

    fclose(output);

    return 0;
  }

  int writerhos()
  {
    writetimeseries(matrixfname + std::string(".etcc"),
                    nstate,
                    statenames,
                    nstep,
                    tstep,
                    rho,
                    rho0);
    return 0;
  }

  int eigenstates()
  {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
    E = es.eigenvalues();
    psi = es.eigenvectors();
    return 0;
  }

  // The quantum time evolution function
  int timeevolution(size_t istate, size_t nstep, double tstep)
  {
    this->istate = istate;
    this->nstep = nstep;
    this->tstep = tstep;
    const size_t upsize = nstate * (nstate - 1) / 2;

    // Actually +2 is good enough, but I think +5 is no harm.
    rho.clear();
    rho.reserve(nstep + 5);
    rho0 = Eigen::VectorXd::Zero(nstate);

    Eigen::MatrixXd param(nstate, upsize);
    Eigen::VectorXd DeltaE(upsize);

    for (size_t fstate = 0; fstate != nstate; ++fstate)
    {
      for (size_t i = 0; i != nstate; ++i)
      {
        rho0(fstate) +=
          std::pow(psi(istate, i), 2) * std::pow(psi(fstate, i), 2);
      }
      size_t iii = 0;
      for (size_t i = 0; i != nstate; ++i)
      {
        for (size_t j = i + 1; j != nstate; ++j)
        {
          param(fstate, iii) =
            psi(istate, i) * psi(istate, j) * psi(fstate, i) * psi(fstate, j);
          ++iii;
        }
      }
    }
    param *= 2.0;
    size_t iii = 0;
    for (size_t i = 0; i != nstate; ++i)
    {
      for (size_t j = i + 1; j != nstate; ++j)
      {
        DeltaE(iii) = E(j) - E(i);
        ++iii;
      }
    }

    // N.B.
    // Here I use the Chebyshev propagator, because the cosine function is
    // extremely slow.
    // Ref: Chen & Guo, Comput. Phys. Comm. 119, 19-31 (1999)
    Eigen::VectorXd omega =
      (tstep * DeltaE).unaryExpr<double (*)(double)>(&std::cos);

    // The follwoing is for k = 0...
    // Eigen::VectorXd chebyshev_k(upsize);
    // Eigen::VectorXd chebyshev_k_1(upsize);
    Eigen::VectorXd chebyshev_k = Eigen::VectorXd::Ones(upsize);
    Eigen::VectorXd chebyshev_k_1 = omega;
    auto* chb_k = &chebyshev_k;
    auto* chb_k_1 = &chebyshev_k_1;
    for (size_t step = 0; step != nstep; ++step)
    {
      Eigen::VectorXd rho = param * (*chb_k) + rho0;
      this->rho.push_back(std::move(rho));
      *chb_k_1 = (2 * omega.array() * chb_k->array()).matrix() - (*chb_k_1);
      std::swap(chb_k, chb_k_1);
    }

    return 0;
  }

  isosys subsystem(std::vector<size_t>& states_selected)
  {
    isosys subsys = *this;
    subsys.nstate = states_selected.size();
    subsys.H = this->H(states_selected, states_selected);
    subsys.statenames.clear();
    for (auto&& i : states_selected)
    {
      subsys.statenames.push_back(this->statenames[i]);
    }
    return subsys;
  }
};

class ensemble
{
private:
  // General settings
  size_t nstate, nsys, istate, nstep, nstatep;
  std::string prefix, suffix;
  std::vector<std::string> statenames;
  std::vector<std::string> statenamesp;
  std::vector<size_t> states_selected;
  double tstep;

  // Workspace
  std::vector<isosys> systems;
  std::vector<Eigen::VectorXd> averrho;
  Eigen::VectorXd averrho0;

public:
  int writerhos()
  {
    for (auto&& sys : systems)
    {
      sys.writerhos();
    }
    return 0;
  }

  // Read the input file, a sample of which is provided.
  int readsettings(const std::string& fname)
  {
    auto settings = toml::parse_file(fname);
    nstate = settings["etqtesettings"]["number_of_states_selected"].value_or(0);
    nstatep =
      settings["etqtesettings"]["number_of_states_presented"].value_or(0);
    prefix = settings["etqtesettings"]["prefix_of_matrix_files"].value_or("");
    suffix = settings["etqtesettings"]["suffix_of_matrix_files"].value_or("");
    nsys = settings["etqtesettings"]["number_of_systems"].value_or(0);
    nstep = settings["etqtesettings"]["number_of_steps"].value_or(0);
    tstep = settings["etqtesettings"]["step_length"].value_or(0.0);
    istate = settings["etqtesettings"]["initial_state"].value_or(0) - 1;
    statenamesp.clear();
    std::cout << "Ln. " << __LINE__ << std::endl;
    for (auto&& i :
         *settings["etqtesettings"]["state_names_presented"].as_array())
    {
      statenamesp.push_back(i.as_string()->value_or(""));
    }
    statenames.clear();
    for (auto&& i :
         *settings["etqtesettings"]["state_names_selected"].as_array())
    {
      statenames.push_back(i.as_string()->value_or(""));
    }

    std::cout << "Ln. " << __LINE__ << std::endl;
    // Health check
    bool health = true;
    do
    {
      if (statenames.size() != nstate)
      {
        std::cout << "*** The size of state_names does not consist with nstate."
                  << std::endl;
        health = false;
        break;
      }

      if (statenamesp.size() != nstatep)
      {
        std::cout << "*** The size of state_names does not consist with nstate."
                  << std::endl;
        health = false;
        break;
      }
    } while (false);

    if (!health)
    {
      exit(-1);
    }

    std::cout << "Ln. " << __LINE__ << std::endl;
    std::map<std::string, size_t> statemap;
    for (size_t i = 0; i != nstatep; ++i)
    {
      statemap[statenamesp[i]] = i;
    }
    std::cout << "Ln. " << __LINE__ << std::endl;

    for (auto&& i : statenames)
    {
      states_selected.push_back(statemap[i]);
    }
    std::cout << "Ln. " << __LINE__ << std::endl;

    return 0;
  }

  // This function read in the elements of effective Hamiltonian
  // Here I follow Jia-Rui's convention.
  int readhamiltonian()
  {
    if (nstate == 0)
    {
      for (size_t sys = 0; sys != nsys; ++sys)
      {
        std::string fname = fmt::format("{}{}{}", prefix, sys + 1, suffix);
        systems.push_back(isosys(nstatep, statenamesp, fname));
      }
    }
    else
    {
      for (size_t sys = 0; sys != nsys; ++sys)
      {
        std::string fname = fmt::format("{}{}{}", prefix, sys + 1, suffix);
        systems.push_back(
          isosys(nstate, states_selected, nstatep, statenamesp, fname));
      }
    }
    return 0;
  }

  // The quantum time evolution function
  int timeevolution()
  {
    for (auto&& sys : systems)
    {
      sys.eigenstates();
      sys.timeevolution(istate, nstep, tstep);
    }
    return 0;
  }

  ensemble(const std::string& input)
  {
    readsettings(input);
    return;
  }

  int driver2023()
  {
    std::cout << "*** Read Hamiltonian matrices in...\n";
    readhamiltonian();
    std::cout << "*** Done.\n" << std::flush;

    {
      std::cout << "*** Start time evolution...\n";
      auto&& t0 = std::chrono::high_resolution_clock::now();
      timeevolution();
      auto&& t1 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double>&& difft = t1 - t0;
      std::cout << "*** End time evolution.\n"
                << "*** Job is done in " << difft.count() << " sec."
                << std::endl;
    }
    {
      auto&& t0 = std::chrono::high_resolution_clock::now();
      writerhos();
      auto&& t1 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double>&& difft = t1 - t0;
      std::cout << "*** Time evolution of each sampled system has been written "
                   "to .etcc files.\n"
                << "*** Job is done in " << difft.count() << " sec."
                << std::endl;
    }

    std::cout << "*** Start ensemble average..." << std::endl;
    averrho.resize(nstep);
    for (size_t step = 0; step != nstep; ++step)
    {
      averrho[step] = Eigen::VectorXd::Zero(nstate);
      for (size_t sys = 0; sys != nsys; ++sys)
      {
        averrho[step] += systems[sys].rho[step];
      }
      averrho[step] /= nsys;
    }
    averrho0 = Eigen::VectorXd::Zero(nstate);
    for (size_t sys = 0; sys != nsys; ++sys)
    {
      averrho0 += systems[sys].rho0;
    }
    averrho0 /= nsys;

    isosys::writetimeseries(fmt::format("{}{}{}.etcc", prefix, 0, suffix),
                            nstate,
                            statenames,
                            nstep,
                            tstep,
                            averrho,
                            averrho0);
    std::cout << "*** Time evolution of the whole ensemble has been written to "
                 "the zeroth .etcc files."
              << std::endl;

    return 0;
  }
};

static int
printbanner(std::ostream& io)
{
  io << R"foo(
 ===============================================================================
 =      MM\ MMMMMMMM\ MMMMMMMM\ MM\  MMMMMM\  MM\ MMMMMMMM\ MMMMMMMM\ MM\      =
 =     MM  |MM  _____|\__MM  __|MM |MM  __MM\ MM |\__MM  __|MM  _____|\MM\     =
 =    MM  / MM |         MM |   MM |MM /  MM |MM |   MM |   MM |       \MM\    =
 =   MM  /  MMMMM\       MM |   MM |MM |  MM |MM |   MM |   MMMMM\      \MM\   =
 =   \MM<   MM  __|      MM |   MM |MM |  MM |MM |   MM |   MM  __|     MM  |  =
 =    \MM\  MM |         MM |   MM |MM MM\MM |MM |   MM |   MM |       MM  /   =
 =     \MM\ MMMMMMMM\    MM |   MM |\MMMMMM / MM |   MM |   MMMMMMMM\ MM  /    =
 =      \__|\________|   \__|   \__| \___MMM\ \__|   \__|   \________|\__/     =
 =                                       \___|                                 =
 ===============================================================================
 ===============================================================================
 =========== Energy    Transfer       Quantum        Time   Evolution ==========
 ================================= Version 1.0.0 ===============================
 =================================== April 2023 ================================
 ===============================================================================
 ====================== J.-R. Li, Y. Zhai, Z. Qu, and H. Li ====================
 ===============================================================================
    )foo"
     << std::endl;
  return 0;
}
} // namespace etqte

int
main(int argc, char** argv)
{
  using namespace etqte;
  etqte::printbanner(std::cout);
  ensemble impl{ std::string(argv[1]) };
  impl.driver2023();
  return 0;
}
