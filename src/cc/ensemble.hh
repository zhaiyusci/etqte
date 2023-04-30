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

#ifndef __ETQTE_ENSEMBLE_HH
#define __ETQTE_ENSEMBLE_HH
#include "isosys.hh"
#include <eigen3/Eigen/Dense>
#include <string>
#include <vector>

namespace etqte {
// inline double sign_func(double x) { return x >= 0 ? 1 : -1; }

class Ensemble {
private:
  // General settings
  size_t nstate_, nsys_, istate_, nstep_, nstatep_;
  std::string prefix_, suffix_;
  std::vector<std::string> statenames_;
  std::vector<std::string> statenamesp_;
  std::vector<size_t> states_selected_;
  double tstep_;

  // Workspace
  std::vector<Isosys> systems_;
  std::vector<Eigen::VectorXd> averrho_;
  Eigen::VectorXd averrho0_;
  bool write_isosys_;

public:
  inline size_t nstate() { return nstate_; }
  inline size_t nsys() { return nsys_; }
  inline size_t nstep() { return nstep_; }
  inline size_t istate() { return istate_; }
  inline double tstep() { return tstep_; }
  inline const std::vector<std::string> &statenames() { return statenames_; }
  int writerhos();

  // Read the input file, a sample of which is provided.
  int readsettings(const std::string &fname);

  // This function read in the elements of effective Hamiltonian
  // Here I follow Jia-Rui's convention.
  int preparesystems();

  // The quantum time evolution function
  int timeevolution();

  Ensemble(const std::string &input);

  Eigen::MatrixXd readhamiltonian(const std::string &matrixfname, size_t nstate,
                                  const std::vector<size_t> &states_selected);

  int ensemblerho();

  Eigen::MatrixXd ensembleH();

  int driver2023();
};

} // namespace etqte
#endif
