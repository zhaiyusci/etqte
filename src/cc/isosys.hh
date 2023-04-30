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

#ifndef __ETQTE_ISOSYS_HH
#define __ETQTE_ISOSYS_HH
#include <eigen3/Eigen/Dense>
#include <vector>
#include <string>

namespace etqte {
// inline double sign_func(double x) { return x >= 0 ? 1 : -1; }
class Isosys {
private:
  // General settings
  size_t nstate_;
  std::vector<std::string> statenames_;

  // Workspace
  Eigen::MatrixXd H_;
  Eigen::VectorXd E_;
  Eigen::MatrixXd psi_;
  size_t istate_;
  size_t nstep_;
  double tstep_;
  std::vector<Eigen::VectorXd> rho_;
  Eigen::VectorXd rho0_;

public:
  // Some getters
  inline size_t nstate() { return nstate_; }
  inline const std::vector<std::string> &statenames() { return statenames_; }
  inline const Eigen::MatrixXd &H() { return H_; }
  inline const Eigen::VectorXd &E() { return E_; }
  inline const Eigen::MatrixXd &psi() { return psi_; }
  inline const std::vector<Eigen::VectorXd> &rho() { return rho_; }
  inline const Eigen::VectorXd &rho0() { return rho0_; }

  // Constructor
  Isosys(size_t nstate, const std::vector<std::string> &statenames,
         const Eigen::MatrixXd &H);

  static int writetimeseries(std::string fname, size_t nstate,
                             const std::vector<std::string> &statenames,
                             size_t nstep, double tstep,
                             const std::vector<Eigen::VectorXd> &instval,
                             const Eigen::VectorXd &infval); 

  int writerhos(const std::string &fname) ;

  int eigenstates() ;

  // The quantum time evolution function
  int timeevolution(size_t istate, size_t nstep, double tstep) ;

  Isosys subsystem(std::vector<size_t> &states_selected) ;
};


} // namespace etqte
#endif
