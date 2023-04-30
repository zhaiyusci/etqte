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
#include "isosys.hh"

namespace etqte {
  Isosys::Isosys(size_t nstate, const std::vector<std::string> &statenames,
         const Eigen::MatrixXd &H)
      : nstate_(nstate), statenames_(statenames), H_(H), E_(), psi_(),
        istate_(0), nstep_(0), tstep_(0.0), rho_(), rho0_() {}

  int Isosys::writetimeseries(std::string fname, size_t nstate,
                             const std::vector<std::string> &statenames,
                             size_t nstep, double tstep,
                             const std::vector<Eigen::VectorXd> &instval,
                             const Eigen::VectorXd &infval) {
    FILE *output = fopen(fname.c_str(), "w");
    fmt::print(output, "#             time");
    for (size_t state = 0; state != nstate; ++state) {
      fmt::print(output, FMT_COMPILE("{:>13s}"), statenames[state]);
    }
    fmt::print(output, "\n");
    for (size_t step = 0; step != nstep; ++step) {
      fmt::print(output, FMT_COMPILE("{:18.5f}"), step * tstep);
      for (size_t state = 0; state != nstate; ++state) {
        fmt::print(output, FMT_COMPILE("{:13.5f}"), instval[step][state]);
      }
      fmt::print(output, "\n");
    }
    // We use a long enough time to present an infinite long time
    fmt::print(output, "\n");
    fmt::print(output, FMT_COMPILE("{:18.5f}"), nstep * tstep * 1.2);
    for (size_t state = 0; state != nstate; ++state) {
      fmt::print(output, FMT_COMPILE("{:13.5f}"), infval[state]);
    }
    fmt::print(output, "\n");
    fmt::print(output, FMT_COMPILE("{:18.5f}"), nstep * tstep * 2);
    for (size_t state = 0; state != nstate; ++state) {
      fmt::print(output, FMT_COMPILE("{:13.5f}"), infval[state]);
    }
    fmt::print(output, "\n");

    fclose(output);

    return 0;
  }

  int Isosys::writerhos(const std::string &fname) {
    writetimeseries(fname, nstate_, statenames_, nstep_, tstep_, rho_, rho0_);
    return 0;
  }

  int Isosys::eigenstates() {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H_);
    E_ = es.eigenvalues();
    psi_ = es.eigenvectors();
    return 0;
  }

  // The quantum time evolution function
  int Isosys::timeevolution(size_t istate, size_t nstep, double tstep) {
    this->istate_ = istate;
    this->nstep_ = nstep;
    this->tstep_ = tstep;
    const size_t upsize = nstate_ * (nstate_ - 1) / 2;

    // Actually +2 is good enough, but I think +5 is no harm.
    rho_.clear();
    rho_.reserve(nstep + 5);
    rho0_ = Eigen::VectorXd::Zero(nstate_);

    Eigen::MatrixXd param(nstate_, upsize);
    Eigen::VectorXd DeltaE(upsize);

    for (size_t fstate = 0; fstate != nstate_; ++fstate) {
      for (size_t i = 0; i != nstate_; ++i) {
        rho0_(fstate) +=
            std::pow(psi_(istate, i), 2) * std::pow(psi_(fstate, i), 2);
      }
      size_t iii = 0;
      for (size_t i = 0; i != nstate_; ++i) {
        for (size_t j = i + 1; j != nstate_; ++j) {
          param(fstate, iii) = psi_(istate, i) * psi_(istate, j) *
                               psi_(fstate, i) * psi_(fstate, j);
          ++iii;
        }
      }
    }
    param *= 2.0;
    size_t iii = 0;
    for (size_t i = 0; i != nstate_; ++i) {
      for (size_t j = i + 1; j != nstate_; ++j) {
        DeltaE(iii) = E_(j) - E_(i);
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
    auto *chb_k = &chebyshev_k;
    auto *chb_k_1 = &chebyshev_k_1;
    for (size_t step = 0; step != nstep; ++step) {
      Eigen::VectorXd rho = param * (*chb_k) + rho0_;
      this->rho_.push_back(std::move(rho));
      *chb_k_1 = (2 * omega.array() * chb_k->array()).matrix() - (*chb_k_1);
      std::swap(chb_k, chb_k_1);
    }

    return 0;
  }

  Isosys Isosys::subsystem(std::vector<size_t> &states_selected) {
    Isosys subsys = *this;
    subsys.nstate_ = states_selected.size();
    subsys.H_ = this->H_(states_selected, states_selected);
    subsys.statenames_.clear();
    for (auto &&i : states_selected) {
      subsys.statenames_.push_back(this->statenames_[i]);
    }
    return subsys;
  }


} // namespace etqte

