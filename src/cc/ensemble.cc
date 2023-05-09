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
#include "ensemble.hh"
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

namespace etqte {

int Ensemble::writerhos() {
  for (size_t sys = 0; sys != nsys_; ++sys) {
    // std::cout << fmt::format("{}{}{}.etcc", prefix_, sys + 1, suffix_)
    // << std::endl;
    systems_[sys].writerhos(
        fmt::format("{}{}{}.etcc", prefix_, sys + 1, suffix_));
  }
  return 0;
}

// Read the input file, a sample of which is provided.
int Ensemble::readsettings(const std::string &fname) {
  settingfname_ = fname;
  try {
    auto settings = toml::parse_file(fname);
    nstate_ =
        settings["etqtesettings"]["number_of_states_selected"].value_or(0);
    nstatep_ =
        settings["etqtesettings"]["number_of_states_presented"].value_or(0);
    prefix_ = settings["etqtesettings"]["prefix_of_matrix_files"].value_or("");
    suffix_ = settings["etqtesettings"]["suffix_of_matrix_files"].value_or("");
    nsys_ = settings["etqtesettings"]["number_of_systems"].value_or(0);
    nstep_ = settings["etqtesettings"]["number_of_steps"].value_or(0);
    tstep_ = settings["etqtesettings"]["step_length"].value_or(0.0);
    write_isosys_ = settings["etqtesettings"]["write_isosys"].value_or(false);

    std::string istate =
        settings["etqtesettings"]["initial_state"].value_or("");

    if (auto arr =
            settings["etqtesettings"]["state_names_presented"].as_array()) {

      statenamesp_.clear();
      for (auto &&i : *arr) {
        statenamesp_.push_back(i.as_string()->value_or(""));
      }
    } else {
      throw std::runtime_error("Please provide the states' names.");
    }

    if (auto arr =
            settings["etqtesettings"]["state_names_selected"].as_array()) {
      statenames_.clear();
      for (auto &&i : *arr) {
        statenames_.push_back(i.as_string()->value_or(""));
      }
    }

    if (statenames_.size() != nstate_) {
      throw std::runtime_error("The size of state_names_selected does not "
                               "consist with number_of_states_selected.");
    }

    if (statenamesp_.size() != nstatep_) {
      throw std::runtime_error("The size of state_names_presented does not "
                               "consist with number_of_states_presented.");
    }

    if (nstate_ > nstatep_) {
      throw std::runtime_error(
          "number_of_states_selected exceeds number_of_states_presented.");
    }

    if (nstate_ <= 0) {
      // Inthis case, treat all states as selected.
      nstate_ = nstatep_;
      statenames_ = statenamesp_;
      for (size_t i = 0; i != nstate_; ++i) {
        states_selected_.push_back(i);
        if (statenames_[i] == istate) {
          this->istate_ = i;
        }
      }
    } else if (nstate_ != nstatep_) {
      // Selected the states we are interested in.
      std::map<std::string, size_t> statemap;
      for (size_t i = 0; i != nstatep_; ++i) {
        statemap[statenamesp_[i]] = i;
      }
      for (size_t i = 0; i != nstate_; ++i) {
        states_selected_.push_back(statemap[statenames_[i]]);
        if (statenames_[i] == istate) {
          this->istate_ = i;
        }
      }
    }

  } catch (const toml::parse_error &err) {
    std::cout << "Error caught: " << err.what() << std::endl;
    std::cout << "*** This error is raised by toml++ libaray, which "
                 "indicates that your input file does not meet the syntax "
                 "requirement of toml format."
              << std::endl;
    throw;
  } catch (const std::runtime_error &err) {
    std::cout << "Error caught: " << err.what() << std::endl;
    std::cout << "*** This error is raised by the main logic of <ET|Q|TE>."
              << std::endl;
    throw;
  }

  return 0;
}

// This function read in the elements of effective Hamiltonian
// Here I follow Jia-Rui's convention.
int Ensemble::preparesystems() {
  for (size_t sys = 0; sys != nsys_; ++sys) {
    std::string fname = fmt::format("{}{}{}", prefix_, sys + 1, suffix_);
    systems_.push_back(
        Isosys(nstate_, statenames_,
               readhamiltonian(fname, nstatep_, states_selected_)));
  }
  return 0;
}

// The quantum time evolution function
int Ensemble::timeevolution() {
  for (auto &&sys : systems_) {
    sys.eigenstates();
    sys.timeevolution(istate_, nstep_, tstep_);
  }
  return 0;
}

Ensemble::Ensemble(const std::string &input) {
  readsettings(input);
  return;
}

Eigen::MatrixXd
Ensemble::readhamiltonian(const std::string &matrixfname, size_t nstate,
                          const std::vector<size_t> &states_selected) {
  std::ifstream io(matrixfname);
  Eigen::MatrixXd Hp(nstate, nstate);
  /*
  std::string tmp;
  // This read hamiltonian version only deal with the upper triangle of the
  // Hamiltonian.
  for (size_t i = 0; i != nstate; ++i) {
    for (size_t j = i; j != nstate; ++j) {
      // However, Eigen lib only using the lower triangle
      io >> tmp >> tmp >> Hp(j, i);
    }
  }
  */

  // My new format of Hamiltonian
  // Upper triangle matrix is recorded
  for (size_t i = 0; i != nstate; ++i) {
    for (size_t j = 0; j != nstate; ++j) {
      // However, Eigen lib only using the lower triangle
      io >> Hp(j, i);
    }
  }

  auto &&H = Hp(states_selected, states_selected);

  // Apparently the matrices are in eV, Errrrr...
  // Turn hamiltonian from eV to Hartree. Idiot.
  H /= 27.21138602;

  /*
  std::cout << "H " << matrixfname << std::endl;
  std::cout << H << std::endl;
  std::cout << "states_selected" << std::endl;
  for (auto &&i : states_selected) {
    std::cout << i << std::endl;
  }
  std::cout << H(states_selected, states_selected) << std::endl;
  */

  /*
  Eigen::VectorXd sign = H.col(0).unaryExpr(std::ref(etqte::sign_func));
  H = sign.asDiagonal() * H * sign.asDiagonal();
  std::cout << matrixfname << std::endl;
  std::cout << H << std::endl;
  */
  return H;
}

int Ensemble::ensemblerho() {
  averrho_.resize(nstep_);
  for (size_t step = 0; step != nstep_; ++step) {
    averrho_[step] = Eigen::VectorXd::Zero(nstate_);
    for (size_t sys = 0; sys != nsys_; ++sys) {
      averrho_[step] += systems_[sys].rho()[step];
    }
    averrho_[step] /= nsys_;
  }
  averrho0_ = Eigen::VectorXd::Zero(nstate_);
  for (size_t sys = 0; sys != nsys_; ++sys) {
    averrho0_ += systems_[sys].rho0();
  }
  averrho0_ /= nsys_;
  return 0;
}

Eigen::MatrixXd Ensemble::ensembleH() {
  // Note that this is not right according to our theory...
  // But in our paper we need it
  Eigen::MatrixXd averH = Eigen::MatrixXd::Zero(nstate_, nstate_);
  for (size_t i = 0; i != nsys_; ++i) {
    averH += systems_[i].H();
  }
  averH /= nsys_;
  //  add the environmental effect temporarily, need refactoring later.
  // averH(0,0) += -0.1027/27.21138602;
  // averH(1,1) += -0.1486/27.21138602;
  return averH;
}

int Ensemble::driver2023() {
  std::cout << "*** Read Hamiltonian matrices in...\n";
  preparesystems();
  std::cout << "*** Done.\n" << std::flush;

  {
    std::cout << "*** Start time evolution...\n";
    auto &&t0 = std::chrono::high_resolution_clock::now();
    timeevolution();
    auto &&t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> &&difft = t1 - t0;
    std::cout << "*** End time evolution.\n"
              << "*** Job is done in " << difft.count() << " sec." << std::endl;
  }
  if (write_isosys_) {
    auto &&t0 = std::chrono::high_resolution_clock::now();
    writerhos();
    auto &&t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> &&difft = t1 - t0;
    std::cout << "*** Time evolution of each sampled system has been written "
                 "to .etcc files.\n"
              << "*** Job is done in " << difft.count() << " sec." << std::endl;
  }

  std::cout << "*** Start ensemble average..." << std::endl;
  ensemblerho();

  Isosys::writetimeseries(fmt::format("{}{}{}.{}.etcc", prefix_, 0, suffix_, settingfname_),
                          nstate_, statenames_, nstep_, tstep_, averrho_,
                          averrho0_);
  std::cout << "*** Time evolution of the whole ensemble has been written to "
               "the zeroth .etcc files."
            << std::endl;

  return 0;
}

} // namespace etqte
