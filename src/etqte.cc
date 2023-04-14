#define TOML_ENABLE_FLOAT16 0
#define TOML_ENABLE_UNRELEASED_FEATURES 0
#define TOML_ENABLE_WINDOWS_COMPAT 0
#define TOML_ENABLE_FORMATTERS 0
#include "toml.hh"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h>
#include <fmt/core.h>
#include <fmt/compile.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>

class etqte {
private:
  // General settings
  size_t nstate, nsys, istate, nstep;
  std::string prefix, suffix;
  double tstep;

  // Workspace
  std::vector<Eigen::MatrixXd> H;
  std::vector<Eigen::VectorXd> E;
  std::vector<Eigen::MatrixXd> psi;
  std::vector<std::vector<Eigen::VectorXd>> rho;
  std::vector<Eigen::VectorXd> rho0;
  std::vector<Eigen::VectorXd> averrho;
  Eigen::VectorXd averrho0;

public:
  int writetimeseries(int sys, std::vector<Eigen::VectorXd> &instval,
                      Eigen::VectorXd &infval) {
    FILE *output =
        fopen(fmt::format(FMT_COMPILE("{}{}{}.etcc"), prefix, sys + 1, suffix).c_str(), "w");
    fmt::print(output, "#             time");
    for (size_t state = 0; state != nstate; ++state) {
      fmt::print(output, FMT_COMPILE("      phi_{:3d}"), state + 1);
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
    fmt::print(output, FMT_COMPILE("{:18.5f}"), nstep * tstep * 2);
    for (size_t state = 0; state != nstate; ++state) {
      fmt::print(output, FMT_COMPILE("{:13.5f}"), infval[state]);
    }
    fmt::print(output, "\n");

    return 0;
  }

  // Read the input file, a sample of which is provided.
  int readsettings(const std::string &fname) {
    auto settings = toml::parse_file(fname);
    nstate = settings["etqtesettings"]["number_of_states"].value_or(0);
    prefix = settings["etqtesettings"]["prefix_of_matrix_files"].value_or("");
    suffix = settings["etqtesettings"]["suffix_of_matrix_files"].value_or("");
    nsys = settings["etqtesettings"]["number_of_systems"].value_or(0);
    nstep = settings["etqtesettings"]["number_of_steps"].value_or(0);
    tstep = settings["etqtesettings"]["step_length"].value_or(0.0);
    istate = settings["etqtesettings"]["initial_state"].value_or(0) - 1;

    // Reserve the workspace...
    H.reserve(nsys);
    E.reserve(nsys);
    psi.reserve(nsys);
    rho.reserve(nsys);
    rho0.reserve(nsys);
    return 0;
  }

  // This function read in the elements of effective Hamiltonian
  // Here I follow Jia-Rui's convention.
  int readhamiltonian(size_t sys) {
    std::ifstream io(fmt::format("{}{}{}", prefix, sys + 1, suffix));
    int tmp;
    Eigen::MatrixXd H(nstate, nstate);
    for (size_t i = 0; i != nstate; ++i) {
      for (size_t j = 0; j != nstate; ++j) {
        io >> tmp >> tmp >> H(i, j);
      }
    }
    this->H.push_back(H);
    io.close();
    return 0;
  }

  // The quantum time evolution function
  int timeevolution(size_t sys) {
    const size_t upsize = nstate * (nstate - 1) / 2;
    Eigen::MatrixXd param(nstate, upsize);
    Eigen::VectorXd DeltaE(upsize);
    // Eigen::VectorXd omega(upsize);

    Eigen::VectorXd E(nstate);
    Eigen::MatrixXd psi(nstate, nstate);
    std::vector<Eigen::VectorXd> rho(nstep, Eigen::VectorXd(nstate));
    Eigen::VectorXd rho0(Eigen::VectorXd::Zero(nstate));

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H[sys]);
    E = es.eigenvalues();
    psi = es.eigenvectors();

    for (size_t fstate = 0; fstate != nstate; ++fstate) {
      for (size_t i = 0; i != nstate; ++i) {
        rho0(fstate) +=
            std::pow(psi(istate, i), 2) * std::pow(psi(fstate, i), 2);
      }
      size_t iii = 0;
      for (size_t i = 0; i != nstate; ++i) {
        for (size_t j = i + 1; j != nstate; ++j) {
          param(fstate, iii) =
              psi(istate, i) * psi(istate, j) * psi(fstate, i) * psi(fstate, j);
          ++iii;
        }
      }
    }
    param *= 2.0;
    size_t iii = 0;
    for (size_t i = 0; i != nstate; ++i) {
      for (size_t j = i + 1; j != nstate; ++j) {
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
    for (size_t step = 0; step != nstep; ++step) {
      rho[step] = param * chebyshev_k + rho0;
      chebyshev_k_1 =
          (2 * omega.array() * chebyshev_k.array()).matrix() - chebyshev_k_1;
      chebyshev_k += chebyshev_k_1;
      chebyshev_k_1 = chebyshev_k - chebyshev_k_1;
      chebyshev_k -= chebyshev_k_1;
    }

    this->E.push_back(E);
    this->psi.push_back(psi);
    this->rho0.push_back(rho0);
    this->rho.push_back(rho);
    return 0;
  }

  static int printbanner(std::ostream &io) {
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
  etqte(const std::string &input) {
    readsettings(input);
    return;
  }
  int driver2023() {
    std::cout << "*** Read Hamiltonian matrices in...\n";
    for (size_t sys = 0; sys != nsys; ++sys) {
      readhamiltonian(sys);

      // Apparently the matrices are in eV, Errrrr...
      // Turn hamiltonian from eV to Hartree. Idiot.
      H[sys] /= 27.21138602;
    }
    std::cout << "*** Done.\n" << std::flush;

    std::cout << "*** Start time evolution...\n";
    auto &&t0 = std::chrono::high_resolution_clock::now();
    for (size_t sys = 0; sys != nsys; ++sys) {
      timeevolution(sys);
    }
    auto &&t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> &&difft = t1-t0;
    std::cout << "*** End time evolution. Job is done in " << difft.count() << " sec." << std::endl;
    for (size_t sys = 0; sys != nsys; ++sys) {
      writetimeseries(sys, rho[sys], rho0[sys]);
    }
    std::cout << "*** Time evolution of each sampled system has been written "
                 "to .etout files."
              << std::endl;
    if (false) {
      for (size_t sys = 0; sys != nsys; ++sys) {

        std::cout << "Hamiltonian of system " << sys + 1 << "\n";
        std::cout << H[sys] << '\n';

        std::cout << "Eigenvectors of system " << sys + 1
                  << "(phi representation) \n";
        std::cout << psi[sys] << '\n';
      }
      std::cout << std::endl;
    }

    std::cout << "*** Start ensemble average..." << std::endl;
    averrho.resize(nstep);
    for (size_t step = 0; step != nstep; ++step) {
      averrho[step] = Eigen::VectorXd::Zero(nstate);
      for (size_t sys = 0; sys != nsys; ++sys) {
        averrho[step] += rho[sys][step];
      }
      averrho[step] /= nsys;
    }
    averrho0 = Eigen::VectorXd::Zero(nstate);
    for (size_t sys = 0; sys != nsys; ++sys) {
      averrho0 += rho0[sys];
    }
    averrho0 /= nsys;
    writetimeseries(-1, averrho, averrho0);
    std::cout << "*** Time evolution of the whole ensemble has been written to "
                 "the zeroth .etout files."
              << std::endl;

    return 0;
  }
};

int main(int argc, char **argv) {
  etqte::printbanner(std::cout);
  etqte impl{std::string(argv[1])};
  impl.driver2023();
  return 0;
}
