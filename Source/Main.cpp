#include "StdAfx.hpp"

#include "Clock.hpp"
#include "Configuration.hpp"
#include "MeshsizeFactory.hpp"
#include "ParticleSimulation.hpp"
#include "Simulation.hpp"
#include "TurbulentSimulation.hpp"

#include "ParallelManagers/PetscParallelConfiguration.hpp"

#ifndef NDEBUG
#include <cfenv>
#pragma STDC FENV_ACCESS ON
#endif
int main(int argc, char* argv[]) {
  spdlog::set_level(spdlog::level::info);

  // Parallelisation related. Initialise and identify.
  // ---------------------------------------------------
  int rank  = -1; // This is the processor's identifier
  int nproc = 0;  // Number of processors in the group
#ifdef ENABLE_PETSC
  PetscInitialize(&argc, &argv, "petsc_commandline_arg", PETSC_NULL);
#else
  MPI_Init(&argc, &argv);
#endif

  MPI_Comm_size(PETSC_COMM_WORLD, &nproc);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  spdlog::info("Rank: {}, Nproc: {}", rank, nproc);
  //----------------------------------------------------

#ifdef ENABLE_SINGLE_PRECISION
  spdlog::warn(
    "Using single floating-point precision for benchmarking only; make sure to switch to double precision "
    "for production runs."
  );
#else
  spdlog::info("Using double floating-point precision");
#endif

#ifndef NDEBUG
  spdlog::warn("Running in Debug mode; make sure to switch to Release mode for production/benchmark runs.");
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#else
  spdlog::info("Running in Release mode");
#endif

  if (!argv[1]) {
    spdlog::error("You need to pass a configuration file: mpirun -np 1 ./NS-EOF ExampleCases/Cavity2D.xml.");
    throw std::runtime_error("Argument parsing error");
  }

  // Read configuration and store information in parameters object
  Configuration configuration(argv[1]);
  Parameters    parameters;
  configuration.loadParameters(parameters);
  ParallelManagers::PetscParallelConfiguration parallelConfiguration(parameters);
  MeshsizeFactory::getInstance().initMeshsize(parameters);
  FlowField*          flowField          = NULL;
  Simulation*         simulation         = NULL;
  ParticleSimulation* particleSimulation = NULL;

  spdlog::debug(
    "Processor {} with index {}, {}, {} is computing the size of its subdomain and obtains {}, {} and {}.",
    parameters.parallel.rank,
    parameters.parallel.indices[0],
    parameters.parallel.indices[1],
    parameters.parallel.indices[2],
    parameters.parallel.localSize[0],
    parameters.parallel.localSize[1],
    parameters.parallel.localSize[2]
  );
  spdlog::debug("Left neighbour: {}, right neighbour: {}", parameters.parallel.leftNb, parameters.parallel.rightNb);
  spdlog::debug("Top neighbour: {}, bottom neighbour: {}", parameters.parallel.topNb, parameters.parallel.bottomNb);
  spdlog::debug("Front neighbour: {}, back neighbour: {}", parameters.parallel.frontNb, parameters.parallel.backNb);
  spdlog::debug(
    "Min. meshsizes: {}, {}, {}",
    parameters.meshsize->getDxMin(),
    parameters.meshsize->getDyMin(),
    parameters.meshsize->getDzMin()
  );

  // Initialise simulation
  if (parameters.simulation.type == "turbulence") {
    // WS2: initialise turbulent flow field and turbulent simulation object
    if (rank == 0) {
      spdlog::info("Start Turbulent simulation in {}D", parameters.geometry.dim);
    }
    flowField = new FlowField(parameters);
    if (flowField == NULL) {
      throw std::runtime_error("flowField == NULL!");
    }
    simulation = new TurbulentSimulation(parameters, *flowField);
  } else if (parameters.simulation.type == "dns") {
    if (rank == 0) {
      spdlog::info("Start DNS simulation in {}D", parameters.geometry.dim);
    }
    flowField = new FlowField(parameters);
    if (flowField == NULL) {
      throw std::runtime_error("flowField == NULL!");
    }
    simulation = new Simulation(parameters, *flowField);
  } else {
    throw std::runtime_error("Unknown simulation type! Currently supported: dns, turbulence");
  }

  // Call initialisation of simulation (initialise flow field)
  if (simulation == NULL) {
    throw std::runtime_error("simulation == NULL!");
  }

  simulation->initializeFlowField();

  if (parameters.particles.enable == true) {
    particleSimulation = new ParticleSimulation(parameters, *flowField);
    particleSimulation->initializeParticles();
  }

  // flowField->getFlags().show();

  RealType time       = 0.0;
  RealType timeVtk    = parameters.vtk.interval;
  RealType timeStdOut = parameters.stdOut.interval;
  RealType timeInject = parameters.particles.injectInterval;
  int      timeSteps  = 0;

  // Plot initial state
#ifndef DISABLE_OUTPUT
  simulation->plotVTK(timeSteps, time);
  if (particleSimulation)
    particleSimulation->plot(timeSteps, time);
#endif

  Clock clock;
  // Time loop
  while (time < parameters.simulation.finalTime) {
    simulation->solveTimestep();
    if (particleSimulation) {
      particleSimulation->solveTimestep();
      if (timeInject <= time) {
        particleSimulation->initializeParticles();
        timeInject += parameters.particles.injectInterval;
      }
    }
    timeSteps++;
    time += parameters.timestep.dt;

    if ((rank == 0) && (timeStdOut <= time)) {
      spdlog::info("Current time: {}\tTimestep: {}", time, parameters.timestep.dt);
      timeStdOut += parameters.stdOut.interval;
    }

    if (timeVtk <= time) {
#ifndef DISABLE_OUTPUT
      simulation->plotVTK(timeSteps, time);
      if (particleSimulation) {
        particleSimulation->plot(timeSteps, time);
      }
#endif
      timeVtk += parameters.vtk.interval;
    }
  }

  if (parameters.parallel.rank == 0)
    spdlog::info("Finished simulation with a duration of {}ns", clock.getTime());

    // Plot final solution
#ifndef DISABLE_OUTPUT
  simulation->plotVTK(timeSteps, time);
#endif

  delete simulation;
  simulation = NULL;

  delete particleSimulation;
  particleSimulation = NULL;

  delete flowField;
  flowField = NULL;

#ifdef ENABLE_PETSC
  PetscFinalize();
#else
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}
