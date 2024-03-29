#include "StdAfx.hpp"

#include "Configuration.hpp"

#include <tinyxml2.h>

void readFloatMandatory(RealType& storage, tinyxml2::XMLElement* node, const char* tag) {
  double value; // Use to be able to select precision
  if (node->QueryDoubleAttribute(tag, &value) != tinyxml2::XML_SUCCESS) {
    throw std::runtime_error("Error while reading mandatory argument");
  } else {
    storage = static_cast<RealType>(value);
  }
}

void readFloatOptional(RealType& storage, tinyxml2::XMLElement* node, const char* tag, RealType defaultValue = 0) {
  double value; // Use to be able to select precision
  int    result = node->QueryDoubleAttribute(tag, &value);
  if (result == tinyxml2::XML_NO_ATTRIBUTE) {
    storage = defaultValue;
  } else if (result == tinyxml2::XML_WRONG_ATTRIBUTE_TYPE) {
    throw std::runtime_error("Error while reading optional argument");
  } else {
    storage = static_cast<RealType>(value);
  }
}

void readIntMandatory(int& storage, tinyxml2::XMLElement* node, const char* tag) {
  int value;
  if (node->QueryIntAttribute(tag, &value) != tinyxml2::XML_SUCCESS) {
    throw std::runtime_error("Error while reading mandatory argument");
  } else {
    storage = value;
  }
}

void readIntOptional(int& storage, tinyxml2::XMLElement* node, const char* tag, int defaultValue = 0) {
  int result = node->QueryIntAttribute(tag, &storage);
  if (result == tinyxml2::XML_NO_ATTRIBUTE) {
    storage = defaultValue;
  } else if (result == tinyxml2::XML_WRONG_ATTRIBUTE_TYPE) {
    throw std::runtime_error("Error while reading optional argument");
  }
}

void readBoolMandatory(bool& storage, tinyxml2::XMLElement* node, const char* tag) {
  bool value;
  if (node->QueryBoolAttribute(tag, &value) != tinyxml2::XML_SUCCESS) {
    throw std::runtime_error("Error while reading mandatory argument");
  } else {
    storage = value;
  }
}

void readBoolOptional(bool& storage, tinyxml2::XMLElement* node, const char* tag, bool defaultValue = false) {
  int result = node->QueryBoolAttribute(tag, &storage);
  if (result == tinyxml2::XML_NO_ATTRIBUTE) {
    storage = defaultValue;
  } else if (result == tinyxml2::XML_WRONG_ATTRIBUTE_TYPE) {
    throw std::runtime_error("Error while reading optional argument");
  }
}

void readStringMandatory(std::string& storage, tinyxml2::XMLElement* node) {
  const char* myText = node->GetText();
  if (myText == NULL) {
    const std::string nodename = node->Name();
    spdlog::error("No string specified for this node: {}", nodename);
    throw std::runtime_error("Error while reading mandatory string");
  } else {
    storage = node->GetText();
    if (!storage.compare("")) {
      throw std::runtime_error("Missing mandatory string!");
    }
  }
}

void readWall(tinyxml2::XMLElement* wall, RealType* vector, RealType& scalar) {
  tinyxml2::XMLElement* quantity = wall->FirstChildElement("vector");
  if (quantity != NULL) {
    readFloatOptional(vector[0], quantity, "x");
    readFloatOptional(vector[1], quantity, "y");
    readFloatOptional(vector[2], quantity, "z");
  }
  quantity = wall->FirstChildElement("scalar");
  if (quantity != NULL) {
    readFloatOptional(scalar, quantity, "value");
  }
}

void broadcastString(std::string& target, const MPI_Comm& communicator, int root = 0) {
  int stringSize = 0, rank = -1;
  MPI_Comm_rank(communicator, &rank);
  if (rank == root) {
    stringSize = static_cast<int>(target.size());
  }
  MPI_Bcast(&stringSize, 1, MPI_INT, 0, communicator);
  char* name = new char[stringSize + 1]; // One more for the null character
  if (rank == root) {
    target.copy(name, stringSize, 0);
  }
  name[stringSize] = '\0';
  MPI_Bcast(name, stringSize + 1, MPI_CHAR, 0, communicator);
  if (rank != root) {
    target = name;
  }
  delete[] name;
  name = NULL;
}

Configuration::Configuration() { filename_ = ""; }

Configuration::Configuration(const std::string& filename) { filename_ = filename; }

void Configuration::setFileName(const std::string& filename) { filename_ = filename; }

void Configuration::loadParameters(Parameters& parameters, const MPI_Comm& communicator) {
  tinyxml2::XMLDocument confFile;
  tinyxml2::XMLElement* node;
  tinyxml2::XMLElement* subNode;

  int rank = -1;
  MPI_Comm_rank(communicator, &rank);

  // We only read on rank 0; afterwards, all configuration parameters are broadcasted to all processes.
  // So, if you add new parameters in the configuration, make sure to broadcast them to the other processes!
  if (rank == 0) {
    // Parse the configuration file and check validity
    confFile.LoadFile(filename_.c_str());
    if (confFile.FirstChildElement() == NULL) {
      throw std::runtime_error("Error parsing the configuration file");
    }

    //--------------------------------------------------
    // Load geometric parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("geometry");

    if (node == NULL) {
      throw std::runtime_error("Error loading geometry properties");
    }

    readIntMandatory(parameters.geometry.sizeX, node, "sizeX");
    readIntMandatory(parameters.geometry.sizeY, node, "sizeY");
    readIntOptional(parameters.geometry.sizeZ, node, "sizeZ");

    if (parameters.geometry.sizeX < 2 || parameters.geometry.sizeY < 2 || parameters.geometry.sizeZ < 0) {
      throw std::runtime_error("Invalid size specified in configuration file");
    }

    parameters.geometry.dim = 0;
    if (node->QueryIntAttribute("dim", &(parameters.geometry.dim)) != tinyxml2::XML_WRONG_ATTRIBUTE_TYPE) {
      if (parameters.geometry.dim == 0) {
        if (parameters.geometry.sizeZ == 0) {
          parameters.geometry.sizeZ = 1;
          parameters.geometry.dim   = 2;
        } else {
          parameters.geometry.dim = 3;
        }
      }
    }

    if (parameters.geometry.dim == 3 && parameters.geometry.sizeZ == 1) {
      throw std::runtime_error("Inconsistent data: 3D geometry specified with Z size zero");
    }

    if (parameters.geometry.dim == 2 && parameters.geometry.sizeZ != 1) {
      parameters.geometry.sizeZ = 1;
    }

    // Determine the sizes of the cells
    readFloatMandatory(parameters.geometry.lengthX, node, "lengthX");
    readFloatMandatory(parameters.geometry.lengthY, node, "lengthY");
    readFloatMandatory(parameters.geometry.lengthZ, node, "lengthZ");
    // Read geometry->meshsize parameters
    std::string meshsizeType = "";
    subNode                  = node->FirstChildElement("mesh");
    readStringMandatory(meshsizeType, subNode);
    if (meshsizeType == "uniform") {
      parameters.geometry.meshsizeType = Uniform;
    } else if (meshsizeType == "stretched") {
      parameters.geometry.meshsizeType = TanhStretching;
      bool buffer                      = false;
      readBoolMandatory(buffer, node, "stretchX");
      parameters.geometry.stretchX = static_cast<int>(buffer);
      readBoolMandatory(buffer, node, "stretchY");
      parameters.geometry.stretchY = static_cast<int>(buffer);
      if (parameters.geometry.dim == 3) {
        readBoolMandatory(buffer, node, "stretchZ");
        parameters.geometry.stretchZ = static_cast<int>(buffer);
      } else {
        parameters.geometry.stretchZ = false;
      }
    } else {
      throw std::runtime_error("Unknown 'mesh'!");
    }

    // Now, the size of the elements should be set

    dim_ = parameters.geometry.dim;

    //--------------------------------------------------
    // Timestep parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("timestep");

    if (node == NULL) {
      throw std::runtime_error("Error loading timestep parameters");
    }

    readFloatOptional(parameters.timestep.dt, node, "dt", 1);
    readFloatOptional(parameters.timestep.tau, node, "tau", 0.5);

    //--------------------------------------------------
    // Flow parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("flow");

    if (node == NULL) {
      throw std::runtime_error("Error loading flow parameters");
    }

    readFloatMandatory(parameters.flow.Re, node, "Re");

    //--------------------------------------------------
    // Solver parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("solver");

    if (node == NULL) {
      throw std::runtime_error("Error loading solver parameters");
    }

    readFloatMandatory(parameters.solver.gamma, node, "gamma");
    readIntOptional(parameters.solver.maxIterations, node, "maxIterations");

    //--------------------------------------------------
    // Environmental parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("environment");

    if (node == NULL) {
      throw std::runtime_error("Error loading environmental parameters");
    }

    readFloatOptional(parameters.environment.gx, node, "gx");
    readFloatOptional(parameters.environment.gy, node, "gy");
    readFloatOptional(parameters.environment.gz, node, "gz");

    //--------------------------------------------------
    // Simulation parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("simulation");

    if (node == NULL) {
      throw std::runtime_error("Error loading simulation parameters");
    }

    readFloatMandatory(parameters.simulation.finalTime, node, "finalTime");

    subNode = node->FirstChildElement("type");
    if (subNode != NULL) {
      readStringMandatory(parameters.simulation.type, subNode);
    } else {
      throw std::runtime_error("Missing type in simulation parameters");
    }

    subNode = node->FirstChildElement("scenario");
    if (subNode != NULL) {
      readStringMandatory(parameters.simulation.scenario, subNode);
    } else {
      throw std::runtime_error("Missing scenario in simulation parameters");
    }

    //--------------------------------------------------
    // VTK parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("vtk");

    if (node == NULL) {
      throw std::runtime_error("Error loading VTK parameters");
    }

    readFloatOptional(parameters.vtk.interval, node, "interval");
    readStringMandatory(parameters.vtk.prefix, node);

    //--------------------------------------------------
    // StdOut parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("stdOut");

    if (node == NULL) {
      throw std::runtime_error("Error loading StdOut parameters");
    }

    // If no value given, print every step
    readFloatOptional(parameters.stdOut.interval, node, "interval", 1);

    //--------------------------------------------------
    // Parallel parameters
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("parallel");

    if (node == NULL) {
      throw std::runtime_error("Error loading parallel parameters");
    }

    readIntOptional(parameters.parallel.numProcessors[0], node, "numProcessorsX", 1);
    readIntOptional(parameters.parallel.numProcessors[1], node, "numProcessorsY", 1);
    readIntOptional(parameters.parallel.numProcessors[2], node, "numProcessorsZ", 1);

    // Start neighbors on null in case that no parallel configuration is used later.
    parameters.parallel.leftNb   = MPI_PROC_NULL;
    parameters.parallel.rightNb  = MPI_PROC_NULL;
    parameters.parallel.bottomNb = MPI_PROC_NULL;
    parameters.parallel.topNb    = MPI_PROC_NULL;
    parameters.parallel.frontNb  = MPI_PROC_NULL;
    parameters.parallel.backNb   = MPI_PROC_NULL;

    // Yet more parameters initialized in case that no parallel configuration is applied
    parameters.parallel.localSize[0] = parameters.geometry.sizeX;
    parameters.parallel.localSize[1] = parameters.geometry.sizeY;
    parameters.parallel.localSize[2] = parameters.geometry.sizeZ;

    parameters.parallel.firstCorner[0] = 0;
    parameters.parallel.firstCorner[1] = 0;
    parameters.parallel.firstCorner[2] = 0;

    // VTK output is named after the rank, so we define it here, again, in case that it's not
    // initialized anywhere else.
    parameters.parallel.rank = rank;

    //--------------------------------------------------
    // Walls
    //--------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("walls");

    if (node == NULL) {
      throw std::runtime_error("Error loading wall parameters");
    }

    tinyxml2::XMLElement* wall;
    wall = node->FirstChildElement("left");
    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorLeft, parameters.walls.scalarLeft);
    }

    wall = node->FirstChildElement("right");
    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorRight, parameters.walls.scalarRight);
    }

    wall = node->FirstChildElement("bottom");
    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorBottom, parameters.walls.scalarBottom);
    }

    wall = node->FirstChildElement("top");
    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorTop, parameters.walls.scalarTop);
    }

    wall = node->FirstChildElement("front");
    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorFront, parameters.walls.scalarFront);
    }

    wall = node->FirstChildElement("back");
    if (wall != NULL) {
      readWall(wall, parameters.walls.vectorBack, parameters.walls.scalarBack);
    }

    // Set the scalar values to zero;
    // do not set the left pressure value to zero, if we have a pressure-channel
    // scenario -> in this case, we need a fixed pressure value there.
    if (parameters.simulation.scenario != "pressure-channel") {
      parameters.walls.scalarLeft = 0.0;
    }
    parameters.walls.scalarRight  = 0.0;
    parameters.walls.scalarBottom = 0.0;
    parameters.walls.scalarTop    = 0.0;
    parameters.walls.scalarFront  = 0.0;
    parameters.walls.scalarBack   = 0.0;

    //--------------------------------------------------
    // Backward facing step
    //--------------------------------------------------
    parameters.bfStep.xRatio = -1.0;
    parameters.bfStep.yRatio = -1.0;
    node                     = confFile.FirstChildElement()->FirstChildElement("backwardFacingStep");
    if (node != NULL) {
      readFloatMandatory(parameters.bfStep.xRatio, node, "xRatio");
      readFloatMandatory(parameters.bfStep.yRatio, node, "yRatio");
    }

    //------------------------------------------------------
    // Turbulence parameters
    //------------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("turbulence");
    if (node != NULL) {
      subNode = node->FirstChildElement("boundaryLayer");
      if (subNode != NULL) {
        std::string boundaryLayer;
        readStringMandatory(boundaryLayer, subNode);
        if (boundaryLayer == "none") {
          parameters.turbulence.boundaryLayer = 0;
        } else if (boundaryLayer == "laminar") {
          parameters.turbulence.boundaryLayer = 1;
        } else if (boundaryLayer == "turbulent") {
          parameters.turbulence.boundaryLayer = 2;
        } else {
          throw std::runtime_error("Unidentified boundary layer type!");
        }
      } else {
        throw std::runtime_error("Boundary layer type not defined!");
      }
    }

    //------------------------------------------------------
    // Particle simulation parameters
    //------------------------------------------------------
    node = confFile.FirstChildElement()->FirstChildElement("particles");
    if (node != NULL) {
      parameters.particles.enable = true;
      readIntMandatory(parameters.particles.particleCount, node, "particleCount");
      parameters.particles.injectInterval = parameters.simulation.finalTime + 1.0;
      readFloatOptional(parameters.particles.injectInterval, node, "injectInterval");
    }
  }

  // Broadcasting of the values
  MPI_Bcast(&(parameters.geometry.sizeX), 1, MPI_INT, 0, communicator);
  MPI_Bcast(&(parameters.geometry.sizeY), 1, MPI_INT, 0, communicator);
  MPI_Bcast(&(parameters.geometry.sizeZ), 1, MPI_INT, 0, communicator);

  MPI_Bcast(&(parameters.geometry.dim), 1, MPI_INT, 0, communicator);

  MPI_Bcast(&(parameters.geometry.meshsizeType), 1, MPI_INT, 0, communicator);
  MPI_Bcast(&(parameters.geometry.stretchX), 1, MPI_INT, 0, communicator);
  MPI_Bcast(&(parameters.geometry.stretchY), 1, MPI_INT, 0, communicator);
  MPI_Bcast(&(parameters.geometry.stretchZ), 1, MPI_INT, 0, communicator);
  MPI_Bcast(&(parameters.geometry.lengthX), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.geometry.lengthY), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.geometry.lengthZ), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(&(parameters.timestep.dt), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.timestep.tau), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(&(parameters.flow.Re), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(&(parameters.solver.gamma), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.solver.maxIterations), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(&(parameters.environment.gx), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.environment.gy), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.environment.gz), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(&(parameters.simulation.finalTime), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(&(parameters.vtk.interval), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.stdOut.interval), 1, MPI_INT, 0, communicator);

  broadcastString(parameters.vtk.prefix, communicator);
  broadcastString(parameters.simulation.type, communicator);
  broadcastString(parameters.simulation.scenario, communicator);

  MPI_Bcast(&(parameters.bfStep.xRatio), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.bfStep.yRatio), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(parameters.parallel.numProcessors, 3, MPI_INT, 0, communicator);

  MPI_Bcast(&(parameters.walls.scalarLeft), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.walls.scalarRight), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.walls.scalarBottom), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.walls.scalarTop), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.walls.scalarFront), 1, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(&(parameters.walls.scalarBack), 1, MY_MPI_FLOAT, 0, communicator);

  MPI_Bcast(parameters.walls.vectorLeft, 3, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(parameters.walls.vectorRight, 3, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(parameters.walls.vectorBottom, 3, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(parameters.walls.vectorTop, 3, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(parameters.walls.vectorFront, 3, MY_MPI_FLOAT, 0, communicator);
  MPI_Bcast(parameters.walls.vectorBack, 3, MY_MPI_FLOAT, 0, communicator);

  // TODO WS2: broadcast turbulence parameters
  MPI_Bcast(&(parameters.turbulence.boundaryLayer), 1, MPI_INT, 0, communicator);

  MPI_Bcast(&(parameters.particles.enable), 1, MPI_C_BOOL, 0, communicator);
  MPI_Bcast(&(parameters.particles.particleCount), 1, MPI_INT, 0, communicator);
  MPI_Bcast(&(parameters.particles.injectInterval), 1, MY_MPI_FLOAT, 0, communicator);
}
