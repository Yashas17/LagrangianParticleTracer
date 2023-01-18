#include "ParticleSimulation.hpp"

ParticleSimulation::ParticleSimulation(Parameters& parameters, FlowField& flowField):
  parameters_(parameters),
  flowField_(flowField) {}

void ParticleSimulation::initializeParticles() {
  // Number of particles in each direction (2D: y - 3D: y & z)
  const int particleCount = parameters_.particles.particleCount;
  const int dim           = parameters_.geometry.dim;

  const RealType lengthY = (parameters_.bfStep.yRatio <= 0.0)
                             ? parameters_.geometry.lengthY
                             : parameters_.geometry.lengthY * (1 - parameters_.bfStep.yRatio);
  const RealType lengthZ = parameters_.geometry.lengthZ;

  // Assuming uniform spacing between particles
  double spacingY = lengthY / (particleCount + 1);
  double spacingZ = lengthZ / (particleCount + 1);

  RealType x = parameters_.meshsize->getDx(2, 0) / 2;
  RealType y = (parameters_.bfStep.yRatio >= 0) ? (parameters_.bfStep.yRatio * parameters_.geometry.lengthY) : 0;
  std::array<int, 3> index = {2, 0, 0};

  ASSERTION(dim == 2 || dim == 3);

  if (dim == 2) {
    // Uniform distribution of the particles in the y-direction
    // Avoid having particles at exactly the bottom and top walls
    int j = 0;
    for (int p = 0; p < particleCount; p++) {
      y += spacingY;
      while (parameters_.meshsize->getPosY(0, j) + parameters_.meshsize->getDy(0, j) <= y) {
        j++;
      }
      index[1] = j;
      if(x > parameters_.meshsize->getPosX(2, 0)){
        particles_.push_front(Particle(x, y, index, flowField_, parameters_));
      }
    }
  } else {
    RealType z;

    // Uniform distributution of the particles in the YZ-direction
    // Total of particleCount*particleCount particles in the YZ-plane
    // Avoid having particles at exactly the bottom, top, front, and back walls
    int j = 0;
    int k = 0;
    for (int py = 0; py < particleCount; py++) {
      k = 0;
      y += spacingY;
      z = 0;
      while (parameters_.meshsize->getPosY(0, j, k) + parameters_.meshsize->getDy(0, j, k) <= y) {
        j++;
      }
      index[1] = j;
      for (int pz = 0; pz < particleCount; pz++) {
        z += spacingZ;
        while (parameters_.meshsize->getPosZ(0, 0, k) + parameters_.meshsize->getDz(0, 0, k) <= z) {
          k++;
        }
        index[2] = k;
        if(x > parameters_.meshsize->getPosX(2, 0, 0)){
          particles_.push_front(Particle(x, y, z, index, flowField_, parameters_));
        }
      }
    }
  }
}

void ParticleSimulation::solveTimestep() {

  for (auto& particle : particles_) {
    particle.update(parameters_.timestep.dt);
  }

  communicateParticles();
}

void ParticleSimulation::plot(int timeSteps, RealType time){
    const int dim = parameters_.geometry.dim;

    //open file
    std::string prefix = parameters_.vtk.prefix; //read the prefix
    std::string outputFolder = "Output/" + prefix;
    std::stringstream namestream;
    std::string       name;
    std::string grid;
    char        buffer[256];

    namestream.precision(4);
    namestream
    << "Output"
    << "/" << prefix << "/" << prefix << "particles." << parameters_.parallel.rank << "." << timeSteps << ".vtk";
    name = namestream.str();
    std::ofstream ofile;
    ofile.open(name.c_str()); //open the file
    namestream.str("");

    //write header
    ofile << "# vtk DataFile Version 2.0" << std::endl << "NS-EOF" << std::endl << "ASCII" << std::endl << std::endl;

    grid.reserve((ofile.precision() + 6) * particles_.size() * 3);

    sprintf(
        buffer,
        "DATASET POLYDATA\nPOINTS %d float\n",
        particles_.size()
    );
    grid.append(buffer);

    //loop over particles
    if(dim == 3){
        for(auto particle: particles_){
            sprintf(
                buffer,
                "%f %f %f\n",
                particle.getX(),
                particle.getY(),
                particle.getZ()
            );
            grid.append(buffer);
        }
    }
    else{
        for(auto particle: particles_){
            sprintf(
                buffer,
                "%f %f 0.0\n",
                particle.getX(),
                particle.getY()
            );
            grid.append(buffer);
        }
    }
    
    //close file
    grid.append("\n");
    ofile << grid;
    ofile.close();
}

std::vector<RealType> ParticleSimulation::collectLeftBoundaryParticles(){
  const int dim = parameters_.geometry.dim;
  std::vector<RealType> sendBuffer{};

  for(auto& particle: particles_){
    if(particle.getI() < 2){ //left ghost cell territory
      if(dim == 2){
        sendBuffer.push_back(particle.getX());
        sendBuffer.push_back(particle.getY());
        sendBuffer.push_back(particle.getU());
        sendBuffer.push_back(particle.getV());
        sendBuffer.push_back(static_cast<RealType>(particle.getI()));
        sendBuffer.push_back(static_cast<RealType>(particle.getJ()));
      }
      else{
        sendBuffer.push_back(particle.getX());
        sendBuffer.push_back(particle.getY());
        sendBuffer.push_back(particle.getZ());
        sendBuffer.push_back(particle.getU());
        sendBuffer.push_back(particle.getV());
        sendBuffer.push_back(particle.getW());
        sendBuffer.push_back(static_cast<RealType>(particle.getI()));
        sendBuffer.push_back(static_cast<RealType>(particle.getJ()));
        sendBuffer.push_back(static_cast<RealType>(particle.getK()));
      }

      //TODO: delete the particle
    }
  }

  if(sendBuffer.size() == 0) sendBuffer.push_back(0.0);
  return sendBuffer;
}

std::vector<RealType> ParticleSimulation::collectRightBoundaryParticles(){
  const int dim = parameters_.geometry.dim;
  std::vector<RealType> sendBuffer{};

  for(auto& particle: particles_){
    if(particle.getI() >= parameters_.parallel.localSize[0] + 2){ //right ghost cell territory
      if(dim == 2){
        sendBuffer.push_back(particle.getX());
        sendBuffer.push_back(particle.getY());
        sendBuffer.push_back(particle.getU());
        sendBuffer.push_back(particle.getV());
        sendBuffer.push_back(static_cast<RealType>(particle.getI()));
        sendBuffer.push_back(static_cast<RealType>(particle.getJ()));
      }
      else{
        sendBuffer.push_back(particle.getX());
        sendBuffer.push_back(particle.getY());
        sendBuffer.push_back(particle.getZ());
        sendBuffer.push_back(particle.getU());
        sendBuffer.push_back(particle.getV());
        sendBuffer.push_back(particle.getW());
        sendBuffer.push_back(static_cast<RealType>(particle.getI()));
        sendBuffer.push_back(static_cast<RealType>(particle.getJ()));
        sendBuffer.push_back(static_cast<RealType>(particle.getK()));
      }

      //TODO: delete the particle
    }
  }

  if(sendBuffer.size() == 0) sendBuffer.push_back(0.0);
  return sendBuffer;
}

std::vector<RealType> ParticleSimulation::collectBottomBoundaryParticles(){
  const int dim = parameters_.geometry.dim;
  std::vector<RealType> sendBuffer{};

  for(auto& particle: particles_){
    if(particle.getJ() < 2){ //bottom ghost cell territory
      if(dim == 2){
        sendBuffer.push_back(particle.getX());
        sendBuffer.push_back(particle.getY());
        sendBuffer.push_back(particle.getU());
        sendBuffer.push_back(particle.getV());
        sendBuffer.push_back(static_cast<RealType>(particle.getI()));
        sendBuffer.push_back(static_cast<RealType>(particle.getJ()));
      }
      else{
        sendBuffer.push_back(particle.getX());
        sendBuffer.push_back(particle.getY());
        sendBuffer.push_back(particle.getZ());
        sendBuffer.push_back(particle.getU());
        sendBuffer.push_back(particle.getV());
        sendBuffer.push_back(particle.getW());
        sendBuffer.push_back(static_cast<RealType>(particle.getI()));
        sendBuffer.push_back(static_cast<RealType>(particle.getJ()));
        sendBuffer.push_back(static_cast<RealType>(particle.getK()));
      }

      //TODO: delete the particle
    }
  }

  if(sendBuffer.size() == 0) sendBuffer.push_back(0.0);
  return sendBuffer;
}

std::vector<RealType> ParticleSimulation::collectTopBoundaryParticles(){
  const int dim = parameters_.geometry.dim;
  std::vector<RealType> sendBuffer{};

  for(auto& particle: particles_){
    if(particle.getJ() >= parameters_.parallel.localSize[1] + 2){ //top ghost cell territory
      if(dim == 2){
        sendBuffer.push_back(particle.getX());
        sendBuffer.push_back(particle.getY());
        sendBuffer.push_back(particle.getU());
        sendBuffer.push_back(particle.getV());
        sendBuffer.push_back(static_cast<RealType>(particle.getI()));
        sendBuffer.push_back(static_cast<RealType>(particle.getJ()));
      }
      else{
        sendBuffer.push_back(particle.getX());
        sendBuffer.push_back(particle.getY());
        sendBuffer.push_back(particle.getZ());
        sendBuffer.push_back(particle.getU());
        sendBuffer.push_back(particle.getV());
        sendBuffer.push_back(particle.getW());
        sendBuffer.push_back(static_cast<RealType>(particle.getI()));
        sendBuffer.push_back(static_cast<RealType>(particle.getJ()));
        sendBuffer.push_back(static_cast<RealType>(particle.getK()));
      }

      //TODO: delete the particle
    }
  }

  if(sendBuffer.size() == 0) sendBuffer.push_back(0.0);
  return sendBuffer;
}

std::vector<RealType> ParticleSimulation::collectFrontBoundaryParticles(){
  const int dim = parameters_.geometry.dim;
  std::vector<RealType> sendBuffer{};

  for(auto& particle: particles_){
    if(particle.getK() < 2){ //front ghost cell territory
      if(dim == 2){
        sendBuffer.push_back(particle.getX());
        sendBuffer.push_back(particle.getY());
        sendBuffer.push_back(particle.getU());
        sendBuffer.push_back(particle.getV());
        sendBuffer.push_back(static_cast<RealType>(particle.getI()));
        sendBuffer.push_back(static_cast<RealType>(particle.getJ()));
      }
      else{
        sendBuffer.push_back(particle.getX());
        sendBuffer.push_back(particle.getY());
        sendBuffer.push_back(particle.getZ());
        sendBuffer.push_back(particle.getU());
        sendBuffer.push_back(particle.getV());
        sendBuffer.push_back(particle.getW());
        sendBuffer.push_back(static_cast<RealType>(particle.getI()));
        sendBuffer.push_back(static_cast<RealType>(particle.getJ()));
        sendBuffer.push_back(static_cast<RealType>(particle.getK()));
      }

      //TODO: delete the particle
    }
  }

  if(sendBuffer.size() == 0) sendBuffer.push_back(0.0);
  return sendBuffer;
}

std::vector<RealType> ParticleSimulation::collectBackBoundaryParticles(){
  const int dim = parameters_.geometry.dim;
  std::vector<RealType> sendBuffer{};

  for(auto& particle: particles_){
    if(particle.getK() >= parameters_.parallel.localSize[2] + 2){ //back ghost cell territory
      if(dim == 2){
        sendBuffer.push_back(particle.getX());
        sendBuffer.push_back(particle.getY());
        sendBuffer.push_back(particle.getU());
        sendBuffer.push_back(particle.getV());
        sendBuffer.push_back(static_cast<RealType>(particle.getI()));
        sendBuffer.push_back(static_cast<RealType>(particle.getJ()));
      }
      else{
        sendBuffer.push_back(particle.getX());
        sendBuffer.push_back(particle.getY());
        sendBuffer.push_back(particle.getZ());
        sendBuffer.push_back(particle.getU());
        sendBuffer.push_back(particle.getV());
        sendBuffer.push_back(particle.getW());
        sendBuffer.push_back(static_cast<RealType>(particle.getI()));
        sendBuffer.push_back(static_cast<RealType>(particle.getJ()));
        sendBuffer.push_back(static_cast<RealType>(particle.getK()));
      }

      //TODO: delete the particle
    }
  }

  if(sendBuffer.size() == 0) sendBuffer.push_back(0.0);
  return sendBuffer;
}

void ParticleSimulation::communicateParticles(){
  const int dim = parameters_.geometry.dim;

  //Collect buffers
  std::vector<RealType> leftSendBuffer;
  std::vector<RealType> rightSendBuffer;
  std::vector<RealType> bottomSendBuffer;
  std::vector<RealType> topSendBuffer;
  std::vector<RealType> frontSendBuffer;
  std::vector<RealType> backSendBuffer;

  int leftSendCount;
  int rightSendCount;
  int bottomSendCount;
  int topSendCount;
  int frontSendCount;
  int backSendCount;

  int leftRecvCount = 0;
  int rightRecvCount = 0;
  int bottomRecvCount = 0;
  int topRecvCount = 0;
  int frontRecvCount = 0;
  int backRecvCount = 0;

  if(dim == 2){
    leftSendBuffer = collectLeftBoundaryParticles();
    rightSendBuffer = collectRightBoundaryParticles();
    bottomSendBuffer = collectBottomBoundaryParticles();
    topSendBuffer = collectTopBoundaryParticles();

    leftSendCount = leftSendBuffer.size();
    rightSendCount = rightSendBuffer.size();
    bottomSendCount = bottomSendBuffer.size();
    topSendCount = topSendBuffer.size();
  }
  else{
    leftSendBuffer = collectLeftBoundaryParticles();
    rightSendBuffer = collectRightBoundaryParticles();
    bottomSendBuffer = collectBottomBoundaryParticles();
    topSendBuffer = collectTopBoundaryParticles();
    frontSendBuffer = collectFrontBoundaryParticles();
    backSendBuffer = collectBackBoundaryParticles();

    leftSendCount = leftSendBuffer.size();
    rightSendCount = rightSendBuffer.size();
    bottomSendCount = bottomSendBuffer.size();
    topSendCount = topSendBuffer.size();
    frontSendCount = frontSendBuffer.size();
    backSendCount = backSendBuffer.size();
  }

  //MPI communication of number of particles
  if(dim == 2){
    // send from left, receive on right
    MPI_Sendrecv(
      &leftSendCount,
      1,
      MPI_INT,
      parameters_.parallel.leftNb,
      0,
      &rightRecvCount,
      1,
      MPI_INT,
      parameters_.parallel.rightNb,
      0,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from right, receive on left
    MPI_Sendrecv(
      &rightSendCount,
      1,
      MPI_INT,
      parameters_.parallel.rightNb,
      1,
      &leftRecvCount,
      1,
      MPI_INT,
      parameters_.parallel.leftNb,
      1,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from bottom, receive on top
    MPI_Sendrecv(
      &bottomSendCount,
      1,
      MPI_INT,
      parameters_.parallel.bottomNb,
      2,
      &topRecvCount,
      1,
      MPI_INT,
      parameters_.parallel.topNb,
      2,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from top, receive on bottom
    MPI_Sendrecv(
      &topSendCount,
      1,
      MPI_INT,
      parameters_.parallel.topNb,
      3,
      &bottomRecvCount,
      1,
      MPI_INT,
      parameters_.parallel.bottomNb,
      3,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

  }
  else{
    // send from left, receive on right
    MPI_Sendrecv(
      &leftSendCount,
      1,
      MPI_INT,
      parameters_.parallel.leftNb,
      0,
      &rightRecvCount,
      1,
      MPI_INT,
      parameters_.parallel.rightNb,
      0,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from right, receive on left
    MPI_Sendrecv(
      &rightSendCount,
      1,
      MPI_INT,
      parameters_.parallel.rightNb,
      1,
      &leftRecvCount,
      1,
      MPI_INT,
      parameters_.parallel.leftNb,
      1,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from bottom, receive on top
    MPI_Sendrecv(
      &bottomSendCount,
      1,
      MPI_INT,
      parameters_.parallel.bottomNb,
      2,
      &topRecvCount,
      1,
      MPI_INT,
      parameters_.parallel.topNb,
      2,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from top, receive on bottom
    MPI_Sendrecv(
      &topSendCount,
      1,
      MPI_INT,
      parameters_.parallel.topNb,
      3,
      &bottomRecvCount,
      1,
      MPI_INT,
      parameters_.parallel.bottomNb,
      3,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from front, receive on back
    MPI_Sendrecv(
      &frontSendCount,
      1,
      MPI_INT,
      parameters_.parallel.frontNb,
      4,
      &backRecvCount,
      1,
      MPI_INT,
      parameters_.parallel.backNb,
      4,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from back, receive on front
    MPI_Sendrecv(
      &backSendCount,
      1,
      MPI_INT,
      parameters_.parallel.backNb,
      5,
      &frontRecvCount,
      1,
      MPI_INT,
      parameters_.parallel.frontNb,
      5,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );
  }

  //MPI communication of buffers
  std::vector<RealType> leftRecvBuffer;
  leftRecvBuffer.reserve(leftRecvCount);

  std::vector<RealType> rightRecvBuffer;
  rightRecvBuffer.reserve(rightRecvCount);

  std::vector<RealType> bottomRecvBuffer;
  bottomRecvBuffer.reserve(bottomRecvCount);

  std::vector<RealType> topRecvBuffer;
  topRecvBuffer.reserve(topRecvCount);

  std::vector<RealType> frontRecvBuffer;
  frontRecvBuffer.reserve(frontRecvCount);

  std::vector<RealType> backRecvBuffer;
  backRecvBuffer.reserve(backRecvCount);

  if(dim == 2){
    // send from left, receive on right
    MPI_Sendrecv(
      leftSendBuffer.data(),
      leftSendCount,
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      0,
      rightRecvBuffer.data(),
      rightRecvCount,
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      0,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from right, receive on left
    MPI_Sendrecv(
      rightSendBuffer.data(),
      rightSendCount,
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      1,
      leftRecvBuffer.data(),
      leftRecvCount,
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      1,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from bottom, receive on top
    MPI_Sendrecv(
      bottomSendBuffer.data(),
      bottomSendCount,
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      2,
      topRecvBuffer.data(),
      topRecvCount,
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      2,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from top, receive on bottom
    MPI_Sendrecv(
      topSendBuffer.data(),
      topSendCount,
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      3,
      bottomRecvBuffer.data(),
      bottomRecvCount,
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      3,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

  }
  else{
    // send from left, receive on right
    MPI_Sendrecv(
      leftSendBuffer.data(),
      leftSendCount,
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      0,
      rightRecvBuffer.data(),
      rightRecvCount,
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      0,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from right, receive on left
    MPI_Sendrecv(
      rightSendBuffer.data(),
      rightSendCount,
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      1,
      leftRecvBuffer.data(),
      leftRecvCount,
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      1,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from bottom, receive on top
    MPI_Sendrecv(
      bottomSendBuffer.data(),
      bottomSendCount,
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      2,
      topRecvBuffer.data(),
      topRecvCount,
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      2,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from top, receive on bottom
    MPI_Sendrecv(
      topSendBuffer.data(),
      topSendCount,
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      3,
      bottomRecvBuffer.data(),
      bottomRecvCount,
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      3,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from front, receive on back
    MPI_Sendrecv(
      frontSendBuffer.data(),
      frontSendCount,
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      4,
      backRecvBuffer.data(),
      backRecvCount,
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      4,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from back, receive on front
    MPI_Sendrecv(
      backSendBuffer.data(),
      backSendCount,
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      5,
      frontRecvBuffer.data(),
      frontRecvCount,
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      5,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );
  }

  //Particle generation
  if(dim == 2){
    if(leftRecvCount >= 6){
      for(int i; i<leftRecvCount/6; i++){
        Particle particle(&leftRecvBuffer[i*6], flowField_, parameters_);
        particle.getI() = 2;
        particles_.push_back(particle);
      }
    }
    if(rightRecvCount >= 6){
      for(int i; i<rightRecvCount/6; i++){
        Particle particle(&rightRecvBuffer[i*6], flowField_, parameters_);
        particle.getI() = parameters_.parallel.localSize[0] + 1;
        particles_.push_back(particle);
      }
    }
    if(bottomRecvCount >= 6){
      for(int i; i<bottomRecvCount/6; i++){
        Particle particle(&bottomRecvBuffer[i*6], flowField_, parameters_);
        particle.getJ() = 2;
        particles_.push_back(particle);
      }
    }
    if(topRecvCount >= 6){
      for(int i; i<topRecvCount/6; i++){
        Particle particle(&topRecvBuffer[i*6], flowField_, parameters_);
        particle.getJ() = parameters_.parallel.localSize[1] + 1;
        particles_.push_back(particle);
      }
    }
  }
  else{
    //TODO: Initialize the particles in 3D case
  }
}