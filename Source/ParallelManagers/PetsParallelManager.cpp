#include "StdAfx.hpp"

#include "PetscParallelManager.hpp"

ParallelManagers::PetscParallelManagerNonBlocking::PetscParallelManagerNonBlocking(
  Parameters& parameters, FlowField& flowField
):
  parameters_(parameters),
  flowField_(flowField),
  pressureBufferFillStencil_(parameters),
  pressureBufferReadStencil_(parameters),
  pressureBufferFillIterator_(flowField, parameters, pressureBufferFillStencil_),
  pressureBufferReadIterator_(flowField, parameters, pressureBufferReadStencil_),
  velocityBufferFillStencil_(parameters),
  velocityBufferReadStencil_(parameters),
  velocityBufferFillIterator_(flowField, parameters, velocityBufferFillStencil_),
  velocityBufferReadIterator_(flowField, parameters, velocityBufferReadStencil_) {}

ParallelManagers::PetscParallelManagerNonBlocking::~PetscParallelManagerNonBlocking() {}

void ParallelManagers::PetscParallelManagerNonBlocking::communicatePressure() {}

void ParallelManagers::PetscParallelManagerNonBlocking::communicateVelocities() {}
