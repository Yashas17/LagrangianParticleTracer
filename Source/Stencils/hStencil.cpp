#include "StdAfx.hpp"

#include "hStencil.hpp"

Stencils::hStencil::hStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters),
  xLimit_(parameters.bfStep.xRatio * parameters.geometry.lengthX),
  yLimit_(parameters.bfStep.yRatio * parameters.geometry.lengthY) {}

void Stencils::hStencil::apply(FlowField& flowField, int i, int j) {
  RealType&      h   = flowField.getH().getScalar(i, j);
  const RealType dx2 = parameters_.meshsize->getDx(i, j) / 2;
  const RealType dy2 = parameters_.meshsize->getDy(i, j) / 2;

  const RealType posX = parameters_.meshsize->getPosX(i, j) + dx2;
  const RealType posY = parameters_.meshsize->getPosY(i, j) + dy2;

  const RealType lengthX = parameters_.geometry.lengthX;
  const RealType lengthY = parameters_.geometry.lengthY;

  const RealType stepLengthX = parameters_.bfStep.xRatio * lengthX + dx2;
  const RealType stepLengthY = parameters_.bfStep.yRatio * lengthY + dy2;

  // Check if this is the BFS case
  if (parameters_.bfStep.xRatio > 0 && parameters_.bfStep.yRatio > 0) {
    // Case 1: the step is a left wall to the current cell
    if (posY <= stepLengthY && posX > stepLengthX)
      h = std::min(posX - stepLengthX, std::min(posY, lengthY - posY));

    // Case 2: the cell is above the step
    else if (posX < stepLengthX) {
      h = std::min(posY - stepLengthY, lengthY - posY);
    }

    // Case 3: the cell is above the level of the step but not above the step
    else {
      const RealType distanceToStepX      = posX - stepLengthX;
      const RealType distanceToStepY      = posY - stepLengthY;
      const RealType distanceToStepCorner = std::hypot(distanceToStepX, distanceToStepY);

      h = std::min(distanceToStepCorner, std::min(posY, lengthY - posY));
    }

  } else {
    // channel flow case
    h = std::min(lengthY - posY, posY);
  }
}

void Stencils::hStencil::apply(FlowField& flowField, int i, int j, int k) {
  RealType&      h   = flowField.getH().getScalar(i, j, k);
  const RealType dx2 = parameters_.meshsize->getDx(i, j, k) / 2;
  const RealType dy2 = parameters_.meshsize->getDy(i, j, k) / 2;
  const RealType dz2 = parameters_.meshsize->getDz(i, j, k) / 2;

  const RealType posX = parameters_.meshsize->getPosX(i, j, k) + dx2;
  const RealType posY = parameters_.meshsize->getPosY(i, j, k) + dy2;
  const RealType posZ = parameters_.meshsize->getPosZ(i, j, k) + dz2;

  const RealType lengthX = parameters_.geometry.lengthX;
  const RealType lengthY = parameters_.geometry.lengthY;
  const RealType lengthZ = parameters_.geometry.lengthZ;

  const RealType stepLengthX = parameters_.bfStep.xRatio * lengthX + dx2;
  const RealType stepLengthY = parameters_.bfStep.yRatio * lengthY + dy2;

  // Check if this is the BFS case
  if (parameters_.bfStep.xRatio > 0 && parameters_.bfStep.yRatio > 0) {
    // Case 1: the step is a left wall to the current cell
    if (posY <= stepLengthY && posX > stepLengthX)
      h = std::min(posX - stepLengthX, std::min(posY, lengthY - posY));

    // Case 2: the cell is above the step
    else if (posX < stepLengthX)
      h = std::min(posY - stepLengthY, lengthY - posY);

    // Case 3: the cell is above the level of the step but not above the step itself
    else {
      const RealType distanceToStepX      = posX - stepLengthX;
      const RealType distanceToStepY      = posY - stepLengthY;
      const RealType distanceToStepCorner = std::hypot(distanceToStepX, distanceToStepY);

      h = std::min(distanceToStepCorner, std::min(posY, lengthY - posY));
    }

  } else {
    // channel flow case
    h = std::min(lengthY - posY, posY);
  }

  // Considering the distance from the front and back walls (z direction)
  h = std::min(h, std::min(posZ, lengthZ - posZ));
}
