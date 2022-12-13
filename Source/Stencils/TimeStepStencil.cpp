#include "StdAfx.hpp"

#include "TimeStepStencil.hpp"

Stencils::TimeStepStencil::TimeStepStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::TimeStepStencil::apply(FlowField& flowField, int i, int j) {

  const RealType        vt       = flowField.getVt().getScalar(i, j);
  const RealType* const velocity = flowField.getVelocity().getVector(i, j);

  const RealType vTotal = vt + (1 / parameters_.flow.Re);

  const RealType dx = parameters_.meshsize->getDx(i, j);
  const RealType dy = parameters_.meshsize->getDy(i, j);

  RealType cell_dt = parameters_.timestep.dt;

  ASSERTION(parameters_.geometry.dim == 2);
  RealType factor = 1.0 / (dx * dx) + 1.0 / (dy * dy);

  /**
   *We need the following 'if' condition to ensure that we do not encounter division by zero. In case the 'if'
   *condition is not satisfied, the value of dt will be same as previous time step or if it is the first time step,
   *the user defined time-step will be used. If user has not defined the time step, dt will be set to 1 by default.
   **/

  // If the cell is NOT an obstacle
  if (flowField.getFlags().getValue(i, j) % 2 == 0) {
    if (velocity[0] != 0) {
      cell_dt = dx / fabs(velocity[0]);
    }

    if (velocity[0] != 0 && velocity[1] != 0) {
      cell_dt = std::min(
        1 / (2 * vTotal * factor), std::min(cell_dt, std::min(dx / fabs(velocity[0]), dy / fabs(velocity[1])))
      );

    } else {
      cell_dt = std::min(1 / (2 * vTotal * factor), cell_dt);
    }

    if (cell_dt < dt)
      dt = cell_dt;
  }
}

void Stencils::TimeStepStencil::apply(FlowField& flowField, int i, int j, int k) {
  const RealType        vt       = flowField.getVt().getScalar(i, j, k);
  const RealType* const velocity = flowField.getVelocity().getVector(i, j, k);

  const RealType vTotal = vt + (1 / parameters_.flow.Re);

  const RealType dx = parameters_.meshsize->getDx(i, j, k);
  const RealType dy = parameters_.meshsize->getDy(i, j, k);
  const RealType dz = parameters_.meshsize->getDz(i, j, k);

  RealType cell_dt = parameters_.timestep.dt;

  ASSERTION(parameters_.geometry.dim == 3);
  RealType factor = 1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz);

  /*We need the following 'if' condition to ensure that we do not encounter division by zero. In case the 'if'
   *condition is not satisfied, the value of dt will be same as previous time step or if it is the first time step,
   *the user defined time-step will be used. If user has not defined the time step, dt will be set to 1 by default.
   **/

  if (velocity[2] != 0)
    cell_dt = dz / fabs(velocity[2]);

  if (velocity[0] != 0 && velocity[1] != 0) {
    cell_dt = std::min(
      1 / (2 * vTotal * factor), std::min(cell_dt, std::min(dx / fabs(velocity[0]), dy / fabs(velocity[1])))
    );
  } else {
    cell_dt = std::min(1 / (2 * vTotal * factor), cell_dt);
  }

  if (cell_dt < dt)
    dt = cell_dt;
}