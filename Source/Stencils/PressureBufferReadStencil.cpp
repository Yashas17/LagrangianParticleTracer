#include "StdAfx.hpp"

#include "PressureBufferReadStencil.hpp"

#include "Parameters.hpp"

#include "../FlowField.hpp"

Stencils::PressureBufferReadStencil::PressureBufferReadStencil(const Parameters& parameters):
  GhostLayerStencil(parameters),
  parameters_(parameters) {}

Stencils::PressureBufferReadStencil::~PressureBufferReadStencil() {}

void Stencils::PressureBufferReadStencil::applyLeftWall2D(FlowField& flowField) {
  assert(left_buffer.size() == flowField.getPressure().getNy() - 3);
  for (int j = 2; j < flowField.getPressure().getNy() - 1; j++) {
    flowField.getPressure().getScalar(1, j) = left_buffer[j - 2];
  }
}

void Stencils::PressureBufferReadStencil::applyRightWall2D(FlowField& flowField) {
  assert(right_buffer.size() == flowField.getPressure().getNy() - 3);
  for (int j = 2; j < flowField.getPressure().getNy() - 1; j++) {
    flowField.getPressure().getScalar(flowField.getPressure().getNx() - 1, j) = right_buffer[j - 2];
  }
}

void Stencils::PressureBufferReadStencil::applyBottomWall2D(FlowField& flowField) {
  assert(bottom_buffer.size() == flowField.getPressure().getNx() - 1);
  RealType* bottom_ghost_line_begin = flowField.getPressure().getScalarOffset(1, 1);
  std::copy(bottom_buffer.begin(), bottom_buffer.end(), bottom_ghost_line_begin);
}

void Stencils::PressureBufferReadStencil::applyTopWall2D(FlowField& flowField) {
  assert(top_buffer.size() == flowField.getPressure().getNx() - 1);
  RealType* top_ghost_line_begin = flowField.getPressure().getScalarOffset(1, flowField.getPressure().getNy() - 1);
  std::copy(top_buffer.begin(), top_buffer.end(), top_ghost_line_begin);
}

void Stencils::PressureBufferReadStencil::applyLeftWall3D(FlowField& flowField) {
  assert(left_buffer.size() == (flowField.getPressure().getNy() - 3) * (flowField.getPressure().getNz() - 3));
  for (int j = 2; j < flowField.getPressure().getNy() - 1; j++) {
    for (int k = 2; k < flowField.getPressure().getNz() - 1; j++) {
      flowField.getPressure().getScalar(1, j, k) = left_buffer[(j - 2) * (k - 2)];
    }
  }
}

void Stencils::PressureBufferReadStencil::applyRightWall3D(FlowField& flowField) {
  assert(right_buffer.size() == (flowField.getPressure().getNy() - 3) * (flowField.getPressure().getNz() - 3));
  for (int j = 2; j < flowField.getPressure().getNy() - 1; j++) {
    for (int k = 2; k < flowField.getPressure().getNz() - 1; k++) {
      flowField.getPressure().getScalar(flowField.getPressure().getNx() - 1, j, k) = right_buffer[(j - 2) * (k - 2)];
    }
  }
}

void Stencils::PressureBufferReadStencil::applyBottomWall3D(FlowField& flowField) {
  assert(bottom_buffer.size() == (flowField.getPressure().getNx() - 1) * (flowField.getPressure().getNz() - 3));
  for (int k = 2; k < flowField.getPressure().getNz() - 1; k++) {
    std::copy_n(
      bottom_buffer.data() + ((k - 2) * (flowField.getPressure().getNx() - 1)),
      flowField.getPressure().getNx() - 1,
      flowField.getPressure().getScalarOffset(1, 1, k)
    );
  }
}

void Stencils::PressureBufferReadStencil::applyTopWall3D(FlowField& flowField) {
  assert(top_buffer.size() == (flowField.getPressure().getNx() - 1) * (flowField.getPressure().getNz() - 3));
  for (int k = 2; k < flowField.getPressure().getNz() - 1; k++) {
    std::copy_n(
      top_buffer.data() + ((k - 2) * (flowField.getPressure().getNx() - 1)),
      flowField.getPressure().getNx() - 1,
      flowField.getPressure().getScalarOffset(1, flowField.getPressure().getNy() - 1, k)
    );
  }
}

void Stencils::PressureBufferReadStencil::applyFrontWall3D(FlowField& flowField) {
  // Copy the whole face so we dont deal with access iterators
  assert(front_buffer.size() == (flowField.getPressure().getNx()) * (flowField.getPressure().getNy() - 1));
  std::copy_n(
    front_buffer.data(),
    (flowField.getPressure().getNx()) * (flowField.getPressure().getNy() - 1),
    flowField.getPressure().getScalarOffset(0, 1, 1)
  );
}

void Stencils::PressureBufferReadStencil::applyBackWall3D(FlowField& flowField) {
  assert(back_buffer.size() == (flowField.getPressure().getNx()) * (flowField.getPressure().getNy() - 1));
  std::copy_n(
    back_buffer.data(),
    (flowField.getPressure().getNx()) * (flowField.getPressure().getNy() - 1),
    flowField.getPressure().getScalarOffset(0, 1, flowField.getPressure().getNz() - 1)
  );
}