#include "StdAfx.hpp"

#include "PressureBufferFillStencil.hpp"

#include "Parameters.hpp"

#include "../FlowField.hpp"

Stencils::PressureBufferFillStencil::PressureBufferFillStencil(const Parameters& parameters):
  GhostLayerStencil(parameters),
  parameters_(parameters) {}

Stencils::PressureBufferFillStencil::~PressureBufferFillStencil() {}

void Stencils::PressureBufferFillStencil::applyLeftWall2D(FlowField& flowField) {
  left_buffer.resize(flowField.getPressure().getNy() - 3);
  for (int j = 2; j < flowField.getPressure().getNy() - 1; j++) {
    left_buffer[j - 2] = flowField.getPressure().getScalar(2, j);
  }
}

void Stencils::PressureBufferFillStencil::applyRightWall2D(FlowField& flowField) {
  right_buffer.resize(flowField.getPressure().getNy() - 3);
  for (int j = 2; j < flowField.getPressure().getNy() - 1; j++) {
    right_buffer[j - 2] = flowField.getPressure().getScalar(flowField.getPressure().getNx() - 2, j);
  }
}

void Stencils::PressureBufferFillStencil::applyBottomWall2D(FlowField& flowField) {
  bottom_buffer.resize(flowField.getPressure().getNx() - 1);
  RealType* bottom_line_begin = flowField.getPressure().getScalarOffset(1, 2);
  std::copy_n(bottom_line_begin, flowField.getPressure().getNx() - 1, bottom_buffer.data());
}

void Stencils::PressureBufferFillStencil::applyTopWall2D(FlowField& flowField) {
  top_buffer.resize(flowField.getPressure().getNx() - 1);
  RealType* top_line_begin = flowField.getPressure().getScalarOffset(1, flowField.getPressure().getNy() - 2);
  std::copy_n(top_line_begin, flowField.getPressure().getNx() - 1, top_buffer.data());
}

void Stencils::PressureBufferFillStencil::applyLeftWall3D(FlowField& flowField) {
  left_buffer.resize((flowField.getPressure().getNy() - 3) * (flowField.getPressure().getNz() - 3));
  for (int j = 2; j < flowField.getPressure().getNy() - 1; j++) {
    for (int k = 2; k < flowField.getPressure().getNz() - 1; j++) {
      left_buffer[(j - 2) * (k - 2)] = flowField.getPressure().getScalar(2, j, k);
    }
  }
}

void Stencils::PressureBufferFillStencil::applyRightWall3D(FlowField& flowField) {
  right_buffer.resize((flowField.getPressure().getNy() - 3) * (flowField.getPressure().getNz() - 3));
  for (int j = 2; j < flowField.getPressure().getNy() - 1; j++) {
    for (int k = 2; k < flowField.getPressure().getNz() - 1; j++) {
      right_buffer[(j - 2) * (k - 2)] = flowField.getPressure().getScalar(flowField.getPressure().getNx() - 2, j, k);
    }
  }
}

void Stencils::PressureBufferFillStencil::applyBottomWall3D(FlowField& flowField) {
  bottom_buffer.resize((flowField.getPressure().getNx() - 1) * (flowField.getPressure().getNz() - 3));
  for (int k = 2; k < flowField.getPressure().getNz() - 1; k++) {
    std::copy_n(
      flowField.getPressure().getScalarOffset(1, 2, k),
      flowField.getPressure().getNx() - 1,
      bottom_buffer.data() + ((k - 2) * (flowField.getPressure().getNx() - 1))
    );
  }
}

void Stencils::PressureBufferFillStencil::applyTopWall3D(FlowField& flowField) {
  bottom_buffer.resize((flowField.getPressure().getNx() - 1) * (flowField.getPressure().getNz() - 3));
  for (int k = 2; k < flowField.getPressure().getNz() - 1; k++) {
    std::copy_n(
      flowField.getPressure().getScalarOffset(1, flowField.getPressure().getNy() - 2, k),
      flowField.getPressure().getNx() - 1,
      bottom_buffer.data() + ((k - 2) * (flowField.getPressure().getNx() - 1))
    );
  }
}

void Stencils::PressureBufferFillStencil::applyFrontWall3D(FlowField& flowField) {
  // Copy the whole face so we dont deal with access iterators
  front_buffer.resize((flowField.getPressure().getNx()) * (flowField.getPressure().getNy() - 1));
  std::copy_n(
    flowField.getPressure().getScalarOffset(0, 1, 2),
    (flowField.getPressure().getNx()) * (flowField.getPressure().getNy() - 1),
    front_buffer.data()
  );
}

void Stencils::PressureBufferFillStencil::applyBackWall3D(FlowField& flowField) {
  back_buffer.resize((flowField.getPressure().getNx()) * (flowField.getPressure().getNy() - 1));
  std::copy_n(
    flowField.getPressure().getScalarOffset(0, 1, flowField.getPressure().getNz() - 2),
    (flowField.getPressure().getNx()) * (flowField.getPressure().getNy() - 1),
    back_buffer.data()
  );
}