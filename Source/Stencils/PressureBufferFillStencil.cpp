#include "StdAfx.hpp"

#include "PressureBufferFillStencil.hpp"

#include "Parameters.hpp"

#include "../FlowField.hpp"

Stencils::PressureBufferFillStencil::PressureBufferFillStencil(const Parameters& parameters):
  GhostLayerStencil(parameters),
  parameters_(parameters) {}

Stencils::PressureBufferFillStencil::~PressureBufferFillStencil() {}

void Stencils::PressureBufferFillStencil::applyLeftWall2D(FlowField& flowField) {
  left_buffer.resize(flowField.getPressure().getNy());
  for (int j = 0; j < flowField.getPressure().getNy(); j++) {
    left_buffer[j] = flowField.getPressure().getScalar(2, j);
  }
#ifndef NDEBUG
  std::stringstream ss;
  ss << parameters_.parallel.rank << "sending to left neighbor (pressure): \n";
  for (auto& el : left_buffer) {
    ss << el << ", ";
  }
  ss << "end";
  ss << std::endl;
  std::cout << ss.str();
#endif
}

void Stencils::PressureBufferFillStencil::applyRightWall2D(FlowField& flowField) {
  right_buffer.resize(flowField.getPressure().getNy());
  for (int j = 0; j < flowField.getPressure().getNy(); j++) {
    right_buffer[j] = flowField.getPressure().getScalar(flowField.getPressure().getNx() - 2, j);
  }
#ifndef NDEBUG
  std::stringstream ss;
  ss << parameters_.parallel.rank << "sending to right neighbor (pressure): \n";
  for (auto& el : right_buffer) {
    ss << el << ", ";
  }
  ss << "end";
  ss << std::endl;
  std::cout << ss.str();
#endif
}

void Stencils::PressureBufferFillStencil::applyBottomWall2D(FlowField& flowField) {
  bottom_buffer.resize(flowField.getPressure().getNx());
  RealType* bottom_line_begin = flowField.getPressure().getScalarOffset(0, 2);
  std::copy_n(bottom_line_begin, flowField.getPressure().getNx(), bottom_buffer.data());
}

void Stencils::PressureBufferFillStencil::applyTopWall2D(FlowField& flowField) {
  top_buffer.resize(flowField.getPressure().getNx());
  RealType* top_line_begin = flowField.getPressure().getScalarOffset(0, flowField.getPressure().getNy() - 2);
  std::copy_n(top_line_begin, flowField.getPressure().getNx(), top_buffer.data());
}

void Stencils::PressureBufferFillStencil::applyLeftWall3D(FlowField& flowField) {}

void Stencils::PressureBufferFillStencil::applyRightWall3D(FlowField& flowField) {}

void Stencils::PressureBufferFillStencil::applyBottomWall3D(FlowField& flowField) {}

void Stencils::PressureBufferFillStencil::applyTopWall3D(FlowField& flowField) {}

void Stencils::PressureBufferFillStencil::applyFrontWall3D(FlowField& flowField) {}

void Stencils::PressureBufferFillStencil::applyBackWall3D(FlowField& flowField) {}