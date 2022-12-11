#pragma once
#include "StdAfx.hpp"

#include "PressureBufferReadStencil.hpp"

#include "Parameters.hpp"

#include "../FlowField.hpp"

Stencils::PressureBufferReadStencil::PressureBufferReadStencil(const Parameters& parameters):
  GhostLayerStencil(parameters),
  parameters_(parameters) {}

Stencils::PressureBufferReadStencil::~PressureBufferReadStencil() {}

void Stencils::PressureBufferReadStencil::applyLeftWall2D(FlowField& flowField) {
#ifndef NDEBUG
  std::stringstream ss;
  ss << parameters_.parallel.rank << " received from left neighbor (pressure): \n";
  for (auto& el : left_buffer) {
    ss << el << ", ";
  }
  ss << "end";
  ss << std::endl;
  std::cout << ss.str();
#endif
  assert(left_buffer.size() == flowField.getPressure().getNy() - 3);
  for (size_t j = 2; j < flowField.getPressure().getNy() - 1; j++) {
    flowField.getPressure().getScalar(1, j) = left_buffer[j - 2];
  }
}

void Stencils::PressureBufferReadStencil::applyRightWall2D(FlowField& flowField) {
#ifndef NDEBUG
  std::stringstream ss;
  ss << parameters_.parallel.rank << " received from right neighbor (pressure): \n";
  for (auto& el : right_buffer) {
    ss << el << ", ";
  }
  ss << "end";
  ss << std::endl;
  std::cout << ss.str();
#endif
  assert(right_buffer.size() == flowField.getPressure().getNy() - 3);
  for (size_t j = 2; j < flowField.getPressure().getNy() - 1; j++) {
    flowField.getPressure().getScalar(flowField.getPressure().getNx() - 1, j) = right_buffer[j - 2];
  }
}

void Stencils::PressureBufferReadStencil::applyBottomWall2D(FlowField& flowField) {
  assert(bottom_buffer.size() == flowField.getPressure().getNx() - 3);
  RealType* bottom_ghost_line_begin = flowField.getPressure().getScalarOffset(2, 1);
  std::copy(bottom_buffer.begin(), bottom_buffer.end(), bottom_ghost_line_begin);
}

void Stencils::PressureBufferReadStencil::applyTopWall2D(FlowField& flowField) {
  assert(top_buffer.size() == flowField.getPressure().getNx() - 3);
  RealType* top_ghost_line_begin = flowField.getPressure().getScalarOffset(2, flowField.getPressure().getNy() - 1);
  std::copy(top_buffer.begin(), top_buffer.end(), top_ghost_line_begin);
}

void Stencils::PressureBufferReadStencil::applyLeftWall3D(FlowField& flowField) {}

void Stencils::PressureBufferReadStencil::applyRightWall3D(FlowField& flowField) {}

void Stencils::PressureBufferReadStencil::applyBottomWall3D(FlowField& flowField) {}

void Stencils::PressureBufferReadStencil::applyTopWall3D(FlowField& flowField) {}

void Stencils::PressureBufferReadStencil::applyFrontWall3D(FlowField& flowField) {}

void Stencils::PressureBufferReadStencil::applyBackWall3D(FlowField& flowField) {}