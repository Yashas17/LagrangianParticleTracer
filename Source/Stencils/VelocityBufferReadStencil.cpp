#include "StdAfx.hpp"

#include "VelocityBufferReadStencil.hpp"

#include <format>
#include <string_view>

#include "Parameters.hpp"

#include "../FlowField.hpp"

Stencils::VelocityBufferReadStencil::VelocityBufferReadStencil(const Parameters& parameters):
  GhostLayerStencil(parameters),
  parameters_(parameters) {}

Stencils::VelocityBufferReadStencil::~VelocityBufferReadStencil() {}

void Stencils::VelocityBufferReadStencil::applyLeftWall2D(FlowField& flowField) {
  assert(left_buffer.size() == (flowField.getVelocity().getNy() - 3) * 4);
  for (int j = 2; j < flowField.getVelocity().getNy() - 1; j++) {
    std::copy_n(&left_buffer[4 * (j - 2)], 4, flowField.getVelocity().getVector(0, j));
  }
}

void Stencils::VelocityBufferReadStencil::applyRightWall2D(FlowField& flowField) {
  for (int j = 2; j < flowField.getVelocity().getNy() - 1; j++) {
    std::copy_n(
      &right_buffer[2 * (j - 2)], 2, flowField.getVelocity().getVector(flowField.getVelocity().getNx() - 1, j)
    );
  }
}

void Stencils::VelocityBufferReadStencil::applyBottomWall2D(FlowField& flowField) {
  RealType* bottom_line_begin = flowField.getVelocity().getVector(0, 0);
  std::copy_n(bottom_buffer.begin(), (flowField.getVelocity().getNx()) * 4, bottom_line_begin);
}

void Stencils::VelocityBufferReadStencil::applyTopWall2D(FlowField& flowField) {
  RealType* top_line_begin = flowField.getVelocity().getVector(0, flowField.getVelocity().getNy() - 1);
  std::copy_n(top_buffer.begin(), (flowField.getVelocity().getNx()) * 2, top_line_begin);
}

void Stencils::VelocityBufferReadStencil::applyLeftWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyRightWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyBottomWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyTopWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyFrontWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyBackWall3D(FlowField& flowField) {}