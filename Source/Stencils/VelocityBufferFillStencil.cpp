#include "StdAfx.hpp"

#include "VelocityBufferFillStencil.hpp"

#include "Parameters.hpp"

#include "../FlowField.hpp"

Stencils::VelocityBufferFillStencil::VelocityBufferFillStencil(const Parameters& parameters):
  GhostLayerStencil(parameters),
  parameters_(parameters) {}

Stencils::VelocityBufferFillStencil::~VelocityBufferFillStencil() {}

void Stencils::VelocityBufferFillStencil::applyLeftWall2D(FlowField& flowField) {
  left_buffer.resize((flowField.getVelocity().getNy() - 3) * 2);
  for (int j = 2; j < flowField.getVelocity().getNy() - 1; j++) {
    std::copy_n(flowField.getVelocity().getVector(2, j), 2, &left_buffer[2 * (j - 2)]);
  }
}

void Stencils::VelocityBufferFillStencil::applyRightWall2D(FlowField& flowField) {
  right_buffer.resize((flowField.getVelocity().getNy() - 3) * 4);
  for (int j = 2; j < flowField.getVelocity().getNy() - 1; j++) {
    std::copy_n(
      flowField.getVelocity().getVector(flowField.getVelocity().getNx() - 3, j), 4, &right_buffer[4 * (j - 2)]
    );
  }
}

void Stencils::VelocityBufferFillStencil::applyBottomWall2D(FlowField& flowField) {
  bottom_buffer.resize((flowField.getVelocity().getNx()) * 2);
  RealType* bottom_line_begin = flowField.getVelocity().getVector(0, 2);
  std::copy_n(bottom_line_begin, (flowField.getVelocity().getNx()) * 2, bottom_buffer.data());
}

void Stencils::VelocityBufferFillStencil::applyTopWall2D(FlowField& flowField) {
  top_buffer.resize((flowField.getVelocity().getNx()) * 4);
  RealType* outer_top_line_begin = flowField.getVelocity().getVector(0, flowField.getVelocity().getNy() - 3);
  std::copy_n(outer_top_line_begin, (flowField.getVelocity().getNx()) * 4, top_buffer.data());
}

void Stencils::VelocityBufferFillStencil::applyLeftWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyRightWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyBottomWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyTopWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyFrontWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyBackWall3D(FlowField& flowField) {}