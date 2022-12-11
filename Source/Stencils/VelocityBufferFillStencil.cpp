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
    // left_buffer[2 * (j - 2)]     = flowField.getVelocity().getVector(2, j)[0];
    // left_buffer[2 * (j - 2) + 1] = flowField.getVelocity().getVector(2, j)[1];
    std::copy_n(flowField.getVelocity().getVector(2, j), 2, &left_buffer[2 * (j - 2)]);
  }
#ifndef NDEBUG
  std::stringstream ss;
  ss << parameters_.parallel.rank << " is sending to its left neighbor (velocity): \n";
  for (auto& el : left_buffer) {
    ss << el << ", ";
  }
  ss << "end";
  ss << std::endl;
  std::cout << ss.str();
#endif
}

void Stencils::VelocityBufferFillStencil::applyRightWall2D(FlowField& flowField) {
  right_buffer.resize((flowField.getVelocity().getNy() - 3) * 4);
  for (int j = 2; j < flowField.getVelocity().getNy() - 1; j++) {
    // Equivalent code below
    // right_buffer[4 * (j - 2)]     = flowField.getVelocity().getVector(flowField.getVelocity().getNx() - 3, j)[0];
    // right_buffer[4 * (j - 2) + 1] = flowField.getVelocity().getVector(flowField.getVelocity().getNx() - 3, j)[1];
    // right_buffer[4 * (j - 2) + 2] = flowField.getVelocity().getVector(flowField.getVelocity().getNx() - 2, j)[0];
    // right_buffer[4 * (j - 2) + 3] = flowField.getVelocity().getVector(flowField.getVelocity().getNx() - 2, j)[1];
    std::copy_n(
      flowField.getVelocity().getVector(flowField.getVelocity().getNx() - 3, j), 4, &right_buffer[4 * (j - 2)]
    );
  }
#ifndef NDEBUG
  std::stringstream ss;
  ss << parameters_.parallel.rank << " is sending to its right neighbor (velocity): \n";
  for (auto& el : right_buffer) {
    ss << el << ", ";
  }
  ss << "end";
  ss << std::endl;
  std::cout << ss.str();
#endif
}

void Stencils::VelocityBufferFillStencil::applyBottomWall2D(FlowField& flowField) {
  bottom_buffer.resize((flowField.getVelocity().getNx()) * 2);
  RealType* bottom_line_begin = flowField.getVelocity().getVector(0, 2);
  std::copy_n(bottom_line_begin, (flowField.getVelocity().getNx()) * 2, bottom_buffer.data());
}

void Stencils::VelocityBufferFillStencil::applyTopWall2D(FlowField& flowField) {
  // We first write the outer line uvuvuvuv
  // And then write the v from the inner line
  //
  // Copy the first line
  top_buffer.resize((flowField.getVelocity().getNx()) * 4);
  RealType* outer_top_line_begin = flowField.getVelocity().getVector(0, flowField.getVelocity().getNy() - 3);
  std::copy_n(outer_top_line_begin, (flowField.getVelocity().getNx()) * 4, top_buffer.data());

  // Copy the second, for easier access we save the pointer of the offset
  /*
  RealType* offset_for_inner_line = top_buffer.data() + (flowField.getVelocity().getNx()) * 2;
  for (int i = 0; i < flowField.getVelocity().getNx(); i++) {
    offset_for_inner_line[i] = flowField.getVelocity().getVector(i, flowField.getVelocity().getNy() - 3)[1];
  }
  */
}

void Stencils::VelocityBufferFillStencil::applyLeftWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyRightWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyBottomWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyTopWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyFrontWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyBackWall3D(FlowField& flowField) {}