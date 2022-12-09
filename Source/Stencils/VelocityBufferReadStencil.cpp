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
#ifndef NDEBUG
  std::stringstream ss;
  ss << parameters_.parallel.rank << " received from left neighbor (velocity): \n";
  for (auto& el : left_buffer) {
    ss << el << ", ";
  }
  ss << "end";
  ss << std::endl;
  std::cout << ss.str();
#endif

#ifndef NDEBUG
  std::stringstream ss2;
  ss2 << left_buffer.size() << " == " << ((flowField.getVelocity().getNy() - 3) * 3);
  ss2 << " doesn't hold.";
  if (left_buffer.size() != ((flowField.getVelocity().getNy() - 3) * 3)) {
    throw std::runtime_error(ss2.str());
  }
#endif
  assert(left_buffer.size() == (flowField.getVelocity().getNy() - 3) * 3);
  for (size_t j = 2; j < flowField.getVelocity().getNy() - 1; j++) {
    flowField.getVelocity().getVector(0, j)[0] = left_buffer[3 * (j - 2)];
    flowField.getVelocity().getVector(1, j)[0] = left_buffer[3 * (j - 2) + 1];
    flowField.getVelocity().getVector(1, j)[1] = left_buffer[3 * (j - 2) + 2];
  }
}

void Stencils::VelocityBufferReadStencil::applyRightWall2D(FlowField& flowField) {
#ifndef NDEBUG
  std::stringstream ss;
  ss << parameters_.parallel.rank << " received from right neighbor (velocity): \n";
  for (auto& el : right_buffer) {
    ss << el << ", ";
  }
  ss << "end";
  ss << std::endl;
  std::cout << ss.str();
#endif

#ifndef NDEBUG
  std::stringstream ss2;
  ss2 << right_buffer.size() << " == " << ((flowField.getVelocity().getNy() - 3) * 2);
  ss2 << " doesn't hold.";
  if (right_buffer.size() != (flowField.getVelocity().getNy() - 3) * 2) {
    throw std::runtime_error(ss2.str());
  }
#endif
  for (size_t j = 2; j < flowField.getVelocity().getNy() - 1; j++) {
    flowField.getVelocity().getVector(flowField.getVelocity().getNx() - 1, j)[0] = right_buffer[2 * (j - 2)];
    flowField.getVelocity().getVector(flowField.getVelocity().getNx() - 1, j)[1] = right_buffer[2 * (j - 2) + 1];
  }
}

void Stencils::VelocityBufferReadStencil::applyBottomWall2D(FlowField& flowField) {
#ifndef NDEBUG
  std::stringstream ss;
  ss << bottom_buffer.size() << " == " << ((flowField.getVelocity().getNx() - 3) * 3);
  ss << " doesn't hold.";
  if (bottom_buffer.size() != (flowField.getVelocity().getNx() - 3) * 3) {
    throw std::runtime_error(ss.str());
  }
#endif
  // The outer line of the received buffer will be brought to line 2
  RealType* bottom_line_begin = flowField.getVelocity().getVector(2, 1);
  std::copy_n(bottom_buffer.begin(), (flowField.getVelocity().getNx() - 3) * 2, bottom_line_begin);

  RealType* buffer_offset = bottom_buffer.data() + ((flowField.getVelocity().getNx() - 3) * 2);
  for (size_t i = 2; i < flowField.getVelocity().getNx() - 1; i++) {
    flowField.getVelocity().getVector(i, 0)[1] = buffer_offset[i - 2];
  }
}

void Stencils::VelocityBufferReadStencil::applyTopWall2D(FlowField& flowField) {
#ifndef NDEBUG
  std::stringstream ss;
  ss << top_buffer.size() << " == " << ((flowField.getVelocity().getNx() - 3) * 2);
  ss << " doesn't hold.";
  if (top_buffer.size() != (flowField.getVelocity().getNx() - 3) * 2) {
    throw std::runtime_error(ss.str());
  }
#endif
  RealType* top_line_begin = flowField.getVelocity().getVector(2, flowField.getVelocity().getNy() - 1);
  std::copy_n(top_buffer.begin(), (flowField.getVelocity().getNx() - 3) * 2, top_line_begin);
}

void Stencils::VelocityBufferReadStencil::applyLeftWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyRightWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyBottomWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyTopWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyFrontWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyBackWall3D(FlowField& flowField) {}