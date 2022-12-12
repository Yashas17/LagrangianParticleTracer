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

/*
void Stencils::VelocityBufferReadStencil::applyBottomWall2D(FlowField& flowField) {
  assert(bottom_buffer.size() == flowField.getVelocity().getNx() - 1);
  RealType* bottom_ghost_line_begin = flowField.getVelocity().getVectorOffset(1, 1);
  std::copy(bottom_buffer.begin(), bottom_buffer.end(), bottom_ghost_line_begin);
}

void Stencils::VelocityBufferReadStencil::applyTopWall2D(FlowField& flowField) {
  assert(top_buffer.size() == flowField.getVelocity().getNx() - 1);
  RealType* top_ghost_line_begin = flowField.getVelocity().getVectorOffset(1, flowField.getVelocity().getNy() - 1);
  std::copy(top_buffer.begin(), top_buffer.end(), top_ghost_line_begin);
}
*/

void Stencils::VelocityBufferReadStencil::applyLeftWall3D(FlowField& flowField) {
  // assert(left_buffer.size() == (flowField.getVelocity().getNy() - 3) * (flowField.getVelocity().getNz() - 3) * 4);
  // for (int j = 2; j < flowField.getVelocity().getNy() - 1; j++) {
  //   for (int k = 2; k < flowField.getVelocity().getNz() - 1; j++) {
  //     std::copy_n(flowField.getVelocity().getVector(0, j, k), 4, left_buffer[(j - 2) * (k - 2) * 4]);
  //   }
  // }
}

void Stencils::VelocityBufferReadStencil::applyRightWall3D(FlowField& flowField) {
  /*
  assert(right_buffer.size() == (flowField.getVelocity().getNy() - 3) * (flowField.getVelocity().getNz() - 3) * 2);
  for (int j = 2; j < flowField.getVelocity().getNy() - 1; j++) {
    for (int k = 2; k < flowField.getVelocity().getNz() - 1; j++) {
      std::copy_n(
        flowField.getVelocity().getVector(flowField.getVelocity().getNx() - 1, j, k),
        2,
        left_buffer[(j - 2) * (k - 2) * 2]
      );
    }
  }
  */
}

void Stencils::VelocityBufferReadStencil::applyBottomWall3D(FlowField& flowField) {
  /*
  assert(bottom_buffer.size() == (flowField.getVelocity().getNx()) * (flowField.getVelocity().getNz() - 3) * 4);
  for (int k = 2; k < flowField.getVelocity().getNz() - 1; j++) {
    std::copy_n(
      bottom_buffer.data() + ((k - 2) * (flowField.getVelocity().getNx()) * 4),
      flowField.getVelocity().getNx() * 4,
      flowField.getVelocity().getVector(0, 0, k),
    );
  }
  */
}

void Stencils::VelocityBufferReadStencil::applyTopWall3D(FlowField& flowField) {
  assert(top_buffer.size() == (flowField.getVelocity().getNx()) * (flowField.getVelocity().getNz() - 3) * 2);
  /*
  for (int k = 2; k < flowField.getVelocity().getNz() - 1; j++) {
    std::copy_n(
      top_buffer.data() + ((k - 2) * (flowField.getVelocity().getNx() - 1)),
      flowField.getVelocity().getNx() * 2,
      flowField.getVelocity().getVector(0, flowField.getVelocity().getNy() - 1, k)
    );
  }
  */
}

void Stencils::VelocityBufferReadStencil::applyFrontWall3D(FlowField& flowField) {
  // Copy the whole face so we dont deal with access iterators
  /*
  assert(front_buffer.size() == (flowField.getVelocity().getNx()) * (flowField.getVelocity().getNy() - 1) * 4);
  std::copy_n(
    front_buffer.data(),
    (flowField.getVelocity().getNx()) * (flowField.getVelocity().getNy() - 1) * 4,
    flowField.getVelocity().getVector(0, 0, 0),
  );
  */
}

void Stencils::VelocityBufferReadStencil::applyBackWall3D(FlowField& flowField) {
  /*
  assert(back_buffer.size() == (flowField.getVelocity().getNx()) * (flowField.getVelocity().getNy() - 1));
  std::copy_n(
    back_buffer.data(),
    (flowField.getVelocity().getNx()) * (flowField.getVelocity().getNy() - 1) * 2,
    flowField.getVelocity().getVector(0, 0, flowField.getVelocity().getNz() - 1),
  );
  */
}