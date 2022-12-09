#include "StdAfx.hpp"

#include "VelocityBufferReadStencil.hpp"

#include "Parameters.hpp"

#include "../FlowField.hpp"

Stencils::VelocityBufferReadStencil::VelocityBufferReadStencil(const Parameters& parameters):
  GhostLayerStencil(parameters),
  parameters_(parameters) {}

Stencils::VelocityBufferReadStencil::~VelocityBufferReadStencil() {}

void Stencils::VelocityBufferReadStencil::applyLeftWall2D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyRightWall2D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyBottomWall2D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyTopWall2D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyLeftWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyRightWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyBottomWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyTopWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyFrontWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferReadStencil::applyBackWall3D(FlowField& flowField) {}