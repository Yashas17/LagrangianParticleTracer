#include "StdAfx.hpp"

#include "VelocityBufferFillStencil.hpp"

#include "Parameters.hpp"

#include "../FlowField.hpp"

Stencils::VelocityBufferFillStencil::VelocityBufferFillStencil(const Parameters& parameters):
  GhostLayerStencil(parameters),
  parameters_(parameters) {}

Stencils::VelocityBufferFillStencil::~VelocityBufferFillStencil() {}

void Stencils::VelocityBufferFillStencil::applyLeftWall2D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyRightWall2D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyBottomWall2D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyTopWall2D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyLeftWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyRightWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyBottomWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyTopWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyFrontWall3D(FlowField& flowField) {}

void Stencils::VelocityBufferFillStencil::applyBackWall3D(FlowField& flowField) {}