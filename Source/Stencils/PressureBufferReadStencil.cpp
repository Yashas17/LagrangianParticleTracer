#pragma once
#include "StdAfx.hpp"

#include "PressureBufferReadStencil.hpp"

#include "Parameters.hpp"

#include "../FlowField.hpp"

Stencils::PressureBufferReadStencil::PressureBufferReadStencil(const Parameters& parameters):
  GhostLayerStencil(parameters),
  parameters_(parameters) {}

Stencils::PressureBufferReadStencil::~PressureBufferReadStencil() {}

void Stencils::PressureBufferReadStencil::applyLeftWall2D(FlowField& flowField) {}

void Stencils::PressureBufferReadStencil::applyRightWall2D(FlowField& flowField) {}

void Stencils::PressureBufferReadStencil::applyBottomWall2D(FlowField& flowField) {}

void Stencils::PressureBufferReadStencil::applyTopWall2D(FlowField& flowField) {}

void Stencils::PressureBufferReadStencil::applyLeftWall3D(FlowField& flowField) {}

void Stencils::PressureBufferReadStencil::applyRightWall3D(FlowField& flowField) {}

void Stencils::PressureBufferReadStencil::applyBottomWall3D(FlowField& flowField) {}

void Stencils::PressureBufferReadStencil::applyTopWall3D(FlowField& flowField) {}

void Stencils::PressureBufferReadStencil::applyFrontWall3D(FlowField& flowField) {}

void Stencils::PressureBufferReadStencil::applyBackWall3D(FlowField& flowField) {}