#pragma once
#include "StdAfx.hpp"

#include "PressureBufferFillStencil.hpp"

#include "Parameters.hpp"

#include "../FlowField.hpp"

Stencils::PressureBufferFillStencil::PressureBufferFillStencil(const Parameters& parameters):
  GhostLayerStencil(parameters),
  parameters_(parameters) {}

Stencils::PressureBufferFillStencil::~PressureBufferFillStencil() {}

void Stencils::PressureBufferFillStencil::applyLeftWall2D(FlowField& flowField) {}

void Stencils::PressureBufferFillStencil::applyRightWall2D(FlowField& flowField) {}

void Stencils::PressureBufferFillStencil::applyBottomWall2D(FlowField& flowField) {}

void Stencils::PressureBufferFillStencil::applyTopWall2D(FlowField& flowField) {}

void Stencils::PressureBufferFillStencil::applyLeftWall3D(FlowField& flowField) {}

void Stencils::PressureBufferFillStencil::applyRightWall3D(FlowField& flowField) {}

void Stencils::PressureBufferFillStencil::applyBottomWall3D(FlowField& flowField) {}

void Stencils::PressureBufferFillStencil::applyTopWall3D(FlowField& flowField) {}

void Stencils::PressureBufferFillStencil::applyFrontWall3D(FlowField& flowField) {}

void Stencils::PressureBufferFillStencil::applyBackWall3D(FlowField& flowField) {}