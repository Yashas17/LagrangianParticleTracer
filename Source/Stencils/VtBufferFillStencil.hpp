#pragma once

#include <vector>

#include "BoundaryStencil.hpp"

#include "../FlowField.hpp"
#include "../Parameters.hpp"

namespace Stencils {

  /** Stencil to set periodic boundary conditions for velocity
   */
  class VtBufferFillStencil: public BoundaryStencil<FlowField> {
  public:
    int                    cellsX;
    int                    cellsY;
    int                    cellsZ;
    std::vector<RealType> leftBuffer;
    std::vector<RealType> rightBuffer;
    std::vector<RealType> bottomBuffer;
    std::vector<RealType> topBuffer;
    std::vector<RealType> frontBuffer;
    std::vector<RealType> backBuffer;

    VtBufferFillStencil(const Parameters& parameters);
    ~VtBufferFillStencil() override = default;

    void applyLeftWall(FlowField& flowField, int i, int j) override;
    void applyRightWall(FlowField& flowField, int i, int j) override;
    void applyBottomWall(FlowField& flowField, int i, int j) override;
    void applyTopWall(FlowField& flowField, int i, int j) override;

    void applyLeftWall(FlowField& flowField, int i, int j, int k) override;
    void applyRightWall(FlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(FlowField& flowField, int i, int j, int k) override;
    void applyTopWall(FlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(FlowField& flowField, int i, int j, int k) override;
    void applyBackWall(FlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
