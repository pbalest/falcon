#pragma once

#include "GenericRayKernel.h"

// Forward declarations
class Function;

template <bool is_ad>
class ConvectiveLineSourceRayKernelTempl : public GenericRayKernel<is_ad>
{
public:
  ConvectiveLineSourceRayKernelTempl(const InputParameters & params);

  static InputParameters validParams();

protected:
  virtual GenericReal<is_ad> computeQpResidual() override;

  /// Heated perimeter to scale the linear heat source.
  const Real & _heated_perimeter;

  /// Convective values.
  const Function & _T_infinity;
  const Function & _htc;

  /// Indices into the Ray data that we want to scale the residual by (may be empty)
  const std::vector<RayDataIndex> _ray_data_factor_indices;
  /// Indices into the Ray aux data that we want to scale the residual by (may be empty)
  const std::vector<RayDataIndex> _ray_aux_data_factor_indices;

  usingGenericRayKernelMembers;
};

typedef ConvectiveLineSourceRayKernelTempl<false> ConvectiveLineSourceRayKernel;
typedef ConvectiveLineSourceRayKernelTempl<true> ADConvectiveLineSourceRayKernel;
