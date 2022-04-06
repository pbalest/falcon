//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ConvectiveLineSourceRayKernel.h"

// MOOSE includes
#include "Function.h"

registerMooseObject("FalconApp", ConvectiveLineSourceRayKernel);
registerMooseObject("FalconApp", ADConvectiveLineSourceRayKernel);

template <bool is_ad>
InputParameters
ConvectiveLineSourceRayKernelTempl<is_ad>::validParams()
{
  InputParameters params = GenericRayKernel<is_ad>::validParams();

  params.addClassDescription(
      "Demonstrates the multiple ways that scalar values can be introduced "
      "into RayKernels, e.g. (controllable) constants, functions, "
      "postprocessors, and data on rays. Implements the weak form $(\\psi_i, -f)$ along a line.");

  params.addParam<Real>("heated_perimeter", 1.0, "Coefficient to multiply by the line source term");
  params.addParam<FunctionName>("T_infinity", "1", "A function that describes the line source");
  params.addParam<FunctionName>("heat_transfer_coefficient", "1", "A function that describes the line source");
  params.addParam<PostprocessorName>(
      "postprocessor", 1, "A postprocessor whose value is multiplied by the line source");
  params.addParam<std::vector<std::string>>(
      "ray_data_factor_names", "The names of the Ray data to scale the source by (if any)");
  params.addParam<std::vector<std::string>>(
      "ray_aux_data_factor_names", "The names of the Ray aux data to scale the source by (if any)");

  params.declareControllable("value");

  return params;
}

template <bool is_ad>
ConvectiveLineSourceRayKernelTempl<is_ad>::ConvectiveLineSourceRayKernelTempl(const InputParameters & params)
  : GenericRayKernel<is_ad>(params),
    _heated_perimeter(this->template getParam<Real>("heated_perimeter")),
    _T_infinity(getFunction("T_infinity")),
    _htc(getFunction("heat_transfer_coefficient")),
    _postprocessor(getPostprocessorValue("postprocessor")),
    _ray_data_factor_indices(this->_study.getRayDataIndices(
        this->template getParam<std::vector<std::string>>("ray_data_factor_names"))),
    _ray_aux_data_factor_indices(this->_study.getRayAuxDataIndices(
        this->template getParam<std::vector<std::string>>("ray_aux_data_factor_names")))
{
}

template <bool is_ad>
GenericReal<is_ad>
ConvectiveLineSourceRayKernelTempl<is_ad>::computeQpResidual()
{
  Real factor = _heated_perimeter * _postprocessor * _htc.value(_t, _q_point[_qp]);

  // Scale by any Ray data and aux data if given
  for (const auto index : _ray_data_factor_indices)
    factor *= currentRay()->data(index);
  for (const auto index : _ray_aux_data_factor_indices)
    factor *= currentRay()->auxData(index);

  return _test[_i][_qp] * -factor*(_T_infinity.value(_t, _q_point[_qp]) -_u[_qp]);
}

template class ConvectiveLineSourceRayKernelTempl<false>;
template class ConvectiveLineSourceRayKernelTempl<true>;
