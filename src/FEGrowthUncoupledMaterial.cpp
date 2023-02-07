
#include <FECore/stdafx.h>
#include <FECore/log.h>
#include "FEGrowthUncoupledMaterial.h"
#include <iostream>

//-----------------------------------------------------------------------------
// Material parameters for FEGrowthUncoupledMaterial
BEGIN_FECORE_CLASS(FEGrowthUncoupledMaterial, FEGrowthMaterial)
    ADD_PARAMETER(m_K, FE_RANGE_GREATER_OR_EQUAL(0.0), "k")->setUnits(UNIT_PRESSURE)->MakeTopLevel(true)->setLongName("bulk modulus");
    ADD_PARAMETER(m_npmodel, "pressure_model")->setEnums("default\0NIKE3D\0Abaqus\0Abaqus (GOH)\0")->MakeTopLevel(true);
    ADD_PROPERTY(m_mat, "base_elastic_material");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEGrowthUncoupledMaterial::FEGrowthUncoupledMaterial(FEModel* pfem) : FEGrowthMaterial(pfem)
{
	m_mat = nullptr;
}

//-----------------------------------------------------------------------------
// This (optional) function is used to validate material parameters.
// It is recommended to provide a valid range during material parameter definition
// (using ADD_PARAMETER2), but any additional verification can be done here.
// For instance, here it is done to set the m_lam and m_mu parameters which depend on the 
// material parameters. 
// Make sure to always call the base class (usually first).
bool FEGrowthUncoupledMaterial::Validate()
{
	if (GetBaseMaterial()->Validate() == false) return false;
    return FEGrowthMaterial::Validate();
}

//-----------------------------------------------------------------------------
//! The stress function calculates the total Cauchy stress as a sum of 
//! two terms, namely the deviatoric stress and the pressure. 
mat3ds FEGrowthUncoupledMaterial::Stress(FEMaterialPoint &mp)
{
    std::function<mat3ds(FEMaterialPoint&)> f = [this](FEMaterialPoint &mp){
        return this->GetBaseMaterial()->Stress(mp);
    };
    return WithProjectedDeformation(mp, f);
}

//------------------------------------------------------------------------------
//! The tangent function calculates the total spatial tangent, that is it calculates
//! the push-forward of the derivative of the 2ndPK stress with respect to C. However,
//! for an uncoupled material, the 2ndPK stress decouples in a deviatoric and a 
//! dilatational component. The deviatoric tangent is provided by the particular
//! material and the dilatational component is added here.
//!
tens4ds FEGrowthUncoupledMaterial::Tangent(FEMaterialPoint &mp)
{
    std::function<tens4ds(FEMaterialPoint&)> f = [this](FEMaterialPoint &mp){
        return this->GetBaseMaterial()->Tangent(mp);
    };
    return WithProjectedDeformation(mp, f);
}

//-----------------------------------------------------------------------------
//! The strain energy density function calculates the total sed as a sum of
//! two terms, namely the deviatoric sed and U(J).
double FEGrowthUncoupledMaterial::StrainEnergyDensity(FEMaterialPoint &mp)
{
    std::function<double(FEMaterialPoint&)> f = [this](FEMaterialPoint &mp){
        return this->GetBaseMaterial()->StrainEnergyDensity(mp);
    };
    return WithProjectedDeformation(mp, f);
}

//-----------------------------------------------------------------------------
//! The strain energy density function calculates the total sed as a sum of
//! two terms, namely the deviatoric sed and U(J).
double FEGrowthUncoupledMaterial::StrongBondSED(FEMaterialPoint &mp)
{
    std::function<double(FEMaterialPoint&)> f = [this](FEMaterialPoint &mp){
        return this->GetBaseMaterial()->StrongBondSED(mp);
    };
    return WithProjectedDeformation(mp, f);
}

//-----------------------------------------------------------------------------
//! The strain energy density function calculates the total sed as a sum of
//! two terms, namely the deviatoric sed and U(J).
double FEGrowthUncoupledMaterial::WeakBondSED(FEMaterialPoint &mp)
{
    std::function<double(FEMaterialPoint&)> f = [this](FEMaterialPoint &mp){
        return this->GetBaseMaterial()->WeakBondSED(mp);
    };
    return WithProjectedDeformation(mp, f);
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEGrowthUncoupledMaterial::CreateMaterialPointData()
{
    return GetBaseMaterial()->CreateMaterialPointData();
}
