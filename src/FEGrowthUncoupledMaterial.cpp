
#include <FECore/stdafx.h>
#include <FECore/log.h>
#include "FEGrowthUncoupledMaterial.h"
#include <iostream>

//-----------------------------------------------------------------------------
// Material parameters for FEGrowthUncoupledMaterial
BEGIN_FECORE_CLASS(FEGrowthUncoupledMaterial, FEElasticMaterial)
    ADD_PARAMETER(m_K, FE_RANGE_GREATER_OR_EQUAL(0.0), "k")->setUnits(UNIT_PRESSURE)->MakeTopLevel(true)->setLongName("bulk modulus");
    ADD_PARAMETER(m_npmodel, "pressure_model")->setEnums("default\0NIKE3D\0Abaqus\0Abaqus (GOH)\0")->MakeTopLevel(true);
    ADD_PROPERTY(m_mat, "base_elastic_material");
    ADD_PARAMETER(m_g_iso, FE_RANGE_GREATER(0.0), "isotropic_growth_factor");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEGrowthUncoupledMaterial::FEGrowthUncoupledMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_mat = nullptr;
    m_g_iso = 1.0;
}

//-------------------------------------------------≈----------------------------
bool FEGrowthUncoupledMaterial::Init()
{
	return FEElasticMaterial::Init();
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

	// Define growth tensors
	Fg = m_g_iso * mat3ds(1,1,1,0,0,0);
	Fgi = Fg.inverse();
	detFg = Fg.det();
	detFgi = 1 / detFg;

    if (detFg <= 0.0) {
        feLogError("detG must be positive.");
        return false;
    }

	return FEElasticMaterial::Validate();
}


// Execute callable with projected deformation gradient
template <typename T>
T FEGrowthUncoupledMaterial::WithProjectedDeformation(FEMaterialPoint& mp, std::function<T(FEMaterialPoint&)> f)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3d F = pt.m_F;
    double J = pt.m_J;
    mat3d Fe = pt.m_F * Fgi;
    double Je = pt.m_J * detFgi;
    pt.m_F = Fe;
    pt.m_J = Je;
    T result = f(mp);
    pt.m_F = F;
    pt.m_J = J;
    return result;
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
