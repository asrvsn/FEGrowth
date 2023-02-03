
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


// Project deformation gradient
void FEGrowthUncoupledMaterial::projectDeformation(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    F_true = pt.m_F;
    detF_true = pt.m_J;
    pt.m_F = pt.m_F * Fg;
    pt.m_J = pt.m_J * detFgi;
}

// Reset deformation gradient
void FEGrowthUncoupledMaterial::resetDeformation(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    pt.m_F = F_true;
    pt.m_J = detF_true;
}

//-----------------------------------------------------------------------------
//! The stress function calculates the total Cauchy stress as a sum of 
//! two terms, namely the deviatoric stress and the pressure. 
mat3ds FEGrowthUncoupledMaterial::Stress(FEMaterialPoint &mp)
{
    projectDeformation(mp);
    mat3ds S = GetBaseMaterial()->Stress(mp);
    resetDeformation(mp);
	return S;
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
    projectDeformation(mp);
    tens4ds C = GetBaseMaterial()->Tangent(mp);
    resetDeformation(mp);
    return C;
}

//-----------------------------------------------------------------------------
//! The strain energy density function calculates the total sed as a sum of
//! two terms, namely the deviatoric sed and U(J).
double FEGrowthUncoupledMaterial::StrainEnergyDensity(FEMaterialPoint &mp)
{
    projectDeformation(mp);
    double sed = GetBaseMaterial()->StrainEnergyDensity(mp);
    resetDeformation(mp);
    return sed;
}

//-----------------------------------------------------------------------------
//! The strain energy density function calculates the total sed as a sum of
//! two terms, namely the deviatoric sed and U(J).
double FEGrowthUncoupledMaterial::StrongBondSED(FEMaterialPoint &mp)
{
    projectDeformation(mp);
    double sed = GetBaseMaterial()->StrongBondSED(mp);
    resetDeformation(mp);
    return sed;
}

//-----------------------------------------------------------------------------
//! The strain energy density function calculates the total sed as a sum of
//! two terms, namely the deviatoric sed and U(J).
double FEGrowthUncoupledMaterial::WeakBondSED(FEMaterialPoint &mp)
{
    projectDeformation(mp);
    double sed = GetBaseMaterial()->WeakBondSED(mp);
    resetDeformation(mp);
    return sed;
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEGrowthUncoupledMaterial::CreateMaterialPointData()
{
    return GetBaseMaterial()->CreateMaterialPointData();
}
