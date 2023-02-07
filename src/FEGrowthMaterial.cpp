#include "FEGrowthMaterial.h"

BEGIN_FECORE_CLASS(FEGrowthMaterial, FEElasticMaterial)
    ADD_PARAMETER(m_Fg, "growth_tensor");
    ADD_PARAMETER(m_Fg_ramp, FE_RANGE_GREATER_OR_EQUAL(0.0), "growth_ramp");
END_FECORE_CLASS();

FEGrowthMaterial::FEGrowthMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_Fg = mat3d(1,0,0,0,1,0,0,0,1);
    m_Fg_ramp = 1.0;
}

bool FEGrowthMaterial::Init()
{
    return FEElasticMaterial::Init();
}

bool FEGrowthMaterial::Validate()
{
    // Define growth 
	detFg = m_Fg.det();
    if (detFg <= 0.0) {
        feLogError("Growth tensor must have positive determinant.");
        return false;
    }
	Fgi = m_Fg.inverse();
	detFgi = 1 / detFg;

	return FEElasticMaterial::Validate();
}