#include "FEGrowthMaterial.h"

BEGIN_FECORE_CLASS(FEGrowthMaterial, FEElasticMaterial)
    ADD_PARAMETER(m_Fg_initial, "initial_growth_tensor");
    ADD_PARAMETER(m_Fg_final, "final_growth_tensor");
END_FECORE_CLASS();

FEGrowthMaterial::FEGrowthMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_Fg_final = mat3d(1,0,0,0,1,0,0,0,1);
    m_Fg_initial = mat3d(1,0,0,0,1,0,0,0,1);
}

bool FEGrowthMaterial::Init()
{
    return FEElasticMaterial::Init();
}

bool FEGrowthMaterial::Validate()
{
    if (m_Fg_final.det() <= 0.0 || m_Fg_initial.det() <= 0.0) {
        feLogError("Growth tensors must have positive determinant.");
        return false;
    }
	return FEElasticMaterial::Validate();
}

FEMaterialPointData* FEGrowthMaterial::CreateMaterialPointData()
{
    return new FEGrowthMaterialPoint(
        GetBaseMaterial()->CreateMaterialPointData(),
        m_Fg_initial,
        m_Fg_final
    );
}