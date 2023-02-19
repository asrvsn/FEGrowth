#include <iostream>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>

#include "FEGrowthMaterial.h"

BEGIN_FECORE_CLASS(FEGrowthMaterial, FEElasticMaterial)
    ADD_PARAMETER(m_Fg_final, "final_growth_tensor");
END_FECORE_CLASS();

FEGrowthMaterial::FEGrowthMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_Fg_final.unit();
}

bool FEGrowthMaterial::Init()
{
    return FEElasticMaterial::Init();
}

bool FEGrowthMaterial::Validate()
{
    if (m_Fg_final.det() <= 0.0) {
        feLogError("Growth tensors must have positive determinant.");
        return false;
    }
	return FEElasticMaterial::Validate();
}

FEMaterialPointData* FEGrowthMaterial::CreateMaterialPointData()
{
    FEAnalysis* currentStep = GetFEModel()->GetCurrentStep();
    double t0 = currentStep->m_tstart;
    double tf = (currentStep->m_dt0 * currentStep->m_ntime) + t0;
    auto pt = new FEGrowthMaterialPoint(GetBaseMaterial()->CreateMaterialPointData());
    pt->SetTarget(m_Fg_final, t0, tf);
    return pt;
}