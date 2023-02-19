#include "FEGrowthCoupledMaterial.h"
#include <FEBioMech/FEUncoupledMaterial.h>

BEGIN_FECORE_CLASS(FEGrowthCoupledMaterial, FEGrowthMaterial)
    ADD_PROPERTY(m_mat, "base_elastic_material");
END_FECORE_CLASS();

FEGrowthCoupledMaterial::FEGrowthCoupledMaterial(FEModel* pfem) : FEGrowthMaterial(pfem)
{
    m_mat = nullptr;
}

bool FEGrowthCoupledMaterial::Init()
{
    FEUncoupledMaterial* m_umat = dynamic_cast<FEUncoupledMaterial*>((FEElasticMaterial*)m_mat);
    if (m_umat != nullptr) {
        feLogError("Base elastic material cannot be uncoupled.");
        return false;
    }
    return FEGrowthMaterial::Init();
}

bool FEGrowthCoupledMaterial::Validate()
{
    if (GetBaseMaterial()->Validate() == false) return false;
    return FEGrowthMaterial::Validate();
}

mat3ds FEGrowthCoupledMaterial::Stress(FEMaterialPoint& mp)
{
    std::function<mat3ds(FEMaterialPoint&)> f = [this](FEMaterialPoint &mp){
        return this->GetBaseMaterial()->Stress(mp);
    };
    return WithProjectedDeformation(mp, f);
}

tens4ds FEGrowthCoupledMaterial::Tangent(FEMaterialPoint &mp) 
{
    std::function<tens4ds(FEMaterialPoint&)> f = [this](FEMaterialPoint &mp) {
        return this->GetBaseMaterial()->Tangent(mp);
    };
    return WithProjectedDeformation(mp, f);
}

double FEGrowthCoupledMaterial::StrainEnergyDensity(FEMaterialPoint &mp)
{
    std::function<double(FEMaterialPoint&)> f = [this](FEMaterialPoint &mp){
        return this->GetBaseMaterial()->StrainEnergyDensity(mp);
    };
    return WithProjectedDeformation(mp, f);
}
