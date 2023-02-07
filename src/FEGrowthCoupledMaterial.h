#pragma once
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/log.h>
#include "FEGrowthMaterial.h"

class FEGrowthCoupledMaterial : public FEGrowthMaterial
{
    public:
        FEGrowthCoupledMaterial(FEModel* pfem);
        bool Init() override;
        bool Validate() override;

    public:
        mat3ds Stress(FEMaterialPoint& mp) override;
        tens4ds Tangent(FEMaterialPoint& mp) override;
        double StrainEnergyDensity(FEMaterialPoint& mp) override;
        FEMaterialPointData* CreateMaterialPointData() override;
        FEElasticMaterial* GetBaseMaterial() { return m_mat; }

    private:
        FEElasticMaterial* m_mat;   // Base elastic material
        
    public:
        DECLARE_FECORE_CLASS();
};