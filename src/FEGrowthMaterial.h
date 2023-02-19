#pragma once
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/log.h>

#include "FEGrowthMaterialPoint.h"

class FEGrowthMaterial : public FEElasticMaterial
{
    public:
        FEGrowthMaterial(FEModel* pfem);
        bool Init() override;
        bool Validate() override;

        // Applicator for deformation gradient projection
        template <typename T> T WithProjectedDeformation(FEMaterialPoint& mp, std::function<T(FEMaterialPoint&)> f)
        {
            // Must be defined inline due to template member inheritance issues
            FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
            FEGrowthMaterialPoint& gpt = *mp.ExtractData<FEGrowthMaterialPoint>();
            mat3d F = pt.m_F;
            double J = pt.m_J;
            mat3d Fgi = gpt.m_Fgi;
            double Jgi = gpt.m_Jgi;
            mat3d Fe = pt.m_F * Fgi;
            double Je = pt.m_J * Jgi;
            pt.m_F = Fe;
            pt.m_J = Je;
            T result = f(mp);
            pt.m_F = F;
            pt.m_J = J;
            return result;
        }
        
        // Projection
        mat3d GetProjection(double t, double dt, vec3d& r0, mat3d& F, mat3ds& S);
        
        // Get base elastic material
        virtual FEElasticMaterial* GetBaseMaterial() = 0;
        
        // Custom material points
        FEMaterialPointData* CreateMaterialPointData() override;

    // Parameters
    private:
        mat3d m_Fg_final;   // Growth tensor

    public:
        DECLARE_FECORE_CLASS();
};