#pragma once
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/log.h>

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

    // Parameters
    private:
        mat3d m_Fg;         // Growth tensor
        double m_Fg_ramp;   // Ramp factor for growth

    // Private instance variables
    private:
        mat3d Fgi;      // Inverse growth tensor
        double detFg;   // 
        double detFgi;  //

    public:
        DECLARE_FECORE_CLASS();
};