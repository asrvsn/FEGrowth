#pragma once
#include <FECore/FEMaterialPoint.h>

// Material point storing projection information
class FEGrowthMaterialPoint : public FEMaterialPointData
{
    public:
        FEGrowthMaterialPoint(FEMaterialPointData *pt, mat3d Fg_initial, mat3d Fg_final);
        FEMaterialPointData* Copy() override;
        void Init() override;
        void Update(const FETimeInfo& timeInfo) override;
        void Serialize(DumpStream& ar) override;

    // Parameters / state
    public:
        mat3d m_Fg;
        double m_Jg;
        mat3d m_Fgi;
        double m_Jgi;
        mat3d m_Fg_final;
};