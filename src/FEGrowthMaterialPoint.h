#pragma once
#include <FEBioMech/FEElasticMaterialPoint.h>

// Material point storing projection information
class FEGrowthMaterialPoint : public FEElasticMaterialPoint
{
    public:
        FEGrowthMaterialPoint(FEMaterialPointData* mp = nullptr);
        FEMaterialPointData* Copy() override;
        void Init() override;
        void Update(const FETimeInfo& timeInfo) override;
        void Serialize(DumpStream& ar) override;
        void SetTarget(mat3d Fg_final, double t0, double tf);

    // Parameters / state
    public:
        mat3d m_Fg_initial; // Initial growth tensor
        mat3d m_Fg;         // Current growth tensor
        double m_Jg;        // Current determinant
        mat3d m_Fgi;        // Current inverse
        double m_Jgi;       // Current inverse determinant
        mat3d m_Fg_final;   // Target growth tensor
        double m_t0;        // Start time
        double m_tf;        // End time
};