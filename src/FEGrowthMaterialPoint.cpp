#include "FEGrowthMaterialPoint.h"

FEGrowthMaterialPoint::FEGrowthMaterialPoint(FEMaterialPointData* pt, mat3d Fg_initial, mat3d Fg_final) : FEMaterialPointData(pt)
{
    m_Fg = Fg_initial;
    m_Jg = m_Fg.det();
    m_Fgi = m_Fg.inverse();
    m_Jgi = 1 / m_Jg;
    m_Fg_final = Fg_final;
}

FEMaterialPointData* FEGrowthMaterialPoint::Copy()
{
    FEGrowthMaterialPoint* pt = new FEGrowthMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

void FEGrowthMaterialPoint::Init()
{
    FEMaterialPointData::Init();
}

// Main function - update new virtual frame
// Can depend on time, position, strain, stress
void FEGrowthMaterialPoint::Update(const FETimeInfo& timeInfo)
{
    // Simple forward Euler
    mat3d dFdt = m_Fg_final - m_Fg;
    m_Fg = m_Fg + dFdt * timeInfo.timeIncrement;

    // Update derived values
    m_Jg = m_Fg.det();
    m_Fgi = m_Fg.inverse();
    m_Jgi = 1 / m_Jg;
}

void FEGrowthMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);
}