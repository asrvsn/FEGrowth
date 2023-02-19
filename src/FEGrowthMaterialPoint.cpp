#include <FECore/log.h>
#include <iostream>
#include "FEGrowthMaterialPoint.h"

FEGrowthMaterialPoint::FEGrowthMaterialPoint(FEMaterialPointData* pt) : FEElasticMaterialPoint(pt)
{
    m_Fg_initial.unit();
    m_Fg.unit();
    m_Jg = 1;
    m_Fgi.unit();
    m_Jgi = 1;
    m_Fg_final.unit();
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
    const double t = timeInfo.currentTime;
    const double t_perc = (t - m_t0) / (m_tf - m_t0);
    
    // // Simple forward Euler
    // mat3d dFdt = m_Fg_final - m_Fg;
    // m_Fg = m_Fg + dFdt * timeInfo.timeIncrement;

    // Linear interpolate from initial to final
    m_Fg = m_Fg_initial * (1 - t_perc) + m_Fg_final * t_perc;

    // Update derived values
    m_Jg = m_Fg.det();
    m_Fgi = m_Fg.inverse();
    m_Jgi = 1 / m_Jg;
}

void FEGrowthMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);
}

void FEGrowthMaterialPoint::SetTarget(mat3d Fg_final, double t0, double tf) {
    m_Fg_final = Fg_final;
    m_t0 = t0;
    m_tf = tf;
}