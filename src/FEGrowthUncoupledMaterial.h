// Uncoupled grown elastic plugin

#pragma once
#include <FEBioMech/FEElasticMaterial.h>
#include <FEBioMech/FEUncoupledMaterial.h>

//-----------------------------------------------------------------------------
//! Base class for uncoupled hyperelastic material formulations.

//! In FEBio, for uncoupled materials it is assumed that the strain energy function 
//! is a sum of two terms, a deviatoric strain energy term, which only depends on the 
//! deviatoric right Cauchy-Green tensor C_tilde, and a volumetric term which only 
//! depends on J, the determinant of the deformation gradient. The total Cauchy stress 
//! is therefore also a sum of two contributions, namely the deviatoric Cauchy stress
//! and a pressure term where p = dU/dJ. 

//! One of the main motivations for the alternative material interface is that some
//! finite element implementations can take advantage of the uncoupling. For example,
//! the three-field formulation integrates the two different material terms differently
//! in order to avoid locking problems for (nearly-) incompressible materials. 

//! When implementing a new material derived from this base class, the developer needs
//! to provide the deviatoric stress function and the pressure function as well as their
//! derivatives. 

class FEGrowthUncoupledMaterial : public FEElasticMaterial
{
public:
	//! constructor
	FEGrowthUncoupledMaterial(FEModel* pfem);

	//! initialization
	bool Init() override;
    bool Validate() override;

public:

//----------------->
	//! total Cauchy stress (do not overload!)
	mat3ds Stress(FEMaterialPoint& mp) final;

	//! total spatial tangent (do not overload!)
	tens4ds Tangent(FEMaterialPoint& mp) final;

	//! calculate strain energy (do not overload!)
	double StrainEnergyDensity(FEMaterialPoint& pt) final;
    double StrongBondSED(FEMaterialPoint& pt) final;
    double WeakBondSED(FEMaterialPoint& pt) final;

	// Create material point data
	FEMaterialPointData* CreateMaterialPointData() override;

    // Get base uncoupled elastic material
    FEUncoupledMaterial* GetBaseMaterial() { return m_mat; }

// Material parameters
private:
    double m_K;			            //!< bulk modulus
	int m_npmodel;                  //!< pressure model for U(J)
    FEUncoupledMaterial* m_mat;     // Base elastic material
    double m_g_iso;                 // Isotropic growth factor
    
private:
    void projectDeformation(FEMaterialPoint& pt);
    void resetDeformation(FEMaterialPoint& pt);
    mat3ds  Fg;             // Growth tensor
    mat3ds  Fgi;            // Inverse growth tensor
    double  detFg;          // 
    double  detFgi;         // 
    mat3d   F_true;         // True deformation gradient
    double  detF_true;      // True volumetric change
    
public:

	DECLARE_FECORE_CLASS();
};
