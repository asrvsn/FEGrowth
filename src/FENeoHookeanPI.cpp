#include "FENeoHookeanPI.h"

//-----------------------------------------------------------------------------
// These macros define the material parameter list for the FENeoHookeanPI material.
//
// The BEGIN_FECORE_CLASS macro takes the material class and its base class as 
// parameters. 
//
// The ADD_PARAMETER macro defines the actual parameter. It takes three parameters:
// - the variable as defined in the class.
// - a range speficier which defines the valid range of the parameter
// - a string that defines the name of the parameter as it will appear in the input file.
//
// The END_PARAMETER_LIST macro just defines the end of the parameter list.
BEGIN_FECORE_CLASS(FENeoHookeanPI, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENeoHookeanPI::FENeoHookeanPI(FEModel* pfem) : FEElasticMaterial(pfem)
{
}

//-----------------------------------------------------------------------------
// This (optional) function is called during initialization. This can be done
// to do one-time initialization. Since for this material this is not required,
// this function could have been ommitted. 
// Make sure to always call the base class (usually first). 
bool FENeoHookeanPI::Init()
{
	// Don't forget the base class initialization first.
	if (FEElasticMaterial::Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
// This (optional) function is used to validate material parameters.
// It is recommended to provide a valid range during material parameter definition
// (using ADD_PARAMETER2), but any additional verification can be done here.
// For instance, here it is done to set the m_lam and m_mu parameters which depend on the 
// material parameters. 
// Make sure to always call the base class (usually first).
bool FENeoHookeanPI::Validate()
{
	if (FEElasticMaterial::Validate() == false) return false;

	// The Stress and Tangent functions are written in terms of the Lame parameters
	// so we calculate these here.
	m_lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	m_mu  = 0.5*m_E/(1+m_v);

	return true;
}

//-----------------------------------------------------------------------------
// This function needs to return the spatial (i.e. Cauchy stress) at the material point
// which is passed as a parameter. The FEMaterialPoint class contains all the state 
// variables of the material (and its base classes).
mat3ds FENeoHookeanPI::Stress(FEMaterialPoint& mp)
{
	// The FEMaterialPoint classes are stored in a linked list. The specific material
	// point data needed by this function can be accessed using the ExtractData member.
	// In this case, we want to FEElasticMaterialPoint data since it stores the deformation
	// information that is needed to evaluate the stress.
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// We'll need the deformation gradient and its determinant in this function.
	// Note that we don't take the determinant of F directly (using mat3d::det)
	// but instead use the m_J member variable of FEElasticMaterialPoint.
	mat3d &F = pt.m_F;
	double detF = pt.m_J;

	// The FEElasticMaterialPoint class defines several useful functions for 
	// evaluating strain measures, such as the left Cauchy-Green tensor.
	mat3ds b = pt.LeftCauchyGreen();

	// This creates the second-order identity tensor which we need
	// to evaluate the Cauchy stress.
	mat3dd I(1);

	// This is the actual computation of the Cauchy stress. 
	mat3ds s = (b - I)*(m_mu/detF) + I*(m_lam*log(detF)/detF);

	// The Cauchy stress is returned.
	return s;
}

//-----------------------------------------------------------------------------
// This function calculates the spatial elasticity tangent tensor. 
// It takes one parameter, the FEMaterialPoint and retursn a tens4ds object
// which is a fourth-order tensor with major and minor symmetries.
tens4ds FENeoHookeanPI::Tangent(FEMaterialPoint& mp)
{
	// As in the Stress function, we need the data from the FEElasticMaterialPoint
	// class to calculate the tangent.
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// Get the deformation gradient and its determinant
	mat3d &F = pt.m_F;
	double detF = pt.m_J;

	// Calculate the modified Lame parameters
	double lam1 = m_lam / detF;
	double mu1  = (m_mu - m_lam*log(detF)) / detF;

	// define identity tensor and some useful
	// dyadic products of the identity tensor.
	mat3dd I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds I4 = dyad4s(I);

	// evaluate the elasticity tensor
	tens4ds c = I4*(2.0*mu1) + IxI*lam1;

	// return the elasticity tensor
	return c;
}
