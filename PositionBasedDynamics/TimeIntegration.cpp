#include "TimeIntegration.h"

using namespace PBD;


// ----------------------------------------------------------------------------------------------
void TimeIntegration::semiImplicitEuler(
	const Real h, 
	const Real mass, 
	Vector3r &position,
	Vector3r &velocity,
	const Vector3r &acceleration)
{				
	if (mass != 0.0)
	{
		velocity += acceleration * h;
		position += velocity * h;
	}
}

#include <iostream>

// ----------------------------------------------------------------------------------------------
void TimeIntegration::semiImplicitEulerRotation(
	const Real h,
	const Real mass,
	const Matrix3r& invertiaW,
	const Matrix3r &invInertiaW,
	Quaternionr &rotation,
	Vector3r &angularVelocity,	
	const Vector3r &torque)
{
	if (mass != 0.0)
	{
		auto angCross = angularVelocity.cross( invertiaW * angularVelocity );
		angularVelocity += h * invInertiaW * (torque - ( angCross ));

		Quaternionr angVelQ(0.0, angularVelocity[0], angularVelocity[1], angularVelocity[2]);
		if ( ( angVelQ.x() != 0 || angVelQ.y() != 0 || angVelQ.z() != 0 || angVelQ.w() != 0 )
			&& ( rotation.x() != 0 || rotation.y() != 0 || rotation.z() != 0 || rotation.w() != 0 ) )
				
		{
			static bool bprint = false;
			if ( bprint )
			{
				std::cout << "angVelQ: " << angVelQ.coeffs().transpose() << std::endl;
				std::cout << "rotation: " << rotation.coeffs().transpose() << std::endl;

				auto tmp = angVelQ * rotation;
				auto tmp2 = rotation * angVelQ;
				std::cout << "tmp: " << tmp.coeffs().transpose() << std::endl;
				std::cout << "tmp2: " << tmp2.coeffs().transpose() << std::endl;
			}
		}
		rotation.coeffs() += h * 0.5 * (angVelQ * rotation).coeffs();
		rotation.normalize();
	}
}

// ----------------------------------------------------------------------------------------------
void TimeIntegration::velocityUpdateFirstOrder(
	const Real h,
	const Real mass,
	const Vector3r &position,
	const Vector3r &oldPosition,
	Vector3r &velocity)
{
	if (mass != 0.0)
		velocity = (1.0 / h) * (position - oldPosition);
}

// ----------------------------------------------------------------------------------------------
void TimeIntegration::angularVelocityUpdateFirstOrder(
	const Real h,
	const Real mass,
	const Quaternionr &rotation,
	const Quaternionr &oldRotation,
	Vector3r &angularVelocity)
{
	if (mass != 0.0)
	{
		const Quaternionr relRot = (rotation * oldRotation.conjugate());
		angularVelocity = relRot.vec() *(2.0 / h);
	}
}

// ----------------------------------------------------------------------------------------------
void TimeIntegration::velocityUpdateSecondOrder(
	const Real h,
	const Real mass,
	const Vector3r &position,
	const Vector3r &oldPosition,
	const Vector3r &positionOfLastStep,
	Vector3r &velocity)
{
	if (mass != 0.0)
		velocity = (1.0 / h) * (1.5*position - 2.0*oldPosition + 0.5*positionOfLastStep);
}

// ----------------------------------------------------------------------------------------------
void TimeIntegration::angularVelocityUpdateSecondOrder(
	const Real h,
	const Real mass,
	const Quaternionr &rotation,				
	const Quaternionr &oldRotation,			
	const Quaternionr &rotationOfLastStep,	
	Vector3r &angularVelocity)
{
	// ToDo: is still first order
	if (mass != 0.0)
	{
		const Quaternionr relRot = (rotation * oldRotation.conjugate());
		angularVelocity = relRot.vec() *(2.0 / h);
	}
}