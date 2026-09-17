/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of turbinesFoam, which is based on OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "turbineALSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(turbineALSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        turbineALSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::turbineALSource::rotateVector
(
    vector& vectorToRotate,
    vector rotationPoint,
    vector axis,
    scalar radians
)
{
    // Declare and define the rotation matrix (from SOWFA)
    tensor RM;
    scalar angle = radians;
    RM.xx() = Foam::sqr(axis.x())
            + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle);
    RM.xy() = axis.x() * axis.y()
            * (1.0 - Foam::cos(angle)) - axis.z() * Foam::sin(angle);
    RM.xz() = axis.x() * axis.z()
            * (1.0 - Foam::cos(angle)) + axis.y() * Foam::sin(angle);
    RM.yx() = axis.x() * axis.y()
            * (1.0 - Foam::cos(angle)) + axis.z() * Foam::sin(angle);
    RM.yy() = Foam::sqr(axis.y())
            + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y() * axis.z()
            * (1.0 - Foam::cos(angle)) - axis.x() * Foam::sin(angle);
    RM.zx() = axis.x() * axis.z()
            * (1.0 - Foam::cos(angle)) - axis.y() * Foam::sin(angle);
    RM.zy() = axis.y() * axis.z()
            * (1.0 - Foam::cos(angle)) + axis.x() * Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z())
            + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);

    // Rotation matrices make a rotation about the origin, so need to subtract
    // rotation point off the point to be rotated.
    vectorToRotate -= rotationPoint;

    // Perform the rotation.
    vectorToRotate = RM & vectorToRotate;

    // Return the rotated point to its new location relative to the rotation
    // point
    vectorToRotate += rotationPoint;
}


void Foam::fv::turbineALSource::createCoordinateSystem()
{
    // Should be unique for each type of turbine
}


void Foam::fv::turbineALSource::createBlades()
{
    // Should be unique for each type of turbine
}


void Foam::fv::turbineALSource::createOutputFile()
{
    fileName dir;

    if (Pstream::parRun())
    {
        dir = time_.path()/"../postProcessing/turbines"
            / time_.timeName();
    }
    else
    {
        dir = time_.path()/"postProcessing/turbines"
            / time_.timeName();
    }

    if (not isDir(dir))
    {
        mkDir(dir);
    }

    outputFile_ = new OFstream(dir/name_ + ".csv");

    *outputFile_<< "time,angle_deg,tsr,cp,cd,ct";

    forAll(blades_, i)
    {
        *outputFile_<< ",cd_" << bladeNames_[i];
        *outputFile_<< ",ct_" << bladeNames_[i];
    }

    *outputFile_<< endl;
}


void Foam::fv::turbineALSource::updateTSROmega()
{
    // Update tip speed ratio and omega
    scalar theta = degToRad(angleDeg_);
    tipSpeedRatio_ = meanTSR_ + tsrAmplitude_*cos(nBlades_*(theta - tsrPhase_));
    omega_ = tipSpeedRatio_*mag(freeStreamVelocity_)/rotorRadius_;
}


void Foam::fv::turbineALSource::updateRestartState()
{
    state_.set("angleDeg", angleDeg_);
    state_.set("checkpointTime", time_.value());
}


Foam::scalar Foam::fv::turbineALSource::wrappedAngleDeg() const
{
    scalar angle = fmod(angleDeg_, 360.0);

    if (angle < 0.0)
    {
        angle += 360.0;
    }

    return angle;
}


void Foam::fv::turbineALSource::rotate()
{
    scalar deltaT = time_.deltaT().value();
    scalar radians = omega_*deltaT;
    rotate(radians);
    angleDeg_ += radToDeg(radians);
    lastRotationTime_ = time_.value();
    updateRestartState();
    updateTSROmega();
}


void Foam::fv::turbineALSource::rotate(scalar radians)
{
    // Should be defined for each turbine type
}


void Foam::fv::turbineALSource::printPerf()
{
    Info<< "Azimuthal angle (degrees) of " << name_ << ": " << angleDeg_
        << endl;
    Info<< "Tip speed ratio of " << name_ << ": " << tipSpeedRatio_ << endl;
    Info<< "Power coefficient from " << name_ << ": " << powerCoefficient_
        << endl;
    Info<< "Rotor drag coefficient from " << name_ << ": " << dragCoefficient_
        << endl << endl;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::turbineALSource::turbineALSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    time_(mesh.time()),
    state_
    (
        IOobject
        (
            name + "State",
            time_.timeName(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    azimuthStateInitialized_(false),
    lastRotationTime_(time_.value()),
    rhoRef_(1.0),
    omega_(0.0),
    angleDeg_(state_.lookupOrDefault<scalar>("angleDeg", 0.0)),
    nBlades_(0),
    freeStreamVelocity_(vector::zero),
    forceField_
    (
        IOobject
        (
            "force." + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "force",
            dimForce/dimVolume/dimDensity,
            vector::zero
        )
    ),
    frontalArea_(0.0),
    powerCoefficient_(0.0),
    dragCoefficient_(0.0),
    torqueCoefficient_(0.0)
{
    forceField_.write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::turbineALSource::~turbineALSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::turbineALSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Should be unique for each turbine type
}


void Foam::fv::turbineALSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Should be unique for each turbine type
}


void Foam::fv::turbineALSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    // Should be unique for each turbine type
}


void Foam::fv::turbineALSource::printCoeffs() const
{
    Info<< "Number of blades: " << nBlades_ << endl;
}


void Foam::fv::turbineALSource::writePerf()
{
    *outputFile_<< time_.value() << "," << angleDeg_ << ","
                << tipSpeedRatio_ << "," << powerCoefficient_ << ","
                << dragCoefficient_ << "," << torqueCoefficient_;

    // Write power, drag, and torque coefficients for each blade
    forAll(blades_, i)
    {
        // Write drag (thrust) coefficient contribution from blade
        scalar bladeCd = blades_[i].force() & freeStreamDirection_
            / (0.5*frontalArea_*magSqr(freeStreamVelocity_));
        *outputFile_<< "," << bladeCd;
        // Write torque coefficient contribution from blade
        scalar bladeTorque = bladeMoments_[i] & axis_;
        scalar bladeCt = bladeTorque
            / (0.5*frontalArea_*rotorRadius_* magSqr(freeStreamVelocity_));
        *outputFile_<< "," << bladeCt;
    }

    *outputFile_<< endl;
}


void Foam::fv::turbineALSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::turbineALSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        // Read coordinate system/geometry invariant properties
        coeffs_.lookup("origin") >> origin_;
        coeffs_.lookup("axis") >> axis_;
        axis_ /= mag(axis_);
        coeffs_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
        coeffs_.lookup("tipSpeedRatio") >> meanTSR_;
        coeffs_.lookup("rotorRadius") >> rotorRadius_;
        tsrAmplitude_ = coeffs_.lookupOrDefault("tsrAmplitude", 0.0);
        tsrPhase_ = coeffs_.lookupOrDefault("tsrPhase", 0.0);

        // Initialize azimuth only once. Runtime dictionary re-reads must not
        // overwrite the angle accumulated during the simulation.
        if (!azimuthStateInitialized_)
        {
            if (state_.found("angleDeg"))
            {
                Info<< "Restored azimuthal rotation of " << name_ << ": "
                    << angleDeg_ << " degrees" << endl;
            }
            else
            {
                angleDeg_ = coeffs_.lookupOrDefault<scalar>
                (
                    "restartAngleDeg",
                    0.0
                );

                if (mag(time_.value()) > SMALL)
                {
                    WarningInFunction
                        << "No saved azimuth state was found for " << name_
                        << " at time " << time_.value() << ". Using "
                        << angleDeg_ << " degrees. For legacy restarts, set "
                        << "restartAngleDeg in the turbine coefficients."
                        << endl;
                }
            }

            azimuthStateInitialized_ = true;
            updateRestartState();
        }

        // Get blade information
        bladesDict_ = coeffs_.subDict("blades");
        nBlades_ = bladesDict_.keys().size();
        bladeNames_ = bladesDict_.toc();
        bladeMoments_.setSize(nBlades_);

        // Set tip speed ratio and omega
        updateTSROmega();

        // Get dynamic stall subdict
        dynamicStallDict_ = coeffs_.subOrEmptyDict("dynamicStall");

	// Get filteredLiftingLine subdict
	filteredLiftingLineDict_ = coeffs_.subOrEmptyDict("filteredLiftingLine"); 
        // Get profiles information
        profileData_ = coeffs_.subDict("profileData");

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
