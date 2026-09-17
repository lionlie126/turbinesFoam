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

#include "actuatorLineElement.H"
#include "addToRunTimeSelectionTable.H"
#include "geometricOneField.H"
#include "fvMatrices.H"
#include "syncTools.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuatorLineElement, 0);
    defineRunTimeSelectionTable(actuatorLineElement, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fv::actuatorLineElement::read()
{
    // Parse dictionary
    dict_.lookup("position") >> position_;
    dict_.lookup("chordLength") >> chordLength_;
    dict_.lookup("chordDirection") >> chordDirection_;
    dict_.lookup("chordRefDirection") >> chordRefDirection_;
    dict_.lookup("chordMount") >> chordMount_;
    dict_.lookup("spanLength") >> spanLength_;
    dict_.lookup("spanDirection") >> spanDirection_;
    dict_.lookup("thickness") >> thickness_;
    dict_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
    freeStreamDirection_ = freeStreamVelocity_/mag(freeStreamVelocity_);
    dict_.lookup("bladeElementPitch") >> bladeElementPitch_; 
    dict_.lookup("bladeChordMax") >> bladeChordMax;
    dict_.lookup("rootDistance") >> rootDistance_;
    dict_.lookup("aspectRatio") >> aspectRatio_;
    dict_.lookup("averageChordLength") >> avegChordLength_; 
    dict_.lookup("bladeRadius") >> bladeRadius_;
    dict_.lookup("bladePointRadius") >> bladePointRadius;
    dict_.lookup("numOfElement") >> numOfElement;
   // dict_.lookup("totalNumberOfElement") >> totalNumberOfElement;
   
    //word defaultVelDragCorrType = "none"; 
    //word velDragCorrType_ = dict_.lookupOrDefault("velDragCorrType", defaultVelDragCorrType);
    
    dict_.lookup("velocitySampleRadius") >> velocitySampleRadius_;
    dict_.lookup("nVelocitySamples") >> nVelocitySamples_;


    if (dict_.found("filteredLiftingLine"))
    {
	dictionary fllDict = dict_.subDict("filteredLiftingLine");
	filteredLiftingLineActive_ = fllDict.lookupOrDefault("active", false);
	filteredLiftingLineWriteOutput_ = fllDict.lookupOrDefault("writeOutputFll", false);
	relaxFactorLC = fllDict.lookupOrDefault("relaxationFactorLC", 0.1);
        relaxFactorIV = fllDict.lookupOrDefault("relaxationFactorIV", 0.1);

    }


    // Create dynamic stall model if found
    if (dict_.found("dynamicStall"))
    {
        dictionary dsDict = dict_.subDict("dynamicStall");
        word dsName;
        dsDict.lookup("dynamicStallModel") >> dsName;
        dynamicStall_ = dynamicStallModel::New
        (
            dsDict,
            dsName,
            mesh_.time(),
            profileData_
        );
        dsDict.lookup("active") >> dynamicStallActive_;
        
    }

    // Read flow curvature correction subdictionary
    if (dict_.found("flowCurvature"))
    {
        dictionary fcDict = dict_.subDict("flowCurvature");
        flowCurvatureActive_ = fcDict.lookupOrDefault("active", false);
        word defaultName = "none";
        flowCurvatureModelName_ = fcDict.lookupOrDefault
        (
            "flowCurvatureModel",
            defaultName
        );
    }

    // Read nu from object registry
    const dictionary& transportProperties = mesh_.lookupObject<IOdictionary>
    (
        "transportProperties"
    );
    dimensionedScalar nu;
    transportProperties.lookup("nu") >> nu;
    nu_ = nu.value();

    // Read velocity evaluation type once (avoids repeated dict lookup per timestep)
    dict_.lookup("velEvalType") >> velEvalType_;

    // Read EVM settings once from disk if EVM is active
    if (velEvalType_ == "EVM")
    {
        IOdictionary EVMsettings
        (
            IOobject
            (
                "EVMsettings",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        EVMsettings.lookup("distance")     >> evmDistance_;
        EVMsettings.lookup("length")       >> evmSamplingLength_;
        EVMsettings.lookup("numEVMPoints") >> evmNumPoints_;
    }

    // Read writePerf switch
    dict_.lookup("writePerf") >> writePerf_;

    if (debug)
    {
        Info<< "actuatorLineElement properties:" << endl;
        Info<< "Position: " << position_ << endl;
        Info<< "chordLength: " << chordLength_ << endl;
        Info<< "chordDirection: " << chordDirection_ << endl;
        Info<< "spanLength: " << spanLength_ << endl;
        Info<< "spanDirection: " << spanDirection_ << endl;
        Info<< "writePerf: " << writePerf_ << endl;
    }
}


void Foam::fv::actuatorLineElement::rotateVector
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




Foam::label Foam::fv::actuatorLineElement::findCell
(
    const point& location
)
{
    if (Pstream::parRun())
    {
        // Use contains() rather than containsInside() so that points lying
        // exactly on a processor subdomain face boundary are not false-rejected.
        // containsInside() is a strict interior test — points on the bounding
        // box surface return false on all processors simultaneously, which is
        // the root cause of systematic misses for tip elements as they rotate.
        if (meshBoundBox_.contains(location))
        {
            if (debug > 1)
            {
                Pout<< "Looking for cell containing " << location
                    << " inside bounding box:" << endl
                    << meshBoundBox_ << endl;
            }
            return mesh_.findCell(location);
        }
        else
        {
            // Point is outside this processor's subdomain — expected in parallel.
            // Only log at debug level 2+ to avoid flooding with O(P*N_EVM) prints.
            if (debug > 1)
            {
                Pout<< "Cell not inside this processor's domain: "
                    << meshBoundBox_ << endl;
            }
            return -1;
        }
    }
    else
    {
        return mesh_.findCell(location);;
    }
}


void Foam::fv::actuatorLineElement::lookupCoefficients(scalar angleOfAttack)
{
    liftCoefficient_ = profileData_.liftCoefficient(angleOfAttack);
    dragCoefficient_ = profileData_.dragCoefficient(angleOfAttack);
    momentCoefficient_ = profileData_.momentCoefficient(angleOfAttack);
}


Foam::vector Foam::fv::actuatorLineElement::calcProjectionEpsilon()
{
    word gaussianRadiusType_;
    dict_.lookup("gaussianRadiusType") >> gaussianRadiusType_; 
	
    epsilon = vector(0.0, 0.0, 0.0);
    
	IOdictionary epsilonSettings
	(
		IOobject
		(
			"epsilonSettings",
			mesh_.time().constant(),
			mesh_,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);
    
    dictionary coeffsEps = epsilonSettings.subDict(gaussianRadiusType_ + "Coeffs");
    scalar bladeEpsilonFactor;
    scalar bladeEpsilonFactorChord;
    scalar bladeEpsilonFactorThickness;
    //Info<< "position:" << position_ << endl; 
    const scalarField& V = mesh_.V();
    label posCellI = findCell(position_);
       
    if (posCellI >= 0)
    {
	scalar cellLen = Foam::cbrt(V[posCellI]);
	
	// epsilon/ grid = constant
        if (gaussianRadiusType_ == "gridBasedGaussianRadius" )
        {
		coeffsEps.lookup("bladeEpsilonFactor") >> bladeEpsilonFactor;
		epsilon[0] = bladeEpsilonFactor * cellLen;
			
	}
	
	// epsilon/c = constant 
	else if (gaussianRadiusType_ == "chordBasedGaussianRadius")  
	{
		coeffsEps.lookup("bladeEpsilonFactor") >> bladeEpsilonFactor;		
		epsilon[0] = bladeEpsilonFactor * chordLength_;
	}
		
	// epsilon/c* = constant. Elliptic planform  
	else if (gaussianRadiusType_ == "ellipticGaussianRadius")
	{
		scalar c0 = 4 * avegChordLength_ / Foam::constant::mathematical::pi ;
		scalar cStar = c0 * Foam::sqrt(1 - ( Foam::sqr(rootDistance_)));
		scalar nMin = 1 ; // since n = 0 leads to numerical instabilities
		scalar epsilonRdivided2 = nMin * cellLen ;
		scalar n = 0.1 * bladeRadius_ * cellLen ; 
		scalar nMax = max(n, nMin) ;
		scalar epsilon0 = nMax * cellLen;
		epsilon[0] = max(cStar*epsilon0/c0, epsilonRdivided2);

	}
		
	else if (gaussianRadiusType_ == "chordThicknessBasedGaussian")
        {
		coeffsEps.lookup("bladeEpsilonFactorChord") >> bladeEpsilonFactorChord;
		coeffsEps.lookup("bladeEpsilonFactorThickness") >> bladeEpsilonFactorThickness;
		epsilon[0] = bladeEpsilonFactorChord * chordLength_ ;
		epsilon[1] = chordLength_* bladeEpsilonFactorThickness* thickness_ ;
		
	}
        
    }
    // Reduce epsilon over all processors
    reduce(epsilon, maxOp<vector>());


    // If epsilon is not reduced, position is not in the mesh
    if (not (mag(epsilon) > 0.0))
    {
        // Raise fatal error since mesh size cannot be detected
        FatalErrorIn("void actuatorLineElement::applyForceField()")
            << "Position of " << name_  << " not found in mesh"
            << abort(FatalError);
    }
	
	if (debug)
	{        
		Info<< "Chosen gaussianRadiusType: " << gaussianRadiusType_ << endl;
		Info<< "Calculated Epsilon" << epsilon << endl;
	}
    return epsilon;
}




void Foam::fv::actuatorLineElement::correctFlowCurvature
(
    scalar& angleOfAttackRad
)
{
    if (debug)
    {
        Info<< "    Correcting for flow curvature with "
            << flowCurvatureModelName_ << " model" << endl;
    }

    if (flowCurvatureModelName_ == "Goude")
    {
        angleOfAttackRad += omega_*chordLength_/(2*mag(relativeVelocity_));
    }
    else if (flowCurvatureModelName_ == "MandalBurton")
    {
        // Calculate relative velocity at leading and trailing edge
        vector relativeVelocityLE = inflowVelocity_ - velocityLE_;
        vector relativeVelocityTE = inflowVelocity_ - velocityTE_;

        // Calculate angle of attack at leading and trailing edge
        scalar alphaLE = asin((planformNormal_ & relativeVelocityLE)
                       / (mag(planformNormal_)*mag(relativeVelocityLE)));
        scalar alphaTE = asin((planformNormal_ & relativeVelocityTE)
                       / (mag(planformNormal_)*mag(relativeVelocityTE)));

        scalar beta = alphaTE - alphaLE;

        angleOfAttackRad += atan2((1.0 - cos(beta/2.0)), sin(beta/2.0));
    }
    else if (flowCurvatureModelName_ == "constantOffset")
    {
        dictionary fcDict = dict_.subDict("flowCurvature");
        dictionary coeffs = fcDict.subDict(flowCurvatureModelName_ + "Coeffs");
        scalar offsetDeg = 0.0;
        coeffs.lookup("offsetDeg") >> offsetDeg;
        angleOfAttackRad += degToRad(offsetDeg);
    }
}


void Foam::fv::actuatorLineElement::multiplyForceRho
(
    const volScalarField& rho
)
{
    // Lookup local density
    label cellI = findCell(position_);
    scalar localRho = VGREAT;
    if (cellI >= 0)
    {
        localRho = rho[cellI];
    }

    reduce(localRho, minOp<scalar>());
    forceVector_ *= localRho;
}


void Foam::fv::actuatorLineElement::applyForceField
(
    volVectorField& forceField
)
{
	// Calculate projection width
	vector epsilon1 = calcProjectionEpsilon();
	//Info<< "epsilon: " << epsilon1 << endl;
	
	word forceProjType_ ;
	dict_.lookup("forceProjType") >> forceProjType_; 
	                
	if (forceProjType_ == "gaussianIsotropic" || forceProjType_ == "variableGaussianIsotropicChord") 
	{
		scalar projectionRadius = (epsilon1[0] * Foam::sqrt(Foam::log(1.0/0.001)));
	        scalar sphereRadius = chordLength_ + projectionRadius;

		if (forceProjType_ == "variableGaussianIsotropicChord")
		{
			projectionRadius *= bladeChordMax;
		}
		// Apply force to the cells within the element's sphere of influence


		//scalar sphereRadius = chordLength_ + projectionRadius;
	//	sphereRadius = 63.0 + projectionRadius;		

		forAll(mesh_.cells(), cellI)
		{
			scalar dis = mag(mesh_.C()[cellI] - position_);
			                            
			if (dis <= sphereRadius)
			{
                                          
				scalar factor = Foam::exp(-Foam::sqr(dis/epsilon1[0]))
							/ (Foam::pow(epsilon1[0], 3)
							* Foam::pow(Foam::constant::mathematical::pi, 1.5));
				// forceField is opposite forceVector
				forceField[cellI] += -forceVector_*factor;
		
			}
		}
	
	/*	fileName smearing = mesh_.time().path() / "../postProcessing/smearing";
		
		if (not isDir(smearing))
		{
			mkDir(smearing);
		}
		fileName outputFile;
                outputFile = smearing / (name_+".csv");

		if (not isDir(outputFile))
		{
			outputFile = smearing / (name_+".csv");
			outputForce = new OFstream(outputFile);
		}
		
		bool isFirstWrite = !isFile(outputFile);

		std::ofstream outputForce;
		if (isFirstWrite)
		{
			outputForce.open(outputFile.c_str(), std::ios::out);  
  			outputForce << "time,root_dist,forcesmearedx,forcesmearedy,forcesmearedz" << std::endl;
		}
		else
		{
			outputForce.open(outputFile.c_str(), std::ios::app);
		}
		

		scalar currentTime = mesh_.time().value();

		volVectorField forceFieldLocal
		(
			IOobject
			(
				"forceFieldLocal",
				mesh_.time().timeName(),
				mesh_,
				IOobject::NO_READ,
				IOobject::NO_WRITE
		),
		mesh_,
		dimensionedVector("zero",forceField.dimensions(),vector::zero)
	);



		forAll(mesh_.cells(), cellI)
		{
			forceFieldLocal[cellI] = vector::zero;
		}


		forAll(mesh_.cells(), cellI)
		{
		    
		    scalar dis = mag(mesh_.C()[cellI] - position_);
		    if (dis <= sphereRadius) // Only consider cells influenced by the actuator line
		    {
			
                                scalar factor = Foam::exp(-Foam::sqr(dis/epsilon1[0]))
                                                        / (Foam::pow(epsilon1[0], 3)
                                                        * Foam::pow(Foam::constant::mathematical::pi, 1.5));
                                // forceField is opposite forceVector
                                forceFieldLocal[cellI] += -forceVector_*factor;

	
	  	    }
		    
		}

		vector summedForce = vector::zero;

		forAll(mesh_.cells(), cellI)
		{
			summedForce += forceFieldLocal[cellI];
		}
		
		vector globalSummedForce = summedForce;
		reduce(globalSummedForce, sumOp<vector>());

		if (Pstream::master())
		{

			outputForce << currentTime << "," << rootDistance_ << ","
				<< globalSummedForce.x()* Foam::pow(cellLength_,3) << "," << globalSummedForce.y()*Foam::pow(cellLength_,3)
				<< "," << globalSummedForce.z()*Foam::pow(cellLength_,3) << std::endl;
		}

		outputForce.close(); */

	}

	
	// Compute 3D gaussian that has different spreading in each direction(chordwise, thickwise, radial direction)
	// This is 2D gaussian. Since it has no optimum value for radial 
	else if (forceProjType_ == "gaussianNonIsotropic")
	{
		scalar epsilonMax = max(epsilon1[0], epsilon1[1]);
		scalar projectionRadius = (bladeChordMax* epsilonMax * Foam::sqrt(Foam::log(1.0/0.001)));
		scalar sphereRadius = chordLength_ + projectionRadius;

		//vector sphereRadiusVector = vector(0.0, 0.0, 0.0);
		//sphereRadiusVector[0] = chordLength_ + projectionRadius;		
		//sphereRadiusVector[1] = (chordLength_* thickness_)+ projectionRadius;
		//sphereRadiusVector[2] = spanLength_+ projectionRadius;

		
		forAll(mesh_.cells(), cellI)
		{
			vector disVector = mesh_.C()[cellI] - position_;

			vector dirC = chordDirection_; //chordwise(1)
		        vector dirT = planformNormal_; //thickwise(0)
		        vector dirR = spanDirection_; //radial direction (2)
			
			// There is no need to rotate the direction depending on the pitch.
			// Code performs the pitching already before here
			rotateVector(dirC, vector::zero, dirR, -degToRad(bladeElementPitch_+angleOfAttack_));
			rotateVector(dirT, vector::zero, dirR, -degToRad(bladeElementPitch_+angleOfAttack_));

			scalar dir0 = disVector & dirC;
			scalar dir1 = disVector & dirT;

			if (mag(disVector) <= sphereRadius) 
			{						
				scalar c =  (1.0 / (epsilon1[0]*epsilon1[1]*Foam::constant::mathematical::pi));
				scalar g = Foam::exp( - Foam::sqr(dir0/epsilon1[0]) - Foam::sqr(dir1/epsilon1[1]) );
				scalar factor = c*g;					
				
			// forceField is opposite forceVector
				forceField[cellI] += -forceVector_*factor;

			}
		
		}
	}
	
/*    if (debug)
    {
        Info<< "    sphereRadius: " << sphereRadius << endl;
    }*/
}

void Foam::fv::actuatorLineElement::calculateInflowVelocity
(
    const volVectorField& Uin
)
{
    // Find local flow velocity by interpolating to element location
    inflowVelocity_ = vector(VGREAT, VGREAT, VGREAT);
    vector inflowVelocityPoint = position_;
    interpolationCellPoint<vector> UInterp(Uin);

    // velEvalType_ already read once in read() — no per-call dict lookup needed

    if (velEvalType_ == "cellCenter")
    {
        label inflowCellI = findCell(inflowVelocityPoint);
        if (inflowCellI >= 0)
        {
            inflowVelocity_ = UInterp.interpolate(inflowVelocityPoint, inflowCellI);
        }

        // Guard against invalid values before reduce
        if (inflowVelocity_[0] >= VGREAT)
        {
            inflowVelocity_ = vector(GREAT, GREAT, GREAT);
        }

        reduce(inflowVelocity_, minOp<vector>()); // Note: still lexicographical min
    }
    else if (velEvalType_ == "integral")
    {
        vector epsilon = calcProjectionEpsilon();
        scalar sampleRadius = epsilon[0] * velocitySampleRadius_;
        vector chordNormal = chordDirection_ / mag(chordDirection_);

        vector velocitySum = vector::zero;
        int validSampleCount = 0;

        for (label point = 0; point < nVelocitySamples_; point++)
        {
            scalar pointAngle = constant::mathematical::pi * 2.0 * point / nVelocitySamples_;
            scalar chordDist = sampleRadius * cos(pointAngle);
            scalar normalDist = sampleRadius * sin(pointAngle);
            vector samplePoint = inflowVelocityPoint + chordDist * chordNormal + normalDist * planformNormal_;

            vector sampleVelocity = vector(GREAT, GREAT, GREAT);
            label sampleCellI = findCell(samplePoint);
            if (sampleCellI >= 0)
            {
                sampleVelocity = UInterp.interpolate(samplePoint, sampleCellI);
            }

            // Validate before reduction
            if (sampleVelocity[0] < VGREAT)
            {
                velocitySum += sampleVelocity;
                validSampleCount++;
            }
        }

        reduce(velocitySum, sumOp<vector>());
        reduce(validSampleCount, sumOp<int>());

        if (validSampleCount > 0)
        {
            inflowVelocity_ = velocitySum / validSampleCount;
        }
        else
        {
            FatalErrorIn("void actuatorLineElement::calculateInflowVelocity")
                << "No valid inflow velocity samples for " << name_ << abort(FatalError);
        }
    }
    else if (velEvalType_ == "EVM")
    {
        vector cellVelocity = vector(GREAT, GREAT, GREAT);
        vector cellCenterI = vector(GREAT, GREAT, GREAT);
        label airfoilCellID = findCell(inflowVelocityPoint);

        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
        const scalarField& V = mesh_.V();

        scalar cellLengthLocal = GREAT; // Local fallback value

        if (airfoilCellID >= 0)
        {
            cellVelocity = U[airfoilCellID];
            cellCenterI = mesh_.C()[airfoilCellID];
            cellLengthLocal = cbrt(V[airfoilCellID]);
        }

        // Replacing old reduce logic with guarded values
        // reduce(cellVelocity, minOp<vector>());
        // reduce(cellCenterI, minOp<vector>());
        // reduce(cellLength_, minOp<scalar>());

        reduce(airfoilCellID, maxOp<label>());
        reduce(cellVelocity, minOp<vector>()); // Still may want better logic
        reduce(cellCenterI, minOp<vector>());
        reduce(cellLengthLocal, minOp<scalar>());

        cellLength_ = cellLengthLocal;

        if (airfoilCellID >= 0)
        {
            vector relVel_ = cellVelocity - velocity_;
            scalar relVelMag = mag(relVel_);

            if (relVelMag < SMALL)
            {
                // Near-zero relative velocity: sampling direction undefined,
                // fall back to last valid inflow rather than divide by zero.
                if (Pstream::master())
                {
                    Warning << "Near-zero relative velocity for " << name_
                            << " in EVM sampling. Using previous inflow velocity." << endl;
                }
                inflowVelocity_ = prevInflowVelocity_;
            }
            else
            {
                vector EVMsampDir = spanDirection_ ^ (relVel_ / relVelMag);
                EVMsampDir /= mag(EVMsampDir);

                // Use cached EVM parameters (read once in read(), not per timestep)
                scalar distance       = evmDistance_;
                scalar samplingLength = evmSamplingLength_;
                label  numEVMPoints   = evmNumPoints_;   // label, not scalar

                scalar EVMInt = samplingLength / (numEVMPoints - 1);
                vector velocityEVM = vector::zero;
                int nSample = 0;

                for (int m = 0; m < numEVMPoints; m++)
                {
                    point EVMPoint = position_ - distance * cellLength_ * (relVel_ / relVelMag)
                                     + (samplingLength / 2 - EVMInt * m) * cellLength_ * EVMsampDir;

                    label EVMCellID = findCell(EVMPoint);
                    if (EVMCellID >= 0)
                    {
                        //velocityEVM += U[EVMCellID];
                        velocityEVM += UInterp.interpolate(EVMPoint, EVMCellID);
                        nSample++;
                    }
                }

                reduce(velocityEVM, sumOp<vector>());
                reduce(nSample, sumOp<int>());

                if (nSample > 0)
                {
                    inflowVelocity_ = velocityEVM / nSample;
                }
                else
                {
                    // All EVM sample points missed every processor's subdomain.
                    // Fall back to the cell-centre velocity (already globally reduced).
                    if (Pstream::master())
                    {
                        Warning << "EVM sample points all missed domain for " << name_
                                << ". Falling back to cell-centre velocity." << endl;
                    }
                    inflowVelocity_ = cellVelocity;
                }
            }
        }
        else
        {
            // Element not found on any processor — transient domain-boundary
            // situation not resolved even by contains(). Use the last valid
            // inflow rather than the far-field value to minimise prediction error.
            if (Pstream::master())
            {
                Warning << "Element " << name_
                        << " not found in any processor domain during EVM sampling."
                        << " Using previous inflow velocity as fallback." << endl;
            }
            inflowVelocity_ = prevInflowVelocity_;
        }
    }

    // Final validation — should now only trigger for genuinely unhandled velEvalType
    if (inflowVelocity_[0] >= VGREAT)
    {
        FatalErrorIn("void actuatorLineElement::calculateInflowVelocity")
            << "Final inflow velocity still unset for " << name_
            << ". velEvalType = " << velEvalType_
            << " — check that this is a valid type (cellCenter, integral, EVM)."
            << abort(FatalError);
    }

    // Cache for use as fallback in subsequent timesteps
    prevInflowVelocity_ = inflowVelocity_;

    reduce(inflowVelocity_, minOp<vector>()); // Could still use better logic here
}


void Foam::fv::actuatorLineElement::createOutputFile()
{
    fileName dir;

    if (Pstream::parRun())
    {
        dir = mesh_.time().path()/"../postProcessing/actuatorLineElements"
            / mesh_.time().timeName();
    }
    else
    {
        dir = mesh_.time().path()/"postProcessing/actuatorLineElements"
            / mesh_.time().timeName();
    }

    if (not isDir(dir))
    {
        mkDir(dir);
    }

    outputFile_ = new OFstream(dir/name_ + ".csv");

    *outputFile_<< "time,root_dist,x,y,z,rel_vel_mag,inflow_X,inflow_Y,inflow_Z,Re,alpha_deg,alpha_uncorrected,"
                << "alpha_geom_deg,cl,cd,fx,fy,fz,end_effect_factor,"
                << "c_ref_t,c_ref_n,f_ref_t,f_ref_n, FL*, FD*, c/l, relVel_X,relVel_Y,relVel_Z, clUn, cdUn" << endl;
}


void Foam::fv::actuatorLineElement::writePerf()
{
    scalar time = mesh_.time().value();

    // write time,root_dist,x,y,z,rel_vel_mag,Re,alpha_deg,alpha_geom_deg,cl,cd,
    // fx,fy,fz,end_effect_factor,c_ref_t,c_ref_n,f_ref_t,f_ref_n, FL*, FD*
    *outputFile_<< time << "," << rootDistance_ << "," << position_.x() << ","
                << position_.y() << "," << position_.z() << ","
                << mag(relativeVelocity_) << "," << inflowVelocity_.x() << "," << inflowVelocity_.y() << "," 
		<< inflowVelocity_.z() << ","  << Re_ << "," << angleOfAttack_ << ","
                << aoaUncorrected << "," << angleOfAttackGeom_ << "," << liftCoefficient_ << ","
                << dragCoefficient_ << "," << forceVector_.x() << ","
                << forceVector_.y() << "," << forceVector_.z() << ","
                << endEffectFactor_ << "," << tangentialRefCoefficient() << ","
                << normalRefCoefficient() << "," << tangentialRefForce() << ","
                << normalRefForce() << "," <<  normalizedLift << "," << normalizedDrag << "," <<  chordLength_/cellLength_ << ","
		<< relativeVelocity_.x() << "," << relativeVelocity_.y() << "," << relativeVelocity_.z() << ","
		<< liftCoefficientUncorrected << "," << dragCoefficientUncorrected  <<endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineElement::actuatorLineElement
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    name_(name),
    mesh_(mesh),
    meshBoundBox_(mesh_.points(), false),
    planformNormal_(vector::zero),
    velocity_(vector::zero),
    forceVector_(vector::zero),
    relativeVelocity_(vector::zero),
    relativeVelocityGeom_(vector::zero),
    inflowVelocity_(vector::zero),
    prevInflowVelocity_(dict.lookup("freeStreamVelocity")),
    angleOfAttack_(0.0),
    angleOfAttackGeom_(0.0),
    liftCoefficient_(0.0),
    dragCoefficient_(0.0),
    momentCoefficient_(0.0),
    profileName_(dict.lookup("profileName")),
    profileData_(profileName_, dict.subDict("profileData"), debug),
    dynamicStallActive_(false),
    filteredLiftingLineActive_(false),
    filteredLiftingLineWriteOutput_(false),
    omega_(0.0),
    chordMount_(0.25),
    flowCurvatureActive_(false),
    flowCurvatureModelName_("none"),
    velocityLE_(vector::zero),
    velocityTE_(vector::zero),
    writePerf_(false),
    rootDistance_(0.0),
    endEffectFactor_(1.0),
    addedMassActive_(dict.lookupOrDefault("addedMass", false)),
    addedMass_(mesh.time(), dict.lookupOrDefault("chordLength", 1.0), debug),
    prevLiftCoefficientG_(0.0),
    prevTime_(mesh_.time().value()),
    nextGPrevTime_(0.0),
    prevGPrevTime_(0.0),
    G_(0.0),
    isFirstCallInTimeStep_(true),
    prevBladePointPerturbation_(0.0)
{
    meshBoundBox_.inflate(1e-6);
    read();
    //inflowVelocityGList_.setSize(bladePointRadius.size(), vector::zero);
    setPrevG(0.0);
    setNextG(0.0);

    if (writePerf_)
    {
        createOutputFile();
    }
}

// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineElement::~actuatorLineElement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::fv::actuatorLineElement::name() const
{
    return name_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::chordLength() const
{
    return chordLength_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::spanLength()
{
    return spanLength_;
}


const Foam::vector& Foam::fv::actuatorLineElement::position()
{
    return position_;
}


const Foam::vector& Foam::fv::actuatorLineElement::velocity()
{
    return velocity_;
}


const Foam::vector& Foam::fv::actuatorLineElement::relativeVelocity()
{
    return relativeVelocity_;
}


const Foam::vector& Foam::fv::actuatorLineElement::relativeVelocityGeom()
{
    return relativeVelocityGeom_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::angleOfAttack()
{
    return angleOfAttack_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::angleOfAttackGeom()
{
    return angleOfAttackGeom_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::liftCoefficient()
{
    return liftCoefficient_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::dragCoefficient()
{
    return dragCoefficient_;
}



const Foam::scalar& Foam::fv::actuatorLineElement::momentCoefficient()
{
    return momentCoefficient_;
}


Foam::scalar Foam::fv::actuatorLineElement::tangentialRefCoefficient()
{
    return profileData_.convertToCRT
    (
        liftCoefficient_,
        dragCoefficient_,
        calcInflowRefAngle()
    );
}


void Foam::fv::actuatorLineElement::setPrevG(scalar prevG)
{
    prevG_ = prevG;
}


void Foam::fv::actuatorLineElement::setNextG(scalar nextG)
{
    nextG_ = nextG;
}

Foam::scalar Foam::fv::actuatorLineElement::getNextG() const
{
    return nextG_;
}

Foam::scalar Foam::fv::actuatorLineElement::getPrevG() const
{
    return prevG_;
}


void Foam::fv::actuatorLineElement::setNextGPrevTime(scalar nextGPrevTime)
{
    nextGPrevTime_ = nextGPrevTime;
}


void Foam::fv::actuatorLineElement::setPrevGPrevTime(scalar prevGPrevTime)
{
    prevGPrevTime_ = prevGPrevTime;
}

Foam::scalar Foam::fv::actuatorLineElement::getNextGPrevTime() const
{ 
    return nextGPrevTime_; 
}

Foam::scalar Foam::fv::actuatorLineElement::getPrevGPrevTime() const
{ 
    return prevGPrevTime_; 
}

void Foam::fv::actuatorLineElement::setPrevBladePointPerturbation(const scalar perturbation)
{
     prevBladePointPerturbation_ = perturbation;
}

Foam::scalar Foam::fv::actuatorLineElement::getPrevBladePointPerturbation() const
{
    return prevBladePointPerturbation_; 
}


void Foam::fv::actuatorLineElement::setPrevLiftCoefficient(const scalar liftcoeff)
{
     prevLiftCoefficientG_ = liftcoeff;
}

Foam::scalar Foam::fv::actuatorLineElement::getPrevLiftCoefficient() const
{
    return prevLiftCoefficientG_;
}


Foam::scalar Foam::fv::actuatorLineElement::tangentialRefForce()
{
    return 0.5 * chordLength_ * tangentialRefCoefficient()
        * magSqr(relativeVelocity_);
}


Foam::scalar Foam::fv::actuatorLineElement::normalRefCoefficient()
{
    return profileData_.convertToCRN
    (
        liftCoefficient_,
        dragCoefficient_,
        calcInflowRefAngle()
    );
}


Foam::scalar Foam::fv::actuatorLineElement::normalRefForce()
{
    return 0.5 * chordLength_ * normalRefCoefficient()
        * magSqr(relativeVelocity_);
}


Foam::scalar Foam::fv::actuatorLineElement::calcInflowRefAngle()
{
    // Calculate inflow velocity angle in degrees (phi angle)
    
    scalar inflowVelAngleRad = acos
    (
        (-relativeVelocity_ & chordRefDirection_)
        / (mag(relativeVelocity_) * mag(chordRefDirection_))
    );
    
    return radToDeg(inflowVelAngleRad);
}


const Foam::scalar& Foam::fv::actuatorLineElement::rootDistance()
{
    return rootDistance_;
}



/*const Foam::List<scalar>& Foam::fv::actuatorLineElement::getDGList() const
{
    return dGList_;
}*/


Foam::scalar Foam::fv::actuatorLineElement::calcAngleOfAttackUncorrected()
{
    // local angle of attack =  phi - pitch
    scalar angleOfAttackUncorrected;
    scalar asinArg = (planformNormal_ & relativeVelocity_) / (mag(planformNormal_) * mag(relativeVelocity_));

    // Check if asinArg is within the valid range [-1, 1]
    // NaN check must come first: NaN comparisons always return false,
    // so a NaN asinArg would silently pass both bounds checks and reach asin(NaN).
    if (std::isnan(asinArg))
    {
        if (Pstream::master())
        {
            Warning << "asin() argument is NaN for " << name_
                    << " — planformNormal=" << planformNormal_
                    << ", relativeVelocity=" << relativeVelocity_
                    << ". Clamping to 0." << endl;
        }
        asinArg = 0.0;
    }
    else if (asinArg > 1.0)
    {
        if (Pstream::master())
        {
            scalar excess = asinArg - 1.0;
            if (excess < 1e-6)
            {
                // Tiny round-off: only log at debug level to avoid noise
                if (debug)
                {
                    Info << "asin() argument clamped from " << asinArg
                         << " to 1.0 (round-off) for " << name_ << endl;
                }
            }
            else
            {
                Warning << "asin() argument " << asinArg
                        << " exceeds 1.0 by " << excess << " for " << name_
                        << " — check planformNormal and relativeVelocity." << endl;
            }
        }
        asinArg = 1.0;
    }
    else if (asinArg < -1.0)
    {
        if (Pstream::master())
        {
            scalar excess = -1.0 - asinArg;
            if (excess < 1e-6)
            {
                if (debug)
                {
                    Info << "asin() argument clamped from " << asinArg
                         << " to -1.0 (round-off) for " << name_ << endl;
                }
            }
            else
            {
                Warning << "asin() argument " << asinArg
                        << " is below -1.0 by " << excess << " for " << name_
                        << " — check planformNormal and relativeVelocity." << endl;
            }
        }
        asinArg = -1.0;
    }

    angleOfAttackUncorrected = asin(asinArg);
    lookupCoefficients(radToDeg(angleOfAttackUncorrected));
    liftCoefficientUncorrected = liftCoefficient_;
    dragCoefficientUncorrected = dragCoefficient_;
    
    return radToDeg(angleOfAttackUncorrected);
}



// applying the EVM angle of attack correction with robust anti-divergence measures
// and selective underrelaxation when convergence stalls
Foam::scalar Foam::fv::actuatorLineElement::calcAngleOfAttackCorrected()
{
    // Initialize the variables for EVM
    scalar tolerance = 1e-4;
    scalar angleOfAttackUncorrected = calcAngleOfAttackUncorrected();
    if (debug)
    {
        Info << "EVM AoA correction starting, uncorrected AoA: "
             << angleOfAttackUncorrected << " deg" << endl;
    }
    
    // Use fixed-point iteration for first attempt (original method)
    scalar angleOfAttackCorrected = angleOfAttackUncorrected;
    scalar angleOfAttackDifference = VGREAT;
    label counter = 0;
    scalar bestAngle = angleOfAttackUncorrected;
    scalar bestDifference = VGREAT;
    scalar relaxationFactor = 0.5; // Initial relaxation factor
    scalar previousDifference = VGREAT;
    scalar stallAngle = 12.0;
    label  divergenceCount = 0;    // consecutive increases needed before switching
    
    // First try pure fixed-point iteration (the original method) for several iterations
    while (angleOfAttackDifference > tolerance && counter < 100)
    {
        scalar oldAngleOfAttackCorrected = angleOfAttackCorrected;
        
        // Lookup the airfoil tabular data
        lookupCoefficients(angleOfAttackCorrected);
        
        // Calculate correction term
        scalar dAlfa = (chordLength_ / cellLength_) * liftCoefficient_ * (1.2553 - (0.0552 * dragCoefficient_));
        
        // Adaptive relaxation near stall <----- MODIFIED
        if (abs(angleOfAttackCorrected - stallAngle) < 5.0) {
            relaxationFactor = 0.2;
        } else {
            relaxationFactor = 0.5;
        }
        
        // Apply under-relaxation
        angleOfAttackCorrected = relaxationFactor * (angleOfAttackUncorrected - dAlfa) + (1 - relaxationFactor) * oldAngleOfAttackCorrected;
        
        // Calculate the difference
        angleOfAttackDifference = mag(angleOfAttackCorrected - oldAngleOfAttackCorrected);
        
        // Track best solution so far
        scalar currentError = mag(angleOfAttackCorrected - (angleOfAttackUncorrected - dAlfa));
        if (currentError < bestDifference) {
            bestDifference = currentError;
            bestAngle = angleOfAttackCorrected;
        }
        
        // Increment counter
        counter++;
        
        // Debug info every 5 iterations
        if (debug && label(counter) % 5 == 0)
        {
            Info << "FPI Iteration " << counter 
                 << ", Current AOA: " << angleOfAttackCorrected
                 << ", Diff: " << angleOfAttackDifference << endl;
        }
        
        // Require 3 consecutive increases before declaring divergence and
        // switching to bisection — a single upward blip under relaxation is
        // not sufficient evidence of true divergence.
        if (counter > 5)
        {
            if (angleOfAttackDifference > previousDifference)
            {
                divergenceCount++;
                if (divergenceCount >= 3)
                {
                    if (debug)
                    {
                        Info << "FPI diverging for 3 consecutive iterations, "
                             << "switching to bisection method..." << endl;
                    }
                    break;
                }
            }
            else
            {
                divergenceCount = 0; // reset on any improvement
            }
        }
        
        previousDifference = angleOfAttackDifference;
        
        // Check if we're converging well with just fixed-point iteration
        if (angleOfAttackDifference < tolerance) {
            if (debug)
            {
                Info << "Converged with FPI in " << counter << " iterations"
                     << ", Final AOA: " << angleOfAttackCorrected << " deg" << endl;
            }
            return angleOfAttackCorrected;
        }
    }
    
    // If we reach here, try a more robust bisection method
    if (debug)
    {
        Info << "FPI did not converge after " << counter
             << " iterations, switching to bisection." << endl;
    }
    
    // Use best value from fixed-point as starting point
    angleOfAttackCorrected = bestAngle;
    
    // Bisection method parameters
    scalar lowerBound = max(-5.0, angleOfAttackUncorrected - 15.0);
    scalar upperBound = min(25.0, angleOfAttackUncorrected + 15.0);
    label bisectionCounter = 0;
    scalar bisectionTolerance = 1e-3;
    
    // Evaluate function at bounds to ensure we bracket a root
    lookupCoefficients(lowerBound);
    scalar lowerDAlfa = (chordLength_ / cellLength_) * liftCoefficient_ * (1.2553 - (0.0552 * dragCoefficient_));
    scalar lowerFunc = lowerBound - (angleOfAttackUncorrected - lowerDAlfa);
    
    lookupCoefficients(upperBound);
    scalar upperDAlfa = (chordLength_ / cellLength_) * liftCoefficient_ * (1.2553 - (0.0552 * dragCoefficient_));
    scalar upperFunc = upperBound - (angleOfAttackUncorrected - upperDAlfa);
    
    // Check if we have a proper bracket (opposite signs)
    if (sign(lowerFunc) == sign(upperFunc)) {
        // If not bracketed, adjust bounds to ensure opposite signs
        if (debug)
        {
            Info << "Bisection bracket invalid, expanding bounds..." << endl;
        }
        
        // Try expanding the interval
        for (label i = 0; i < 5; i++) {
            lowerBound = max(-10.0, lowerBound - 5.0);
            upperBound = min(30.0, upperBound + 5.0);
            
            lookupCoefficients(lowerBound);
            lowerDAlfa = (chordLength_ / cellLength_) * liftCoefficient_ * (1.2553 - (0.0552 * dragCoefficient_));
            lowerFunc = lowerBound - (angleOfAttackUncorrected - lowerDAlfa);
            
            lookupCoefficients(upperBound);
            upperDAlfa = (chordLength_ / cellLength_) * liftCoefficient_ * (1.2553 - (0.0552 * dragCoefficient_));
            upperFunc = upperBound - (angleOfAttackUncorrected - upperDAlfa);
            
            if (sign(lowerFunc) != sign(upperFunc)) {
                break;  // Found a bracket
            }
        }
    }
    
    // If still not bracketed, fall back to best approximation from fixed-point
    if (sign(lowerFunc) == sign(upperFunc)) {
        if (Pstream::master())
        {
            Warning
                << "Could not bracket EVM AoA correction for " << name_
                << ". Using best fixed-point approximation: " << bestAngle
                << " deg" << endl;
        }
        return bestAngle;
    }
    
    // Bisection method main loop
    while (mag(upperBound - lowerBound) > bisectionTolerance && bisectionCounter < 30)
    {
        // Calculate midpoint
        scalar midPoint = 0.5 * (lowerBound + upperBound);
        
        // Evaluate function at midpoint
        lookupCoefficients(midPoint);
        scalar midDAlfa = (chordLength_ / cellLength_) * liftCoefficient_ * (1.2553 - (0.0552 * dragCoefficient_));
        scalar midFunc = midPoint - (angleOfAttackUncorrected - midDAlfa);
        
        // Bisection: narrow the bracket — no underrelaxation.
        // The bracket itself guarantees convergence; blending midPoint
        // with a previous iterate would corrupt the returned value.
        angleOfAttackCorrected = midPoint;
        
        // Update bounds based on sign
        if (sign(midFunc) == sign(lowerFunc)) {
            lowerBound = midPoint;
            lowerFunc = midFunc;
        } else {
            upperBound = midPoint;
            upperFunc = midFunc;
        }
        
        // Increment counter
        bisectionCounter++;
        
        // Debug info
        if (debug && label(bisectionCounter) % 3 == 0)
        {
            Info << "Bisection Iteration " << bisectionCounter 
                 << ", Current AOA: " << angleOfAttackCorrected
                 << ", Interval: [" << lowerBound << ", " << upperBound << "]" << endl;
        }
    }
    
    // Return the bracket midpoint as the converged answer
    angleOfAttackCorrected = 0.5 * (lowerBound + upperBound);
    
    if (debug)
    {
        Info << "Bisection converged in " << bisectionCounter << " iterations"
             << ", Final AOA: " << angleOfAttackCorrected << " deg" << endl;
    }
    return angleOfAttackCorrected;
}


Foam::List<vector> Foam::fv::actuatorLineElement::relativeVelocityGList_;

void Foam::fv::actuatorLineElement::setRelativeVelocityGList
(
    const volVectorField& Uin
)
{
    // Calculate vector normal to chord--span plane
    planformNormal_ = -chordDirection_ ^ spanDirection_;
    planformNormal_ /= mag(planformNormal_);

    // Find local flow velocity by interpolating to element location
    calculateInflowVelocity(Uin);

    // Subtract spanwise component of inflow velocity
    vector spanwiseVelocity = spanDirection_
                            * (inflowVelocity_ & spanDirection_)
                            / magSqr(spanDirection_);

   inflowVelocityG_ = inflowVelocity_ - spanwiseVelocity;

   if (numOfElement == 0)
   {
	relativeVelocityGList_.clear();
        relativeVelocityGList_.setSize(bladePointRadius.size(), vector::zero);
   }
        

   //inflowVelocityGList_[numOfElement] = inflowVelocityG_;

    // Calculate relative velocity and Reynolds number
    relativeVelocity_ = inflowVelocityG_ - velocity_;

   relativeVelocityGList_[numOfElement] = relativeVelocity_ ;

   scalar angleOfAttackBeforeInducedVel =  asin((planformNormal_ & relativeVelocity_)
                            / (mag(planformNormal_)
                           *  mag(relativeVelocity_)));

   lookupCoefficients(radToDeg(angleOfAttackBeforeInducedVel));
	
   liftCoefficientG_ = liftCoefficient_;
   //Info << "liftCoefficientG:" << liftCoefficientG_ << endl;
   //Info << "liftCoefficientGold:" << prevLiftCoefficientG_ << endl;

   if (numOfElement != 0)
   {
        liftCoefficientG_ = (relaxFactorLC * liftCoefficientG_) + ( (1- relaxFactorLC) * prevLiftCoefficientG_);

   }
   prevLiftCoefficientG_ = liftCoefficientG_;

}

Foam::scalar Foam::fv::actuatorLineElement::computeG()
{
   scalar Gnew = (0.5 * liftCoefficientG_ * mag(relativeVelocity_) * mag(relativeVelocity_)* chordLength_) ;

   if (mesh_.time().value() != prevTime_)
   {
	prevTime_ = mesh_.time().value();
        isFirstCallInTimeStep_ = true; 
   }

   
   if (prevTime_ > 0 && isFirstCallInTimeStep_) // First timestep 
   {
	GPrevTime_ = G_;
	isFirstCallInTimeStep_ = false;
   }
	

   G_ = Gnew;
   
   return G_;
}

Foam::List<scalar> Foam::fv::actuatorLineElement:: dGList_;

void Foam::fv::actuatorLineElement::computeDGandAppend()
{
   label totNumberOfElement = bladePointRadius.size();

   scalar dG;

   if (numOfElement == 0)
   {
	dGList_.clear();
        dGList_.setSize(totNumberOfElement);
	//prevTime_ = mesh_.time().value();
   }


   if (numOfElement == 0)
   {
        dG = GPrevTime_;
   }
   else if (numOfElement == totNumberOfElement-1)
   {
        dG = -GPrevTime_;
   }
   else
   {
        dG = 0.5 * (nextGPrevTime_ - prevGPrevTime_);
        //Info << "PreviousG:" << prevG_ <<endl;
        //Info << "NextG:" << nextG_ << endl;
        //Info << "PreviousGPrevTime:" << prevGPrevTime_ <<endl;
        //Info << "NextGPrevTime:" << nextGPrevTime_ << endl;
   }

   //Info << "CurrentG:" << G_ << endl;
   //Info << "PrevTimeG:" << GPrevTime_ << endl;

   dGList_[numOfElement] = dG;
   
}


Foam::vector Foam::fv::actuatorLineElement::computeInducedVelocity()
{
   label numBladePoints = bladePointRadius.size();
   scalar epsThickness;
 
   if (debug)
   {
       Info<< "computeInducedVelocity: element index " << numOfElement << endl;
   }
   
   if (velEvalType_ == "integral")
   {
   	epsThickness = epsilon[0];
   }
   else 
   {
   	vector epsilonVector = calcProjectionEpsilon();
	epsThickness = epsilonVector[0];
   }


   scalar sumLES = 0;
   scalar sumOpt = 0;
   for (label i= 0; i < numBladePoints;  i++)
   {
	
	if (numOfElement != i)
	{
		sumLES += (-1/mag(relativeVelocityGList_[i])) * dGList_[i] * (1/ (4* Foam::constant::mathematical::pi * (-bladePointRadius[i] + bladePointRadius[numOfElement])) * ( 1 - (Foam::exp(-(Foam::pow((bladePointRadius[numOfElement] - bladePointRadius[i]), 2)) / ( Foam::pow((cellLength_),2))))));

                sumOpt += (-1/mag(relativeVelocityGList_[i])) * dGList_[i] * (1/ (4* Foam::constant::mathematical::pi * (-bladePointRadius[i] + bladePointRadius[numOfElement])) * ( 1 - (Foam::exp(-(Foam::pow((bladePointRadius[numOfElement] - bladePointRadius[i]), 2))  / ( Foam::pow((bladeRadius_ * 1 /32 ),2))))));

	}
	else
	{
		sumLES += 0;
		sumOpt += 0;
	}		

   } 
   scalar bladePointPerturbation = sumOpt - sumLES ;
  
   if (numOfElement != 0)
   { 
	bladePointPerturbation = (relaxFactorIV * bladePointPerturbation) + ( (1 - relaxFactorIV) * prevBladePointPerturbation_);
   }

   prevBladePointPerturbation_ = sumOpt - sumLES;

   // Vector for the correction is computed 

   //vector notInfluencedVelocityParallelVector = chordDirection_ * inflowVelocity_.x() + planformNormal_ * inflowVelocity_.y();
   vector notInfluencedVelocityParallelVector = relativeVelocityGList_[numOfElement] /  mag(relativeVelocityGList_[numOfElement]);
   vector bladePointPerturbationVector = notInfluencedVelocityParallelVector ^ spanDirection_;
   bladePointPerturbationVector = bladePointPerturbation * bladePointPerturbationVector ; 
 
	
   if (filteredLiftingLineWriteOutput_)
   {

	fileName filteredLiftingLineOutput = mesh_.time().path() / "../postProcessing/filteredLiftingLineOutput";

        if (not isDir(filteredLiftingLineOutput))
        {
		mkDir(filteredLiftingLineOutput);
        }
        
        fileName outputFileFll;
        outputFileFll = filteredLiftingLineOutput / (name_ +".csv");

        bool isFirstWrite = !isFile(outputFileFll);

        std::ofstream outputFll;
	
	if (isFirstWrite)
        {
		outputFll.open(outputFileFll.c_str(), std::ios::out);
		outputFll << "time, root_dist,dG, bladePerturbationX, bladePerturbationY, bladePerturbationZ" << std::endl;
	}
               
	else
        { 
        	outputFll.open(outputFileFll.c_str(), std::ios::app);
        }


	scalar currTime = mesh_.time().value();

	if (Pstream::master())
	{

		outputFll << currTime << "," << rootDistance_ << ","
                            << dGList_[numOfElement] << "," <<  bladePointPerturbationVector.x() << "," << bladePointPerturbationVector.y()
                            << "," << bladePointPerturbationVector.z() << std::endl;
        }

	outputFll.close();

   }

    
   return bladePointPerturbationVector;	 

}

void Foam::fv::actuatorLineElement::calculateForce
(
    const volVectorField& Uin
)
{
    scalar pi = Foam::constant::mathematical::pi;
    
    // Calculate vector normal to chord--span plane
    planformNormal_ = -chordDirection_ ^ spanDirection_;
    planformNormal_ /= mag(planformNormal_);

    if (debug)
    {
        Info<< "Calculating force contribution from actuatorLineElement "
            << name_ << endl;
        Info<< "    position: " << position_ << endl;
        Info<< "    chordDirection: " << chordDirection_ << endl;
        Info<< "    spanDirection: " << spanDirection_ << endl;
        Info<< "    elementVelocity: " << velocity_ << endl;
        Info<< "    planformNormal: " << planformNormal_ << endl;
    }

        // Find local flow velocity by interpolating to element location
        //calculateInflowVelocity(Uin);
    
    //word velEval = dict_.lookup("velEvalType");
    
    if (filteredLiftingLineActive_)
    {
        vector inducedVelocity = computeInducedVelocity();
	//Info << "inducedVel:" << inducedVelocity << endl;
        relativeVelocity_ = relativeVelocity_ + inducedVelocity;
	//Info << "relativeVelocity_:" << relativeVelocity_ << endl;
    }
    else 
    {
	// Find local flow velocity by interpolating to element location
	calculateInflowVelocity(Uin);

        // Subtract spanwise component of inflow velocity
        vector spanwiseVelocity = spanDirection_
                            * (inflowVelocity_ & spanDirection_)
                            / magSqr(spanDirection_);
        //Info << "inflowVelocity:" << inflowVelocity_ ;
        inflowVelocity_ -= spanwiseVelocity;
	// Calculate relative velocity and Reynolds number
	relativeVelocity_ = inflowVelocity_ - velocity_;

    }


    aoaUncorrected = calcAngleOfAttackUncorrected();

    // Calculate angle of attack 
    if (velEvalType_ == "EVM")
    {
		angleOfAttack_ = calcAngleOfAttackCorrected(); // in degrees
		// rotate wind vector according to correction
		scalar differenceAoA = angleOfAttack_ - aoaUncorrected;
		rotateVector(relativeVelocity_, vector::zero, spanDirection_, degToRad(-differenceAoA));
	
    }
    else
    {
		angleOfAttack_ = aoaUncorrected; // in degrees
    }


    //angleOfAttack_ = aoaUncorrected;

    Re_ = mag(relativeVelocity_)*chordLength_/nu_;	
    // Calculate geometric angle of attack
    relativeVelocityGeom_ = freeStreamVelocity_ - velocity_;
    angleOfAttackGeom_ = asin((planformNormal_ & relativeVelocityGeom_)
                       / (mag(planformNormal_)*mag(relativeVelocityGeom_)));
    angleOfAttackGeom_ *= 180.0/pi;

    angleOfAttackRad = degToRad(angleOfAttack_);
    
    // Apply flow curvature correction to angle of attack
    if (flowCurvatureActive_)
    {
        correctFlowCurvature(angleOfAttackRad);
    }

    // Update Reynolds number of profile data
    profileData_.updateRe(Re_);

    // Lookup lift and drag coefficients
    lookupCoefficients(angleOfAttack_);



    if (debug)
    {
        Info<< "    inflowVelocity: " << inflowVelocity_ << endl;
        Info<< "    relativeVelocity: " << relativeVelocity_ << endl;
        Info<< "    Reynolds number: " << Re_ << endl;
        Info<< "    Geometric angle of attack (degrees): "
            << angleOfAttackGeom_ << endl;
        Info<< "    Angle of attack (uncorrected, degrees): "
            << aoaUncorrected << endl;
        Info<< "    Angle of attack (corrected, degrees): "
            << angleOfAttack_ << endl;
    }

    // Correct coefficients with dynamic stall model
    if (dynamicStallActive_)
    {
        dynamicStall_->correct
        (
            mag(relativeVelocity_),
            angleOfAttack_,
            liftCoefficient_,
            dragCoefficient_,
            momentCoefficient_
        );
    }

    // Correct for added mass effects
    if (addedMassActive_)
    {
        addedMass_.correct
        (
            liftCoefficient_,
            dragCoefficient_,
            momentCoefficient_,
            degToRad(angleOfAttack_),
            mag(chordDirection_ & relativeVelocity_),
            mag(planformNormal_ & relativeVelocity_)
        );
    }
    
    
    // Apply end effect correction factor to lift coefficient
    liftCoefficient_ *= endEffectFactor_;

    // Calculate force per unit density
    scalar area = chordLength_ * spanLength_;
    scalar magSqrU = magSqr(relativeVelocity_);
    scalar lift = 0.5*area*liftCoefficient_*magSqrU;
    scalar drag = 0.5*area*dragCoefficient_*magSqrU;
    vector liftDirection = relativeVelocity_ ^ spanDirection_;
    liftDirection /= mag(liftDirection);
    vector dragDirection = relativeVelocity_/mag(relativeVelocity_);
    forceVector_ = lift*liftDirection + drag*dragDirection;
	
    //Info<< "lift:" << lift << endl;
    //Info<< "drag:" << drag << endl;
    //Info<< "force (per unit density):" << forceVector_ << endl;
    
    normalizedLift = lift / (spanLength_ * 2 * bladeRadius_ * mag(freeStreamVelocity_) * mag(freeStreamVelocity_));
    normalizedDrag = drag / (spanLength_ * 2 * bladeRadius_ * mag(freeStreamVelocity_) * mag(freeStreamVelocity_));	


    if (debug)
    {
        //Info<< "    liftDirection: " << liftDirection << endl;
        //Info<< "    dragDirection: " << dragDirection << endl;
        //Info<< "    force  (per unit density): " << forceVector_ << endl;
    }
}


void Foam::fv::actuatorLineElement::rotate
(
    vector rotationPoint,
    vector axis,
    scalar radians,
    bool rotateVelocity=true
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

    if (debug)
    {
        Info<< "Rotating actuatorLineElement: " << name_ << endl;
        Info<< "Rotation point: " << rotationPoint << endl;
        Info<< "Rotation axis: " << axis << endl;
        Info<< "Rotation angle (radians): " << radians << endl;
        Info<< "Rotation matrix:" << endl << RM << endl;
        Info<< "Initial position: " << position_ << endl;
        Info<< "Initial chordDirection: " << chordDirection_ << endl;
        Info<< "Initial spanDirection: " << spanDirection_ << endl;
        Info<< "Initial velocity: " << velocity_ << endl;
    }

    // Rotation matrices make a rotation about the origin, so need to subtract
    // rotation point off the point to be rotated.
    vector point = position_;
    point -= rotationPoint;

    // Perform the rotation.
    point = RM & point;

    // Return the rotated point to its new location relative to the rotation
    // point
    point += rotationPoint;

    // Set the position of the element
    position_ = point;

    // Rotate the span and chord vectors of the element
    chordDirection_ = RM & chordDirection_;
    spanDirection_ = RM & spanDirection_;
    //Info<< "Chord direction (after pitching RM): " << chordDirection_<< endl;
    // Rotate the element's velocity vector if specified
    if (rotateVelocity)
    {
        velocity_ = RM & velocity_;
        chordRefDirection_ = RM & chordRefDirection_;
    }

    if (debug)
    {
        Info<< "Final position: " << position_ << endl;
        Info<< "Final chordDirection: " << chordDirection_ << endl;
        Info<< "Final chordRefDirection: " << chordRefDirection_ << endl;
        Info<< "Final spanDirection: " << spanDirection_ << endl;
        Info<< "Final velocity: " << velocity_ << endl << endl;
    }
}


void Foam::fv::actuatorLineElement::pitch
(
    scalar radians,
    scalar chordFraction
)
{
    vector rotationPoint = position_;
    rotationPoint += chordDirection_*(chordMount_ - chordFraction);
    rotate(rotationPoint, spanDirection_, radians, false);
}


void Foam::fv::actuatorLineElement::translate(vector translationVector)
{
    position_ += translationVector;
}


void Foam::fv::actuatorLineElement::setVelocity(vector velocity)
{
    if (debug)
    {
        Info<< "Changing velocity of " << name_ << " from "
            << velocity_ << " to " << velocity << endl << endl;
    }
    velocity_ = velocity;
}


void Foam::fv::actuatorLineElement::setSpeed(scalar speed)
{
    if (mag(velocity_) > 0)
    {
        velocity_ /= mag(velocity_);
        velocity_ *= speed;
    }
}


void Foam::fv::actuatorLineElement::setSpeed
(
    vector point,
    vector axis,
    scalar omega
)
{
    if (debug)
    {
        Info<< "Setting speed of " << name_ << " from rotation" << endl;
        Info<< "    Initial velocity: " << velocity_ << endl;
    }

    // First find radius from axis to element position -- formula from
    // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    vector point2 = point + axis;
    scalar radius = mag((position_ - point) ^ (position_ - point2))
                  / mag(point2 - point);
    scalar speed = omega*radius;
    setSpeed(speed);

    scalar angleLE = 0.0;
    scalar angleTE = 0.0;
    if (radius > 0.0)
    {
        // Set velocity at leading edge
        scalar radiusLE = sqrt(magSqr(0.25*chordLength_) + magSqr(radius));
        angleLE = atan2(0.25*chordLength_, radius);
        velocityLE_ = velocity_*radiusLE/radius;
        rotateVector(velocityLE_, vector::zero, spanDirection_, angleLE);

        // Set velocity at trailing edge
        scalar radiusTE = sqrt(magSqr(0.75*chordLength_) + magSqr(radius));
        angleTE = atan2(-0.75*chordLength_, radius);
        velocityTE_ = velocity_*radiusTE/radius;
        rotateVector(velocityTE_, vector::zero, spanDirection_, angleTE);
    }

    // Also set omega for flow curvature correction
    setOmega(omega);

    if (debug)
    {
        Info<< "    Radius: " << radius << endl;
        Info<< "    Final velocity: " << velocity_ << endl;
        Info<< "    Leading edge velocity: " << velocityLE_ << endl;
        Info<< "    Trailing edge velocity: " << velocityTE_ << endl;
        Info<< "    Leading edge velocity angle (radians): "
            << angleLE << endl;
        Info<< "    Trailing edge velocity angle (radians): "
            << angleTE << endl;
    }
}


void Foam::fv::actuatorLineElement::scaleVelocity(scalar scale)
{
    velocity_ *= scale;
}


const Foam::vector& Foam::fv::actuatorLineElement::force()
{
    return forceVector_;
}

Foam::vector Foam::fv::actuatorLineElement::moment(vector point)
{
    // Calculate radius vector
    vector radius = position_ - point;
    vector moment = radius ^ forceVector_;
    vector pitchingMoment = 0.5*chordLength_*chordLength_*spanLength_
                          * momentCoefficient_*magSqr(relativeVelocity_)
                          * spanDirection_;
    return moment + pitchingMoment;
}


void Foam::fv::actuatorLineElement::addSup
(
    fvMatrix<vector>& eqn,
    volVectorField& forceField
    //volScalarField& factorField
)
{
    volVectorField forceFieldI
    (
        IOobject
        (
            "force." + name_,
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            forceField.dimensions(),
            vector::zero
        )
    );

	/*volScalarField factorFieldI
	(
		IOobject
		(   "gBlade." + name_,
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_
    );
       */ 
    const volVectorField& Uin(eqn.psi());
    calculateForce(Uin);
    applyForceField(forceFieldI);

// Pout<< "Processor" << Pstream::myProcNo() << "forceField:" << forceField << endl;


    // Add force to total actuator line force
    forceField += forceFieldI;
 
	//factorField += factorFieldI;
	
    // Write performance to file
    if (writePerf_ and Pstream::master())
    {
        writePerf();
    }
}


void Foam::fv::actuatorLineElement::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    volVectorField& forceField
)
{
    volVectorField forceFieldI
    (
        IOobject
        (
            "force." + name_,
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            forceField.dimensions()/rho.dimensions(),
            vector::zero
        )
    );

    const volVectorField& Uin(eqn.psi());
    calculateForce(Uin);
    applyForceField(forceFieldI);

    // Multiply force vector by local density
    multiplyForceRho(rho);

    // Multiply this element's force field by density field
    forceFieldI *= rho;

    // Add force to total actuator line force
    forceField += forceFieldI;

    // Write performance to file
    if (writePerf_ and Pstream::master())
    {
        writePerf();
    }
}


void Foam::fv::actuatorLineElement::addTurbulence
(
    fvMatrix<scalar>& eqn,
    word fieldName
)
{
    volScalarField turbulence
    (
        IOobject
        (
            "turbulence." + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "zero",
            eqn.dimensions()/dimVolume,
            0.0
        )
    );

    // Calculate projection radius
    vector epsilon = calcProjectionEpsilon();
    scalar projectionRadius = (epsilon[0]*Foam::sqrt(Foam::log(1.0/0.001)));

    // Calculate TKE injection rate
    scalar k = 0.1*mag(dragCoefficient_);

    // Add turbulence to the cells within the element's sphere of influence
    scalar sphereRadius = chordLength_ + projectionRadius;
    forAll(mesh_.cells(),cellI)
    {
	scalar dis = mag(mesh_.C()[cellI] - position_);
	if (dis <= sphereRadius)
	{
	    scalar factor = Foam::exp(-Foam::sqr(dis/epsilon[0]))
			   / (Foam::pow(epsilon[0],3)
			   * Foam::pow(Foam::constant::mathematical::pi, 1.5));
	    if (fieldName == "k")
	    {
		turbulence[cellI] = factor*k;
	    }
	    else if (fieldName == "epsilon")
	    {
		turbulence[cellI] = factor* Foam::pow(k, 1.5)
				   * 0.09/(chordLength_/10.0);
	    }
	}
   }

   eqn += turbulence;
}

void Foam::fv::actuatorLineElement::setDynamicStallActive(bool active)
{
   dynamicStallActive_ = active;
}

void Foam::fv::actuatorLineElement::setFilteredLiftingLineActive(bool active)
{
   filteredLiftingLineActive_ = active;
}

void Foam::fv::actuatorLineElement::setFilteredLiftingLineWriteOutput(bool active)
{
   filteredLiftingLineWriteOutput_ = active;
}

void Foam::fv::actuatorLineElement::setOmega(scalar omega)
{
   omega_ = omega;
}

void Foam::fv::actuatorLineElement::setEndEffectFactor(scalar factor)
{
   endEffectFactor_ = factor;
}

void Foam::fv::actuatorLineElement::setVelocitySampleRadius(scalar radius)
{
   velocitySampleRadius_ = radius;
}

void Foam::fv::actuatorLineElement::setNVelocitySamples(label nSamples)
{
    nVelocitySamples_ = nSamples;
}


// ************************************************************************* //
