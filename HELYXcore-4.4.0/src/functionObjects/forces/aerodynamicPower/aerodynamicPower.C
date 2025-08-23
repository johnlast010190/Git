/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of HELYXcore.
    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    HELYXcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2015-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "aerodynamicPower/aerodynamicPower.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "cfdTools/general/porosityModel/porosityModel/porosityModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "coordinate/systems/cartesianCS.H"
#include "fields/GeometricFields/GeometricField/GeometricIOBoundaryField.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"
#include "pointPatchFields/derived/modalForcedMotion/modalForcedMotionPointPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(aerodynamicPower, 0);
    addToRunTimeSelectionTable(functionObject, aerodynamicPower, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::tmp<Foam::volScalarField> Foam::functionObjects::aerodynamicPower::getp()
{

    const volScalarField& p = lookupObject<volScalarField>(pName_);

    if (p.dimensions() == dimensionSet(1, -1, -2, 0, 0, 0, 0))
    {
        return p - dimensionedScalar("Pref", dimPressure, pRef_);
    }
    else
    {
        dimensionedScalar dimRho
        (
            "rho",
            dimensionSet(1, -3, 0, 0, 0, 0, 0),
            rhoRef_
        );

        return dimRho*p - dimensionedScalar("Pref", dimPressure, pRef_);
    }
}


Foam::word Foam::functionObjects::aerodynamicPower::fieldName(const word& name) const
{
    return this->name() + "_" + name;
}


void Foam::functionObjects::aerodynamicPower::createFiles()
{
    // Note: Only possible to create bin files after bins have been initialised

    if (writeToFile() && !powerFilePtr_.valid())
    {

        powerFilePtr_ = createFile("aerodynamicPower");
        writeIntegratedHeader();


        if (nBin_ > 1)
        {
            powerBinFilePtr_ = createFile("aerodynamicPower_bins");
            powerAvgBinFilePtr_ = createFile("aerodynamicPowerAvg_bins");
            velBinFilePtr_ = createFile("solidVelocity_bins");
            writeBinHeader();
        }
    }
}


void Foam::functionObjects::aerodynamicPower::writeIntegratedHeader()
{
    Ostream& os = powerFilePtr_();
    writeHeaderValue(os, "CofR", csys().origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeDelimited(os, "TimeAveragedAeroPower");
    writeDelimited(os, "aeroPower");
    writeDelimited(os, "aeroPowerPatch");
    writeDelimited(os, "totalBinAeroPower");

    os  << endl;
}



void Foam::functionObjects::aerodynamicPower::writeBinHeader()
{
    Ostream& os = powerBinFilePtr_();

    os << "# bins      : " << nBin_ << nl
       << "# start     : " << binMin_ << nl
       << "# delta     : " << binDx_ << nl
       << "# direction : " << binDir() << nl
       << "# Time";

    for (label j = 0; j < nBin_; j++)
    {
        const word jn('[' + Foam::name(j) + ']');

        os << tab
           << "aeroPower" << jn;
    }

    os << endl;


    Ostream& os2 = powerAvgBinFilePtr_();

    os2 << "# bins      : " << nBin_ << nl
       << "# start     : " << binMin_ << nl
       << "# delta     : " << binDx_ << nl
       << "# direction : " << binDir() << nl
       << "# Time";

    for (label j = 0; j < nBin_; j++)
    {
        const word jn('[' + Foam::name(j) + ']');

        os2 << tab
           << "timeAveragedAeroPower" << jn;
    }

    os2 << endl;


    Ostream& os3 = velBinFilePtr_();

    os3 << "# bins      : " << nBin_ << nl
       << "# start     : " << binMin_ << nl
       << "# delta     : " << binDx_ << nl
       << "# direction : " << binDir() << nl
       << "# Time";

    for (label j = 0; j < nBin_; j++)
    {
        const word jn('[' + Foam::name(j) + ']');

        os3 << tab
           << "solidVelocity" << jn;
    }

    os3 << endl;

}


void Foam::functionObjects::aerodynamicPower::initialise()
{
    if (initialised_)
    {
        return;
    }

    if (directForceDensity_)
    {
        if (!foundObject<volVectorField>(fDName_))
        {
            FatalErrorInFunction
                << "Could not find " << fDName_ << " in database"
                << exit(FatalError);
        }
    }
    else
    {
        if
        (
            !foundObject<volVectorField>(UName_)
         || !foundObject<volScalarField>(pName_)

        )
        {
            FatalErrorInFunction
                << "Could not find U: " << UName_ << " or p: " << pName_
                << " in database"
                << exit(FatalError);
        }

        if (rhoName_ != "rhoInf" && !foundObject<volScalarField>(rhoName_))
        {
            FatalErrorInFunction
                << "Could not find rho: " << rhoName_
                << exit(FatalError);
        }
    }

    initialiseBins();
    totalIter_ = 1;
    totalTime_ = obr().time().deltaTValue();
    initialised_ = true;
}


void Foam::functionObjects::aerodynamicPower::initialiseBins()
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // Determine extents of patches
    binMin_ = GREAT;
    scalar binMax = -GREAT;
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const polyPatch& pp = pbm[patchi];
        scalarField d(pp.faceCentres() & binDir());
        binMin_ = min(min(d), binMin_);
        binMax = max(max(d), binMax);
    }

    // Include porosity
    if (porosity_)
    {
        const HashTable<const porosityModel*> models =
            obr_.lookupClass<porosityModel>();

        const scalarField dd(mesh_.C() & binDir());

        forAllConstIter(HashTable<const porosityModel*>, models, iter)
        {
            const porosityModel& pm = *iter();
            const labelList& cellZoneIDs = pm.cellZoneIDs();

            forAll(cellZoneIDs, i)
            {
                label zonei = cellZoneIDs[i];
                const cellZone& cZone = mesh_.cellZones()[zonei];

                bool calculateForce = checkZone(cZone);

                if (calculateForce)
                {
                    const scalarField d(dd, cZone);
                    binMin_ = min(min(d), binMin_);
                    binMax = max(max(d), binMax);
                }
            }
        }
    }

    reduce(binMin_, minOp<scalar>());
    reduce(binMax, maxOp<scalar>());

    // Slightly boost binMax so that region of interest is fully
    // within bounds
    binMax = 1.0001*(binMax - binMin_) + binMin_;

    binDx_ = (binMax - binMin_)/scalar(nBin_);

    // Create the bin points used for writing
    binPoints_.setSize(nBin_);
    forAll(binPoints_, i)
    {
        binPoints_[i] = binMin_*binDir() + (i + 0.5)*binDir()*binDx_;
    }

    // Allocate storage for aerodynamicPower and moments
    forAll(force_, i)
    {
        force_[i].setSize(nBin_, vector::zero);
    }
    aeroPowerAvg_.setSize(nBin_, 0);
    aeroPower_.setSize(nBin_, 0);
}


void Foam::functionObjects::aerodynamicPower::resetFields()
{
    force_[0] = Zero;
    force_[1] = Zero;
    force_[2] = Zero;

    aeroPowerAvg_ = Zero;
    aeroPower_ = Zero;
    timePowerPatch_ = Zero; // reset as it sums up over each patches entry

    if (writeFields_)
    {
        volVectorField& force =
            const_cast<volVectorField&>
            (
                lookupObject<volVectorField>(fieldName("force"))
            );

        force.forceAssign(dimensionedVector(force.dimensions(), Zero));
    }
}


bool Foam::functionObjects::aerodynamicPower::checkZone(const cellZone& cZone)
{
    bool calculateForce = true;
    if (porosityZones_.size())
    {
        const word& cZoneName = cZone.name();
        calculateForce = false;
        forAll(porosityZones_, zI)
        {
            if (porosityZones_[zI] == cZoneName)
            {
                calculateForce = true;
            }
        }
    }
    return calculateForce;
}


Foam::tmp<Foam::volSymmTensorField::Boundary>
Foam::functionObjects::aerodynamicPower::bDevRhoReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    tmp<volSymmTensorField> tDevRhoReff;

    if (foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        // Only use caching for this option as the others can be affected
        // by different values of rhoInf in different function objects
        return cachedBDevRhoReff();
    }
    else if (foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const incompressible::turbulenceModel& turb =
            lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        tDevRhoReff = rho()*turb.devReff();
    }
    else if (foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const fluidThermo& thermo =
            lookupObject<fluidThermo>(fluidThermo::dictName);

        const volVectorField& U = lookupObject<volVectorField>(UName_);

        tDevRhoReff = -thermo.mu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const transportModel& laminarT =
            lookupObject<transportModel>("transportProperties");

        const volVectorField& U = lookupObject<volVectorField>(UName_);

        tDevRhoReff = -rho()*laminarT.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
            lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu("nu", transportProperties.lookup("nu"));

        const volVectorField& U = lookupObject<volVectorField>(UName_);

        tDevRhoReff = -rho()*nu*dev(twoSymm(fvc::grad(U)));
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);
        ::abort();
    }
    return
        tmp<volSymmTensorField::Boundary>
        (
            new volSymmTensorField::Boundary
            (
                volSymmTensorField::null(),
                tDevRhoReff.ref().boundaryFieldRef(),
                true
            )
        );
}


Foam::tmp<Foam::volSymmTensorField::Boundary>
Foam::functionObjects::aerodynamicPower::cachedBDevRhoReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;

    const cmpTurbModel& turb =
        lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

    if (devRhoReffCachingActive_)
    {
        typedef GeometricIOBoundaryField<symmTensor, fvPatchField, volMesh>
            GeoBFType;
        GeoBFType* cachedFieldPtr =
            obr_.lookupObjectRefPtr<GeoBFType>(bDevRhoReffName());
        if (!cachedFieldPtr)
        {
            cachedFieldPtr =
                new GeoBFType
                (
                    IOobject
                    (
                        bDevRhoReffName(),
                        mesh_.time().timeName(),
                        obr_
                    ),
                    volSymmTensorField::null(), // dummy internal field ref
                    turb.devRhoReff().ref().boundaryFieldRef(),
                    true
                );
            cachedFieldPtr->rename(bDevRhoReffName());
            cachedFieldPtr->store();
            if (debug)
            {
                Info<< "aerodynamicPower: Calculated and stored devRhoReff field"
                    << endl;
            }
        }
        else
        {
            if (debug)
            {
                Info<< "aerodynamicPower: Used cached devRhoReff field" << endl;
            }
        }
        return tmp<volSymmTensorField::Boundary>(*cachedFieldPtr);
    }
    else
    {
        if (debug)
        {
            Info<< "aerodynamicPower: Calculating devRhoReff field" << endl;
        }
        return
            tmp<volSymmTensorField::Boundary>
            (
                new volSymmTensorField::Boundary
                (
                    volSymmTensorField::null(),
                    turb.devRhoReff().ref().boundaryFieldRef(),
                    true
                )
            );
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::aerodynamicPower::mu() const
{
    if (foundObject<fluidThermo>(basicThermo::dictName))
    {
        const fluidThermo& thermo =
            lookupObject<fluidThermo>(basicThermo::dictName);

        return thermo.mu();
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const transportModel& laminarT =
            lookupObject<transportModel>("transportProperties");

        return rho()*laminarT.nu();
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
            lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu
        (
            "nu",
            dimViscosity,
            transportProperties.lookup("nu")
        );

        return rho()*nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::aerodynamicPower::rho() const
{
    if (rhoName_ == "rhoInf")
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("rho", dimDensity, rhoRef_)
            )
        );
    }
    else
    {
        return(lookupObject<volScalarField>(rhoName_));
    }
}


Foam::scalar Foam::functionObjects::aerodynamicPower::rho(const volScalarField& p) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1.0;
    }
    else
    {
        if (rhoName_ != "rhoInf")
        {
            FatalErrorInFunction
                << "Dynamic pressure is expected but kinematic is provided."
                << exit(FatalError);
        }

        return rhoRef_;
    }
}


void Foam::functionObjects::aerodynamicPower::applyBins
(
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP,
    const vectorField& d
)
{
    if (nBin_ == 1)
    {
        force_[0][0] += sum(fN);
        force_[1][0] += sum(fT);
        force_[2][0] += sum(fP);
    }
    else
    {
        scalarField dd((d & binDir()) - binMin_);

        forAll(dd, i)
        {
            label bini = min
            (
                max(floor(dd[i]/binDx_), 0),
                force_[0].size() - 1
            );

            force_[0][bini] += fN[i];
            force_[1][bini] += fT[i];
            force_[2][bini] += fP[i];
        }
    }

    //- get solid body velocity
    autoPtr<vectorField> solidVelocityPtr = solidVelocityBin();
    const vectorField& Usolid = solidVelocityPtr();

    forAll(force_[0], bi)
    {
        aeroPower_[bi] = (force_[0][bi] +force_[1][bi] +force_[2][bi])
            & Usolid[bi];

        if (debug)
        {
            Info<< "Aero power: " << aeroPower_[bi]
                << ", Usolid " << Usolid[bi]
                << ", force (" << force_[0][bi]
                << ", "<< force_[1][bi]
                << ", "<< force_[2][bi]
                << ")" << endl;
        }
    }

    //- time integration
    timePower_ = calculateBinsIntegral(aeroPower_);

    aeroPowerAvg_ = calculateMeanField(aeroPower_)();

    //- spatial integration of the different bins
    power_ = calculateBinsIntegral(aeroPowerAvg_);
}


void Foam::functionObjects::aerodynamicPower::integrateOnPatch
(
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP,
    label patchID
)
{
    autoPtr<vectorField> solidVelocityPtr = solidVelocity(patchID);
    const vectorField& Usolid = solidVelocityPtr();

    const scalarField& magSf(mesh_.magSf().boundaryField()[patchID]);

    if (areaWeighted_)
    {
        timePowerPatch_ += gSum(((fN +fT +fP) & Usolid) * magSf);
    }
    else
    {
        timePowerPatch_ += gSum((fN +fT +fP) & Usolid);
    }
}


Foam::autoPtr<Foam::vectorField> Foam::functionObjects::aerodynamicPower::solidVelocityBin()
{
    autoPtr<vectorField> Usolid( new vectorField(nBin_, Zero));

    scalar timeNew = time_.value();

    Foam::pointVectorField& pointDisp =
        mesh_.lookupObjectRef<pointVectorField>("pointDisplacement");

    forAll(pointDisp.boundaryField(), pI)
    {
        if
        (
            isA<modalForcedMotionPointPatchVectorField>
            (pointDisp.boundaryField()[pI])
        )
        {
            const modalForcedMotionPointPatchVectorField& ptf =
                refCast<const modalForcedMotionPointPatchVectorField>
                (pointDisp.boundaryField()[pI]);

            forAll(binPoints_, i)
            {
                scalar radius = binPoints_[i].z();
                Usolid()[i] = ptf.radialDisplacement().value(radius)*
                    ptf.amplification().derivative(timeNew);
            }
        }
    }

    return Usolid;
}


Foam::autoPtr<Foam::vectorField> Foam::functionObjects::aerodynamicPower::solidVelocity(label patchID)
{
    autoPtr<vectorField> Usolid
    (
        new vectorField(mesh_.C().boundaryField()[patchID].size(), Zero)
    );

    scalar timeNew = time_.value();

    Foam::pointVectorField& pointDisp =
        mesh_.lookupObjectRef<pointVectorField>("pointDisplacement");

    if
    (
        isA<modalForcedMotionPointPatchVectorField>
        (pointDisp.boundaryField()[patchID])
    )
    {
        const modalForcedMotionPointPatchVectorField& ptf =
            refCast<const modalForcedMotionPointPatchVectorField>
            (pointDisp.boundaryField()[patchID]);

        const vectorField& Cf = mesh_.boundary()[patchID].Cf();
        forAll(Cf, i)
        {
            scalar radius = Cf[i].z();
            Usolid()[i] = ptf.radialDisplacement().value(radius)*
                ptf.amplification().derivative(timeNew);
        }
    }

    return Usolid;
}


Foam::autoPtr<Foam::scalarField> Foam::functionObjects::aerodynamicPower::calculateMeanField
(
    const scalarField& baseField
)
{
    autoPtr<scalarField>  meanField( new scalarField(baseField.size(), 0));

    scalar dt = mesh_.time().deltaTValue();

    if (totalIter_ == 1)
    {
        meanField() = baseField*dt;
    }
    else
    {
        meanField() = meanField() + (baseField*dt)/window_;
    }

    totalTime_ += mesh_.time().deltaTValue();
    totalIter_++;

    return meanField;
}


Foam::scalar Foam::functionObjects::aerodynamicPower::calculateBinsIntegral
(
    const scalarField& baseField
)
{
    scalar meanValue = 0;
    if (nBin_ == 1)
    {
        meanValue = baseField[0];
    }
    else
    {
        // loop all components, i.e. all bins in this case
        for (label bini=1; bini<baseField.size(); bini++)
        {
            //- trapezoidal rule
            meanValue += 0.5*(baseField[bini-1] +baseField[bini])*
                mag(binPoints_[bini] -binPoints_[bini-1]); //dx
        }
    }
    reduce(meanValue, sumOp<scalar>());
    return meanValue;
}

void Foam::functionObjects::aerodynamicPower::addToFields
(
    const label patchi,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP
)
{
    if (!writeFields_)
    {
        return;
    }

    volVectorField& force =
        lookupObjectRef<volVectorField>(fieldName("force"));

    vectorField& pf = force.boundaryFieldRef()[patchi];
    pf += fN + fT + fP;
}


void Foam::functionObjects::aerodynamicPower::addToFields
(
    const labelList& cellIDs,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP
)
{
    if (!writeFields_)
    {
        return;
    }

    volVectorField& force =
        lookupObjectRef<volVectorField>(fieldName("force"));

    forAll(cellIDs, i)
    {
        label celli = cellIDs[i];
        force[celli] += fN[i] + fT[i] + fP[i];
    }
}


void Foam::functionObjects::aerodynamicPower::writeIntegratedAeroPower()
{
    scalar integratedPower = sum(aeroPower_);

    Log << " Time Avg Aero Power =  "
        << power_ << nl
        << "  Aero Power (bins)  =  "
        << timePower_ << nl
        << " Aero Power (patch)  =  "
        << timePowerPatch_ << nl
        << "Total Bin Aero Power =  "
        << integratedPower  << nl
        << nl;

    if (writeToFile())
    {
        Ostream& os = powerFilePtr_();
        writeTime(os);
        os << tab << power_ << ", "
           << tab << timePower_ << ", "
           << tab << timePowerPatch_ << ", "
           << integratedPower;
        os << endl;
    }
}


void Foam::functionObjects::aerodynamicPower::writeAeroPower()
{
    Log << type() << " " << name() << " write:" << nl;

    writeIntegratedAeroPower();

    Log << endl;
}


void Foam::functionObjects::aerodynamicPower::writeBinnedAeroPower()
{
    if ((nBin_ == 1) || !writeToFile())
    {
        return;
    }

    //- write aero power
    scalarField& aeroPower(aeroPower_);

    if (binCumulative_)
    {
        for (label i=1; i<aeroPower.size(); i++)
        {
            aeroPower[i] += aeroPower[i-1];
        }
    }

    Ostream& os = powerBinFilePtr_();

    writeTime(os);

    forAll(aeroPower, i)
    {
        os  << tab
            << aeroPower[i] ;
    }
    os  << nl;

    //- write time averaged power
    scalarField& aeroPowerAvg(aeroPowerAvg_);

    if (binCumulative_)
    {
        for (label i=1; i<aeroPowerAvg.size(); i++)
        {
            aeroPowerAvg[i] += aeroPowerAvg[i-1];
        }
    }

    Ostream& os2 = powerAvgBinFilePtr_();

    writeTime(os2);

    forAll(aeroPowerAvg, i)
    {
        os2  << tab
            << aeroPowerAvg[i] ;
    }
    os2  << nl;

    //- write solid velocity
    vectorField vel(solidVelocityBin()());
    Ostream& os3 = velBinFilePtr_();

    writeTime(os3);

    forAll(vel, i)
    {
        os3  << tab
             << vel[i] ;
    }

    os3  << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::aerodynamicPower::aerodynamicPower
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool baseClass
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name),
    force_(3),
    powerFilePtr_(),
    powerBinFilePtr_(),
    velBinFilePtr_(),
    powerAvgBinFilePtr_(),
    patchSet_(),
    pName_(word::null),
    UName_(word::null),
    rhoName_(word::null),
    directForceDensity_(false),
    fDName_(""),
    devRhoReffCaching_(true),
    devRhoReffCachingActive_(false),
    rhoRef_(VGREAT),
    pRef_(0),
    coordSys_(),
    coorFramePtr_(nullptr),
    definedInFrame_(false),
    localSystem_(false),
    porosity_(false),
    porosityZones_(wordList::null()),
    nBin_(1),
    binDir_(Zero),
    binDx_(0.0),
    binMin_(GREAT),
    binPoints_(),
    binCumulative_(true),
    aeroPower_(nBin_),
    aeroPowerAvg_(nBin_),
    power_(0),
    totalTime_(),
    totalIter_(),
    window_(-1),
    writeFields_(false),
    initialised_(false),
    combineFiles_(false),
    areaWeighted_(false)
{
    if (!baseClass)
    {
        aerodynamicPower::read(dict);
        Log << endl;
    }
}


Foam::functionObjects::aerodynamicPower::aerodynamicPower
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool baseClass
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(obr_, name),
    force_(3),
    powerFilePtr_(),
    powerBinFilePtr_(),
    velBinFilePtr_(),
    powerAvgBinFilePtr_(),
    patchSet_(),
    pName_(word::null),
    UName_(word::null),
    rhoName_(word::null),
    directForceDensity_(false),
    fDName_(""),
    devRhoReffCaching_(true),
    devRhoReffCachingActive_(false),
    rhoRef_(VGREAT),
    pRef_(0),
    coordSys_(),
    coorFramePtr_(nullptr),
    definedInFrame_(false),
    localSystem_(false),
    porosity_(false),
    porosityZones_(wordList::null()),
    nBin_(1),
    binDir_(Zero),
    binDx_(0.0),
    binMin_(GREAT),
    binPoints_(),
    binCumulative_(true),
    aeroPower_(nBin_),
    aeroPowerAvg_(nBin_),
    power_(0),
    totalTime_(),
    totalIter_(),
    window_(-1),
    writeFields_(false),
    initialised_(false),
    combineFiles_(false),
    areaWeighted_(false)
{
    if (!baseClass)
    {
        aerodynamicPower::read(dict);
        Log << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::aerodynamicPower::~aerodynamicPower()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::coordinateSystem& Foam::functionObjects::aerodynamicPower::csys() const
{
    if (coorFramePtr_)
    {
        return coorFramePtr_->coorSys();
    }
    return coordSys_;
}


Foam::vector Foam::functionObjects::aerodynamicPower::binDir() const
{
    if (coorFramePtr_ && definedInFrame_)
    {
        return coorFramePtr_->coorSys().globalVector(binDir_);
    }
    return binDir_;
}


bool Foam::functionObjects::aerodynamicPower::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    initialised_ = false;

    Info<< type() << " " << name() << ":" << nl;

    directForceDensity_ = dict.lookupOrDefault("directForceDensity", false);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    //set patches used by aerodynamicPower
    patchSet_.clear();

    if (dict.found("patches"))
    {
        patchSet_ =
            pbm.patchSet(wordReList(dict.lookup("patches")), false, true);
    }
    //patchSet_ = pbm.patchSet(wordReList(dict.lookup("patches")));

    devRhoReffCaching_ =
        dict.lookupOrDefault<Switch>("devRhoReffCaching", true);

    if (directForceDensity_)
    {
        // Optional entry for fDName
        fDName_ = dict.lookupOrDefault<word>("fD", "fD");
    }
    else
    {
        // Optional entries U and p
        const wordList pNames({"p", "pName"});
        pName_ = dict.lookupOrDefault<word>(pNames, "p");

        const wordList UNames({"U", "UName"});
        UName_ = dict.lookupOrDefault<word>(UNames, "U");

        const wordList rhoNames({"rho", "rhoName"});
        rhoName_ = dict.lookupOrDefault<word>(rhoNames, "rho");

        if (dict.found("rhoInf") && !dict.found("rho") && !dict.found("rhoName"))
        {
            rhoName_ = "rhoInf";
        }

        if (dict.found("rhoName") && dict.found("rho"))
        {
            WarningInFunction
                << "\n    Since both entries 'rho' and 'rhoName' "
                << " were found in aerodynamicPower function object, '"
                << rhoNames[0] << "' value from 'rho' entry "
                << "is being used for aerodynamicPower calculation\n" << endl;
        }

        // Reference density needed for incompressible calculations
        if (rhoName_ == "rhoInf")
        {
            rhoRef_ = dict.lookup<scalar>("rhoInf");
        }

        // Reference pressure, 0 by default
        const wordList pRefNames({"pRef", "Pref"});
        pRef_ = dict.lookupOrDefault<scalar>(pRefNames, 0.0);
    }

    if (dict.found("referenceFrame"))
    {
        coorFramePtr_ = coordinateFrame::lookupNew(mesh_, dict);
        definedInFrame_ =
            dict.lookupOrDefault<Switch>("definedInFrame", false);
        localSystem_ = true;
        if (dict.found("CofR") || dict.found("referencePoint"))
        {
            WarningInFunction
                << "Reference frame is switched on. "
                << "\"CofR\" and \"referencePoint\" have no effect."
                << endl;
        }
    }
    else
    {
        // Centre of rotation for moment calculations
        // specified directly, from coordinate system, or implicitly (0 0 0)
        coordSys_.clear();
        if (!dict.readIfPresent<point>("CofR", coordSys_.origin()))
        {
            if (!dict.readIfPresent<point>("referencePoint", coordSys_.origin()))
            {
                // The 'coordinateSystem' sub-dictionary is optional,
                // but enforce use of a cartesian system.
                if (dict.found(coordinateSystem::typeName_()))
                {
                    // New() for access to indirect (global) coordinate system
                    coordSys_ =
                        coordinateSystem::New
                        (
                            obr_, dict, coordinateSystem::typeName_()
                        );
                }
                else
                {
                    coordSys_ = coordSystem::cartesian(dict);
                }

                localSystem_ = true;
            }
        }
    }

    dict.readIfPresent("porosity", porosity_);
    if (porosity_)
    {
        Info<< "    Including porosity effects" << endl;
        porosityZones_ = dict.lookupOrDefault<wordList>
        (
            "porosityZones", wordList::null()
        );
    }
    else
    {
        Info<< "    Not including porosity effects" << endl;
    }

    if (dict.found("binData"))
    {
        const dictionary& binDict(dict.subDict("binData"));

        // check for new format of liftDrag (uses "nBins")
        if (binDict.found("nBin"))
        {
            nBin_ = binDict.lookup<label>("nBin");

            if (nBin_ < 0)
            {
                FatalIOErrorInFunction(dict)
                    << "Number of bins (nBin) must be zero or greater"
                    << exit(FatalIOError);
            }
            else if (nBin_ == 0)
            {
                // Case of no bins equates to a single bin to collect all data
                nBin_ = 1;
            }
            else
            {
                binCumulative_ = binDict.lookup<bool>("cumulative");
                binDir_ = binDict.lookup<vector>("direction");
                binDir_ /= mag(binDir_);
            }
        }
    }

    writeFields_ = dict.lookupOrDefault("writeFields", false);

    if (writeFields_)
    {
        Info<< "    Fields will be written" << endl;

        if (!mesh_.foundObject<volVectorField>(fieldName("force")))
        {
            volVectorField* forcePtr
            (
                new volVectorField
                (
                    IOobject
                    (
                        fieldName("force"),
                        time_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedVector(dimForce, Zero)
                )
            );

            mesh_.objectRegistry::store(forcePtr);
        }
    }

    dict.readIfPresent("combineFiles", combineFiles_);
    dict.readIfPresent("areaWeighted", areaWeighted_);

    window_ = dict.lookupOrDefault<scalar>("window", -1.0);

    return true;
}


void Foam::functionObjects::aerodynamicPower::calcForcesMoment()
{
    if (coorFramePtr_ && definedInFrame_)
    {
        initialiseBins();
    }

    initialise();

    resetFields();
    if (directForceDensity_)
    {
        const volVectorField& fD = lookupObject<volVectorField>(fDName_);

        const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();

            vectorField Md
            (
                mesh_.C().boundaryField()[patchi] - csys().origin()
            );

            scalarField sA(mag(Sfb[patchi]));

            // Normal force = surfaceUnitNormal*(surfaceNormal & forceDensity)
            vectorField fN
            (
                Sfb[patchi]/sA
               *(
                    Sfb[patchi] & fD.boundaryField()[patchi]
                )
            );

            // Tangential force (total force minus normal fN)
            vectorField fT(sA*fD.boundaryField()[patchi] - fN);

            // Porous force
            vectorField fP(Md.size(), Zero);

            addToFields(patchi, fN, fT, fP);

            applyBins(fN, fT, fP, mesh_.C().boundaryField()[patchi]);

            integrateOnPatch(fN, fT, fP, patchi);
        }
    }
    else
    {
        tmp<volScalarField> p = getp();

        const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

        tmp<volSymmTensorField::Boundary> tdevRhoReffb = bDevRhoReff();

        const volSymmTensorField::Boundary& devRhoReffb = tdevRhoReffb();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();

            vectorField Md
            (
                mesh_.C().boundaryField()[patchi] - csys().origin()
            );

            vectorField fN(Sfb[patchi]*p->boundaryField()[patchi]);

            vectorField fT(Sfb[patchi] & devRhoReffb[patchi]);

            vectorField fP(Md.size(), Zero);

            addToFields(patchi, fN, fT, fP);

            applyBins(fN, fT, fP, mesh_.C().boundaryField()[patchi]);

            integrateOnPatch(fN, fT, fP, patchi);
        }
    }
    if (porosity_)
    {
        const volVectorField& U = lookupObject<volVectorField>(UName_);
        const volScalarField rho(this->rho());
        const volScalarField mu(this->mu());

        const HashTable<const porosityModel*> models =
            obr_.lookupClass<porosityModel>();

        if (models.empty())
        {
            WarningInFunction
                << "Porosity effects requested, but no porosity models found "
                << "in the database"
                << endl;
        }

        forAllConstIter(HashTable<const porosityModel*>, models, iter)
        {
            // Non-const access required if mesh is changing
            porosityModel& pm = const_cast<porosityModel&>(*iter());

            vectorField fPTot(pm.force(U, rho, mu));

            const labelList& cellZoneIDs = pm.cellZoneIDs();

            forAll(cellZoneIDs, i)
            {
                label zonei = cellZoneIDs[i];
                const cellZone& cZone = mesh_.cellZones()[zonei];

                const vectorField d(mesh_.C(), cZone);
                const vectorField fP(fPTot, cZone);
                const vectorField Md(d - csys().origin());

                const vectorField fDummy(Md.size(), Zero);

                bool calculateForce = checkZone(cZone);

                if (calculateForce)
                {
                    addToFields(cZone, fDummy, fDummy, fP);
                    applyBins(fDummy, fDummy, fP, d);
                }
            }
        }
    }

    Pstream::listCombineGather(force_, plusEqOp<vectorField>());
    Pstream::listCombineScatter(force_);
    Pstream::listCombineGather(aeroPower_, plusEqOp<scalar>());
    Pstream::listCombineScatter(aeroPower_);
    Pstream::listCombineGather(aeroPowerAvg_, plusEqOp<scalar>());
    Pstream::listCombineScatter(aeroPowerAvg_);
}


Foam::vector Foam::functionObjects::aerodynamicPower::forceEff() const
{
    return sum(force_[0]) + sum(force_[1]) + sum(force_[2]);
}


bool Foam::functionObjects::aerodynamicPower::execute()
{
    calcForcesMoment();

    if (Pstream::master())
    {
        createFiles();

        writeAeroPower();

        writeBinnedAeroPower();

        Log << endl;
    }

    setResult("aerodynamicPower", sum(aeroPowerAvg_));

    return true;
}


bool Foam::functionObjects::aerodynamicPower::write()
{
    return true;
}


void Foam::functionObjects::aerodynamicPower::timeIncremented()
{
    // If this function is called, we know this is actually being run as a
    // function object and we can safely activate the caching as it will be
    // cleared
    devRhoReffCachingActive_ = devRhoReffCaching_;

    if (devRhoReffCachingActive_)
    {
        typedef GeometricIOBoundaryField<symmTensor, fvPatchField, volMesh>
            GeoBFType;
        GeoBFType* cachedFieldPtr =
            obr_.lookupObjectRefPtr<GeoBFType>(bDevRhoReffName());
        if (cachedFieldPtr)
        {
            cachedFieldPtr->release();
            delete cachedFieldPtr;
            if (debug)
            {
                Info<< "aerodynamicPower: Cleared cached devRhoReff field" << endl;
            }
        }
    }
}


// ************************************************************************* //
