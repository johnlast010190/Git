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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2010-2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.H"
#include "scatterModels/scatterModel/scatterModel.H"
#include "global/constants/constants.H"
#include "finiteVolume/fvm/fvm.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "boundaryRadiationProperties/boundaryRadiationProperties.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiationModels
    {
        defineTypeNameAndDebug(fvDOM, 0);
        addToRadiationRunTimeSelectionTables(fvDOM);
    }

    template<>
    const char* NamedEnum<radiationModels::fvDOM::solarLoadType, 3>::names[] =
    {
        "viewFactor",
        "DOM",
        "none"
    };

    const NamedEnum<radiationModels::fvDOM::solarLoadType, 3>
        radiationModels::fvDOM::solarLoadTypeNames_(0);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiationModels::fvDOM::setMaxIter()
{
    if (!coeffs_.found("initialMaxIter"))
    {
        if (participating())
        {
            label nCells = mesh_.nCells();
            reduce(nCells, sumOp<label>());
            if (mesh_.nSolutionD() == 3)
            {
                initMaxIter_ = int(pow(nCells, 1.0/3.0));
            }
            else if (mesh_.nSolutionD() == 2)
            {
                initMaxIter_ = int(Foam::sqrt(scalar(nCells)));
            }
            else
            {
                initMaxIter_ = nCells;
            }
        }
        else
        {
            initMaxIter_ = max(initMaxIter_, Pstream::nProcs());
        }
    }
}

void Foam::radiationModels::fvDOM::rotateInitialRays(const vector& sunDir)
{
    // Rotate Y spherical cordinates to Sun direction.
    // Solid angles on the equator are better fit for planar radiation
    const tensor coordRot = rotationTensor(vector(0, 1, 0), sunDir);

    forAll(IRay_, rayId)
    {
        IRay_[rayId].dAve() = coordRot & IRay_[rayId].dAve();
        IRay_[rayId].d() = coordRot & IRay_[rayId].d();
    }
}


void Foam::radiationModels::fvDOM::alignClosestRayToSun(const vector& sunDir)
{
    label SunRayId(-1);
    scalar maxSunRay = -GREAT;

    // Looking for the ray closest to the Sun direction
    forAll(IRay_, rayId)
    {
        const vector& iD = IRay_[rayId].d();
        scalar dir = sunDir & iD;
        if (dir > maxSunRay)
        {
            maxSunRay = dir;
            SunRayId = rayId;
        }
    }

    // Second rotation to align colimated radiation with the closest ray
    const tensor coordRot = rotationTensor(IRay_[SunRayId].d(), sunDir);

    forAll(IRay_, rayId)
    {
        IRay_[rayId].dAve() = coordRot & IRay_[rayId].dAve();
        IRay_[rayId].d() = coordRot & IRay_[rayId].d();
    }

    Info<< "Sun direction : " << sunDir << nl << endl;
    Info<< "Sun ray ID : " << SunRayId << nl << endl;
}


void Foam::radiationModels::fvDOM::updateRaysDir()
{
    solarCalculator_->correctSunDirection();
    const vector sunDir = solarCalculator_->direction();

    // First iteration
    if (updateTimeIndex_ == 0)
    {
        rotateInitialRays(sunDir);
        alignClosestRayToSun(sunDir);
    }
    else if (updateTimeIndex_ > 0)
    {
        alignClosestRayToSun(sunDir);
    }
}

void Foam::radiationModels::fvDOM::initialise()
{
    setMaxIter();

    coeffs_.readIfPresent("useExternalBeam", useExternalBeam_);

    if (useExternalBeam_)
    {
        spectralDistributions_.reset
        (
            Function1<scalarField>::New
            (
                "spectralDistribution",
                coeffs_
            )
        );

        spectralDistribution_ =
            spectralDistributions_->value(mesh_.time().timeOutputValue());

        spectralDistribution_ =
            spectralDistribution_/sum(spectralDistribution_);

        const dictionary& solarDict = this->subDict("solarCalculatorCoeffs");
        solarCalculator_.reset(new solarCalculator(solarDict, mesh_));

        if (mesh_.nSolutionD() != 3)
        {
            FatalErrorInFunction
                << "External beam model only available in 3D meshes "
                << abort(FatalError);
        }

        if (solarCalculator_->diffuseSolarRad() > 0)
        {
            FatalErrorInFunction
                << "External beam model does not support Diffuse "
                << "Solar Radiation. Set diffuseSolarRad to zero"
                << abort(FatalError);
        }
        if (spectralDistribution_.size() != nLambda_)
        {
            FatalErrorInFunction
                << "The epectral energy distribution has different bands "
                << "than the absoprtivity model "
                << abort(FatalError);
        }
    }

    // 3D
    if (mesh_.nSolutionD() == 3)
    {
        scalar deltaPhi = 0;
        labelList nPhiOfTheta;
        if (nRayInput_)
        {
            // nTheta and nPhi(Theta)  calculated based on user specified
            // value for nRay_ - See 6.12 and 6.13 FDS Technical Ref.
            nTheta_ = 2.0*round(0.5*1.17*pow(nRay_, (1.0/2.26)));
            nPhiOfTheta.setSize(nTheta_);

            const label nRayUser = nRay_;

            // Reset and compute actual value
            nRay_ = 0;
            for (label i = 1; i <= nTheta_; i++)
            {
                const scalar lowerTheta = pi*(i - 1)/nTheta_;
                const scalar upperTheta = pi*i/nTheta_;
                nPhi_ =
                    4.0*round
                    (
                        0.25*max
                        (
                            4.0,
                            round
                            (
                                0.5*nRayUser
                              *(Foam::cos(lowerTheta) - Foam::cos(upperTheta))
                            )
                        )
                    );

                nPhiOfTheta[i - 1] = nPhi_;
                nRay_ += nPhi_;
            }
        }
        else
        {
            nRay_ = 4*nPhi_*nTheta_;
            deltaPhi = pi/(2.0*nPhi_);
        }

        IRay_.setSize(nRay_);
        const scalar deltaTheta = pi/nTheta_;
        label i = 0;
        for (label n = 1; n <= nTheta_; n++)
        {
            const scalar thetai = (2*n - 1)*deltaTheta/2.0;
            const label nPhi = nRayInput_ ? nPhiOfTheta[n - 1] : 4*nPhi_;
            deltaPhi = nRayInput_ ? 2*pi/nPhiOfTheta[n - 1] : deltaPhi;

            for (label m = 1; m <= nPhi; m++)
            {
                const scalar phii = (2*m - 1)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_,
                        blackBody_,
                        i
                   )
               );
               i++;
            }
        }
    }
    // 2D
    else if (mesh_.nSolutionD() == 2)
    {
        const scalar thetai = piByTwo;
        const scalar deltaTheta = pi;

        label iterSize;
        scalar deltaPhi;

        // FDS code equvivalent
        if (nRayInput_)
        {
            nTheta_ = 1;
            nPhi_ = 4.0*round(0.25*nRay_);
            nRay_ = nTheta_*nPhi_;
            deltaPhi = 2.0*pi/nPhi_;
            iterSize = nPhi_;
        }
        else
        {
            nRay_ = 4*nPhi_;
            deltaPhi = pi/(2.0*nPhi_);
            iterSize = 4*nPhi_;
        }
        IRay_.setSize(nRay_);
        label i = 0;
        for (label m = 1; m <= iterSize; m++)
        {
            const scalar phii = (2*m - 1)*deltaPhi/2.0;
            IRay_.set
            (
                i,
                new radiativeIntensityRay
                (
                    *this,
                    mesh_,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta,
                    nLambda_,
                    absorptionEmission_,
                    blackBody_,
                    i
                )
            );
            i++;
        }
    }
    // 1D
    else
    {
        const scalar thetai = piByTwo;
        const scalar deltaTheta = pi;

        if (nRayInput_)
        {
            NotImplemented;
        }

        nRay_ = 2;
        IRay_.setSize(nRay_);
        const scalar deltaPhi = pi;
        label i = 0;
        for (label m = 1; m <= 2; m++)
        {
            const scalar phii = (2*m - 1)*deltaPhi/2.0;
            IRay_.set
            (
                i,
                new radiativeIntensityRay
                (
                    *this,
                    mesh_,
                    phii,
                    thetai,
                    deltaPhi,
                    deltaTheta,
                    nLambda_,
                    absorptionEmission_,
                    blackBody_,
                    i
                 )
            );
            i++;
        }
    }


    // Construct absorption field for each wavelength
    forAll(aLambda_, lambdaI)
    {
        aLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "aLambda_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                a_
            )
        );
    }

    Info<< "fvDOM : Allocated " << IRay_.size()
        << " rays with average orientation:" << nl;

    if (useExternalBeam_)
    {
        // Rotate rays for Sun direction
        updateRaysDir();
    }

    scalar totalOmega = 0;
    forAll(IRay_, rayId)
    {
        if (omegaMax_ <  IRay_[rayId].omega())
        {
            omegaMax_ = IRay_[rayId].omega();
        }
        totalOmega += IRay_[rayId].omega();
        Info<< '\t' << IRay_[rayId].rayName() << " : " << "dAve : "
            << '\t' << IRay_[rayId].dAve() << " : " << "omega : "
            << '\t' << IRay_[rayId].omega() << " : " << "d : "
            << '\t' << IRay_[rayId].d() << nl;
    }

    Info<< "Total omega : " << totalOmega << endl;

    Info<< endl;

    solarLoadMode_ =
        solarLoadTypeNames_
        [
            this->lookupOrDefault<word>("solarLoadMode", "none")
        ];

    switch (solarLoadMode_)
    {
        case sltViewFactor:
        {
            const dictionary& solarDict = this->subDict("solarLoadCoeffs");
            vfSolarLoad_.reset
            (
                new solarLoad(solarDict, T_, TRef_, externalRadHeatFieldName_)
            );

            if (vfSolarLoad_->nBands() > 1)
            {
                FatalErrorInFunction
                    << "Requested solar radiation with fvDOM. Using "
                    << "more than one band for the solar load is not allowed"
                    << exit(FatalError);
            }

            Info<< "Creating Viewfactor Solar Load Model " << nl;
        }
        break;

        case sltDOM:
        {
            const dictionary& solarDict = this->subDict("domSolarCoeffs");

            domSolarLoad_.reset
            (
                new domSolar(solarDict, T_, externalRadHeatFieldName_, TRef_)
            );

            Info<< "Creating DOM Solar Load Model " << nl;
        }
        break;

        case sltnone:
        {}
        break;

        default:
        {
            FatalErrorInFunction
                << "Unsupported solarLoadMode "
                << solarLoadMode_
                << exit(FatalError);
        }
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::fvDOM::fvDOM
(
    const volScalarField& T,
    const dimensionedScalar& TRef
)
:
    radiationModel(typeName, T, TRef),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    qin_
    (
        IOobject
        (
            "qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
//    qem_
//    (
//        IOobject
//        (
//            "qem",
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh_,
//        dimensionedScalar("qem", dimMass/pow3(dimTime), 0)
//    ),
    qg_
    (
        IOobject
        (
            "QrG",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    nTheta_(coeffs_.lookupOrDefault<label>("nTheta", 4)),
    nPhi_(coeffs_.lookupOrDefault<label>("nPhi", 2)),
    nRayInput_(coeffs_.lookupOrDefault<Switch>("nRayInput", false)),
    nRay_(coeffs_.lookupOrDefault<label>("nRay", 0)),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),
    blackBody_(nLambda_, T),
    IRay_(0),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    initConvergence_
    (
        coeffs_.lookupOrDefault<scalar>("initialConvergence", convergence_)
    ),
    initMaxIter_(coeffs_.lookupOrDefault<label>("initialMaxIter", maxIter_)),
    omegaMax_(0),
    useExternalBeam_(false),
    spectralDistribution_(),
    spectralDistributions_(),
    solarCalculator_(),
    solarLoadMode_
    (
        solarLoadTypeNames_
        [
            this->lookupOrDefault<word>("solarLoadMode", "none")
        ]
    ),
    vfSolarLoad_(),
    domSolarLoad_(),
    meshOrientation_(coeffs_.lookupOrDefault<vector>("meshOrientation", Zero)),
    updateTimeIndex_(0)
{
    initialise();
}


Foam::radiationModels::fvDOM::fvDOM
(
    const dictionary& dict,
    const volScalarField& T,
    const dimensionedScalar& TRef
)
:
    radiationModel(typeName, dict, T, TRef),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    qin_
    (
        IOobject
        (
            "qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
//    qem_
//    (
//        IOobject
//        (
//            "qem",
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh_,
//        dimensionedScalar("qem", dimMass/pow3(dimTime), 0)
//    ),
    qg_
    (
        IOobject
        (
            "QrG",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    nTheta_(coeffs_.lookupOrDefault<label>("nTheta", 4)),
    nPhi_(coeffs_.lookupOrDefault<label>("nPhi", 2)),
    nRayInput_(coeffs_.lookupOrDefault<Switch>("nRayInput", false)),
    nRay_(coeffs_.lookupOrDefault<label>("nRay", 0)),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),
    blackBody_(nLambda_, T),
    IRay_(0),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    initConvergence_
    (
        coeffs_.lookupOrDefault<scalar>("initialConvergence", convergence_)
    ),
    initMaxIter_(coeffs_.lookupOrDefault<label>("initialMaxIter", maxIter_)),
    omegaMax_(0),
    useExternalBeam_(false),
    spectralDistribution_(),
    spectralDistributions_(),
    solarCalculator_(),
    solarLoadMode_
    (
        solarLoadTypeNames_
        [
            this->lookupOrDefault<word>("solarLoadMode", "none")
        ]
    ),
    vfSolarLoad_(),
    domSolarLoad_(),
    meshOrientation_(coeffs_.lookupOrDefault<vector>("meshOrientation", Zero)),
    updateTimeIndex_(0)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::fvDOM::~fvDOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiationModels::fvDOM::read()
{
    if (radiationModel::read())
    {
        // Only reading solution parameters - not changing ray geometry
        coeffs_.readIfPresent("convergence", convergence_);
        coeffs_.readIfPresent("maxIter", maxIter_);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiationModels::fvDOM::calculate()
{
    absorptionEmission_->correct(a_, aLambda_);

    updateBlackBodyEmission();

    if (useExternalBeam_)
    {
        switch (solarCalculator_->sunDirectionModel())
        {
            case solarCalculator::mSunDirConstant:
            {
                break;
            }
            case solarCalculator::mSunDirTracking:
            {
                label updateIndex = label
                (
                    mesh_.time().value()
                    /solarCalculator_->sunTrackingUpdateInterval()
                );

                if (updateIndex > updateTimeIndex_)
                {
                    Info<< "Updating Sun position..." << endl;
                    updateTimeIndex_ = updateIndex;
                    updateRaysDir();
                }
                break;
            }
        }
    }

    switch (solarLoadMode_)
    {
        case sltViewFactor:
        {
            vfSolarLoad_->calculate();
        }
        break;

        case sltDOM:
        {
            domSolarLoad_->calculate();
        }
        break;

        case sltnone:
        {}
        break;

        default:
        {
            FatalErrorInFunction
                << "Unsupported solarLoadMode "
                << solarLoadMode_
                << exit(FatalError);
        }
    }

    // Set rays convergence false
    List<bool> rayIdConv(nRay_, false);

    scalar maxResidual = 0;
    label radIter = 0;

    scalar conv = convergence_;
    scalar maxi = maxIter_;

    if (time_.timeIndex() == 0)
    {
        conv = initConvergence_;
        maxi = initMaxIter_;
        Info<< nl <<  "Initial fvDOM solution " << maxi
            << " iterations." << endl;
    }

    do
    {
        Info<< "Radiation solver iter: " << radIter << endl;

        radIter++;
        maxResidual = 0;
        forAll(IRay_, rayi)
        {
            if (!rayIdConv[rayi])
            {
                scalar maxBandResidual = IRay_[rayi].correct();
                maxResidual = max(maxBandResidual, maxResidual);

                reduce(maxBandResidual, maxOp<scalar>());

                if (maxBandResidual < convergence_)
                {
                    rayIdConv[rayi] = true;
                }
            }
        }

        reduce(maxResidual, maxOp<scalar>());
        Info<< "Global Max initial residual: " <<maxResidual << endl;

    } while (maxResidual > conv && radIter < maxi);

    updateG();
}


Foam::tmp<Foam::volScalarField> Foam::radiationModels::fvDOM::Rp() const
{
    tmp<volScalarField> tTabs(T_ + TRef_);
    volScalarField& Tabs = tTabs.ref();

    // Construct using contribution from first frequency band
    tmp<volScalarField> tRp
    (
        volScalarField::New
        (
            "Rp",
            4
           *physicoChemical::sigma
           *(aLambda_[0] - absorptionEmission_->aDisp(0)())
           *blackBody_.deltaLambdaT(Tabs, absorptionEmission_->bands(0))
        )
    );

    volScalarField& Rp = tRp.ref();

    // Add contributions over remaining frequency bands
    for (label j = 1; j < nLambda_; j++)
    {
        Rp +=
            (
                4
               *physicoChemical::sigma
               *(aLambda_[j] - absorptionEmission_->aDisp(j)())
               *blackBody_.deltaLambdaT(Tabs, absorptionEmission_->bands(j))
            );
    }

    return tRp;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiationModels::fvDOM::Ru() const
{
    tmp<DimensionedField<scalar, volMesh>> tRu
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "Ru",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimensionSet(1, -1, -3, 0, 0), 0)
        )
    );

    DimensionedField<scalar, volMesh>& Ru = tRu.ref();

    // Sum contributions over all frequency bands
    for (label j = 0; j < nLambda_; j++)
    {
        // Compute total incident radiation within frequency band
        tmp<DimensionedField<scalar, volMesh>> Gj
        (
            IRay_[0].ILambda(j)()*IRay_[0].omega()
        );

        for (label rayi=1; rayi < nRay_; rayi++)
        {
            Gj.ref() += IRay_[rayi].ILambda(j)()*IRay_[rayi].omega();
        }

        Ru +=
            (aLambda_[j]() - absorptionEmission_->aDisp(j)()())*Gj
          - absorptionEmission_->ECont(j)()();
    }

    return tRu;
}


void Foam::radiationModels::fvDOM::updateBlackBodyEmission()
{
    for (label j=0; j < nLambda_; j++)
    {
        blackBody_.correct(j, absorptionEmission_->bands(j));
    }
}


void Foam::radiationModels::fvDOM::updateG()
{
    G_ = dimensionedScalar(dimMass/pow3(dimTime), 0);
    qr_ = dimensionedScalar(dimMass/pow3(dimTime), 0);
    qg_ = dimensionedScalar(dimMass/pow3(dimTime), 0);
    qin_ = dimensionedScalar(dimMass/pow3(dimTime), 0);

    forAll(IRay_, rayi)
    {
        IRay_[rayi].addIntensity();
        G_ += IRay_[rayi].I()*IRay_[rayi].omega();
        qr_.boundaryFieldRef() += IRay_[rayi].qr();
        qg_.boundaryFieldRef() += IRay_[rayi].qg();
        qin_.boundaryFieldRef() += IRay_[rayi].qin();
    }

    if
    (
        mesh_.foundObject<volScalarField>(externalRadHeatFieldName_)
     && solarLoadMode_ != sltnone
    )
    {
        const volScalarField& Qext =
            mesh_.lookupObject<volScalarField>(externalRadHeatFieldName_);

        forAll(mesh_.boundaryMesh(), patchi)
        {
            if (!boundaryProperties().radBoundaryProperties()[patchi].empty())
            {
                qin_.boundaryFieldRef()[patchi] -=
                    boundaryProperties().emissivity(patchi)
                   *Qext.boundaryField()[patchi];

                qr_.boundaryFieldRef()[patchi] +=
                    (1 - boundaryProperties().transmissivity(patchi))
                   *Qext.boundaryField()[patchi];

                qg_.boundaryFieldRef()[patchi] -= Qext.boundaryField()[patchi];
            }
        }
    }
}


void Foam::radiationModels::fvDOM::setRayIdLambdaId
(
    const word& name,
    label& rayId,
    label& lambdaId
) const
{
    // Assume name is in the form: <name>_<rayId>_<lambdaId>
    const size_type i1 = name.find_first_of("_");
    const size_type i2 = name.find_last_of("_");

    rayId = readLabel(IStringStream(name.substr(i1 + 1, i2 - 1))());
    lambdaId = readLabel(IStringStream(name.substr(i2 + 1, name.size() - 1))());
}


Foam::tmp<Foam::scalarField> Foam::radiationModels::fvDOM::emittedRadiantIntensity
(
    label patchi,
    const scalarField& TpAbs
) const
{
    tmp<scalarField> efPtr(new scalarField(TpAbs.size(), 0));
    scalarField& ef = efPtr.ref();

    forAll(IRay_, rayI)
    {
        ef += IRay_[rayI].emittedRadiantIntensity(patchi, TpAbs);
    }

    return efPtr;
}


const Foam::radiationModels::domSolar& Foam::radiationModels::fvDOM::getDomSolarObj() const
{
    return domSolarLoad_;
}


const Foam::solarCalculator& Foam::radiationModels::fvDOM::solarCalc() const
{
    return solarCalculator_();
}

// ************************************************************************* //
