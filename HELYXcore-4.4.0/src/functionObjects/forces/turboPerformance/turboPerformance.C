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
    (c) D. Boger, B. Lewis M. Auvinen, H. Nilsson
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "turboPerformance/turboPerformance.H"
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"
#include "db/IOstreams/Pstreams/Pstream.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/unitConversion/unitConversion.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "fields/volFields/volFields.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "materialModels/materialTables/materialTables.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(turboPerformance, 0);
    addToRunTimeSelectionTable(functionObject, turboPerformance, dictionary);
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::turboPerformance::machineType,
    5
>::names[] =
{
    "pump",
    "hydraulicReactionTurbine",
    "fan",
    "compressor",
    "gasTurbine"
};

const Foam::NamedEnum<Foam::functionObjects::turboPerformance::machineType, 5>
    Foam::functionObjects::turboPerformance::machineTypeNames_;


template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::turboPerformance::weightingType,
    2
>::names[] =
{
    "area",
    "flux"
};

const Foam::NamedEnum<Foam::functionObjects::turboPerformance::weightingType, 2>
    Foam::functionObjects::turboPerformance::weightingTypeNames_;

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::turboPerformance::transformSectorQuantitiesTo360
(
    vector& totForce,
    vector& totMoment
) const
{
    if (movingWallSectors_>1)
    {
        //- Compute total (360) force and moments

        vector omega = coorFramePtr_->Omega();
        scalar magOmega = mag(omega);
        if (magOmega<SMALL) return;
        vector axis = omega/magOmega;
        const tensor T(axis*axis);
        scalar rotationalAngle = degToRad(360/scalar(movingWallSectors_));

        const tensor S
        (
            0, -axis.z(), axis.y(),
            axis.z(), 0, -axis.x(),
            -axis.y(), axis.x(), 0
        );

        const tensor revTPos
        (
            T
          + cos(rotationalAngle)*(tensor::I - T)
          + sin(rotationalAngle)*S
        );

        vector forceRev(totForce);
        vector sumForce(forceRev);

        vector momRev(totMoment);
        vector sumMom(momRev);

        for (int iSec=2; iSec<=movingWallSectors_; iSec++)
        {
            forceRev = transform(revTPos, forceRev);
            sumForce += forceRev;

            momRev = transform(revTPos, momRev);
            sumMom += momRev;
        }
        totForce = sumForce;
        totMoment = sumMom;
    }
}


void Foam::functionObjects::turboPerformance::writeFileHeaderPumpOrTurbine
(
    Ostream& os
) const
{
    os << "Density (kg/m^3)" << tab
       << "Gravity (m/s^2)" << tab
       << "Omega (rad/s)" << tab
       << "Head (m)" <<tab
       << "Shaft Power (W)" << tab
       << "Efficiency (%)" << tab
       << "F_x (N)" << tab << "F_y (N)" << tab << "F_z (N)" << tab
       << "M_x (Nm)" << tab << "M_y (Nm)" << tab << "M_z (Nm)"
       << endl;
}


void Foam::functionObjects::turboPerformance::writeFileHeaderFan
(
    Ostream& os
) const
{
    os << "Density (kg/m^3)" << tab
       << "Omega (rad/s)" << tab
       << "DPs (Pa)" << tab
       << "DPt (Pa)" << tab
       << "Shaft Power (W)" << tab
       << "Static Efficiency (%)" << tab
       << "Total Efficiency (%)" << tab
       << "F_x (N)" << tab << "F_y (N)" << tab << "F_z (N)" << tab
       << "M_x (Nm)" << tab << "M_y (Nm)" << tab << "M_z (Nm)"
       << endl;
}


void Foam::functionObjects::turboPerformance::writeFileHeaderCompGasTurbine
(
    Ostream& os
) const
{
    os << "Density (kg/m^3)" << tab
       << "Omega (rad/s)" << tab
       << "pRatio_ts " << tab
       << "pRatio_tt" << tab
       << "Shaft Power (W)" << tab
       << "Efficiency_isen_ts (%)" << tab
       << "Efficiency_isen_tt (%)" << tab
       << "F_x (N)" << tab << "F_y (N)" << tab << "F_z (N)" << tab
       << "M_x (Nm)" << tab << "M_y (Nm)" << tab << "M_z (Nm)"
       << endl;
}

Foam::scalar Foam::functionObjects::turboPerformance::entropy
(
    fluidThermo& thermo, const scalar& pp, const scalar& tt
)
{
    volScalarField& T = thermo.T();
    volScalarField& p = thermo.p();
    scalar Told = T[0];
    scalar pold = p[0];
    T[0] = (tt + thermo.TRefValue());
    p[0] = (pp + thermo.pRefValue());
    const scalar s = thermo.materials()(sModel::typeName)[0];
    T[0] = Told;
    p[0] = pold;
    return s;
}

Foam::scalar Foam::functionObjects::turboPerformance::enthalpy
(
    fluidThermo& thermo, const scalar& pp, const scalar& tt
)
{
    volScalarField& T = thermo.T();
    volScalarField& p = thermo.p();
    scalar Told = T[0];
    scalar pold = p[0];
    T[0] = (tt + thermo.TRefValue());
    p[0] = (pp + thermo.pRefValue());
    const scalar h = thermo.materials()(heModel::typeName)[0];
    T[0] = Told;
    p[0] = pold;
    return h;
}

Foam::scalar Foam::functionObjects::turboPerformance::dSdT
(
    fluidThermo& thermo, const scalar& pp, const scalar& tt, const scalar& ee
)
{
    scalar ttm = tt-ee;
    scalar ttp = tt+ee;
    scalar sm = entropy(thermo, pp, ttm);
    scalar sp = entropy(thermo, pp, ttp);
    return (sp-sm)/(2*ee);
}


void Foam::functionObjects::turboPerformance::calculateStagnation
(
    fluidThermo& thermo,
    const scalar& p,
    const scalar& t,
    const scalar& ht,
    const scalar& s,
    scalar& p0,
    scalar& T0
)
{
    const scalar ep = 1e-3;
    const scalar et = 1e-3;
    const label maxIter = 50;
    const scalar tolp = 1e-3;
    const scalar tolt = 1e-4;
    for (label i = 0; i < maxIter; ++i)
    {
        scalar h_curr = enthalpy(thermo, p0, T0);
        scalar s_curr = entropy(thermo, p0, T0);

        //- FD derivatives
        scalar dh_dp = (enthalpy(thermo, p0 + ep, T0) - h_curr)/ep;
        scalar dh_dT = (enthalpy(thermo, p0, T0 + et) - h_curr)/et;
        scalar ds_dp = (entropy(thermo, p0 + ep, T0) - s_curr)/ep;
        scalar ds_dT = (entropy(thermo, p0, T0 + et) - s_curr)/et;

        scalar F1 = h_curr - ht;
        scalar F2 = s_curr - s;

        scalar det = dh_dp * ds_dT - dh_dT * ds_dp;

        //- delta updates
        scalar dp0 = (-F1 * ds_dT + F2 * dh_dT)/det;
        scalar dT0 = (-dh_dp * F2 + ds_dp * F1)/det;

        //- update (maybe using a step?)
        p0 += dp0;
        T0 += dT0;

        if (false)
        {
            Info<< i << tab << dp0 << tab
                 << dT0 << tab << p0 << tab << T0 << endl;
        }

        if (mag(dp0) < tolp && mag(dT0) < tolt) break;
    }
}


void Foam::functionObjects::turboPerformance::calculateTis
(
    fluidThermo& thermo,
    const scalar& p2,
    const scalar& s1,
    const label& maxIter,
    const scalar& tolerance,
    scalar& T2ist
)
{
    scalar error = 1.0;
    for (label iter = 0; iter < maxIter && error > tolerance; iter++)
    {
        // Evaluate entropy at pressure p2 and current T2is
        // The function may vary depending on the thermo model.
        scalar s02 = entropy(thermo, p2, T2ist);

        // minimize f - enforce isentropic condition
        scalar f = s02 - s1;

        // Evaluate the derivative of entropy with respect to
        // temperature at (p2, T2is) using FD
        scalar dsdT = dSdT(thermo, p2, T2ist, 1e-3);

        // Update T2is using N-R:
        scalar dT = -f/dsdT;
        T2ist += dT;

        error = mag(f);

        if (false)
        {
            Info<< "Iteration: " << iter
                 << " T2ist: " << T2ist
                 << " Residual: " << error << endl;
        }
    }
}


void Foam::functionObjects::turboPerformance::executePumpOrTurbine
(
    const vector& totForce,
    const vector& totMoment,
    const vector& omega,
    const scalar& TOmega
)
{

    //- work and hydroHead
    dEmHead dEmH(scalar(0), scalar(0));

    // Pump Efficiency (%)
    scalar eff(0);

    if (machineType_ == pump)
    {
        dEmH = fluidPower::calcDEmHead(false);
        if (mag(TOmega) > 0) eff = (dEmH.first()/TOmega)*scalar(100);
    }
    else if (machineType_ == hydraulicReactionTurbine)
    {
        dEmH = fluidPower::calcDEmHead(true);
        if (mag(dEmH.first()) > 0) eff = (TOmega/dEmH.first())*scalar(100);
    }

    scalar rhoave(rho()().weightedAverage(mesh_.V()).value());
    scalar grav(9.81);
    if (mesh_.foundObject<uniformDimensionedVectorField>("g"))
    {
        const uniformDimensionedVectorField& g =
            mesh_.lookupObject<uniformDimensionedVectorField>("g");

        grav = mag(g[0]);
    }

    // Tab separated output ... to avoid those irritating parenthesis. -- mikko
    file() << obr_.time().timeName() << tab
           << volflux_ << tab
           << rhoave << tab
           << grav << tab
           << mag(omega) << tab
           << dEmH.second() << tab
           << TOmega << tab
           << eff << tab
           << totForce[0] << tab << totForce[1] << tab << totForce[2] << tab
           << totMoment[0] << tab << totMoment[1] << tab << totMoment[2]
           << endl;


    const volScalarField& p
        = lookupObject<volScalarField>(pName_);
    if (p.dimensions() == dimPressure)
    {
        Log << "    Flow (kg/s)      = " << volflux_ << nl;
    }
    else
    {
        Log << "    Flow (m^3/s)     = " << volflux_ << nl;
    }

        Log << "    Density (kg/m^3) = " << rhoave << nl
            << "    Gravity (m/s^2)  = " << grav << nl
            << "    Omega (rad/s)    = " << mag(omega) << nl
            << "    Head (m)         = " << dEmH.second()  << nl
            << "    TOmega (W)       = " << TOmega << nl
            << "    Eff (%)          = " << eff << nl
            << "    Forces (N)       = " << totForce << nl
            << "    Moments (Nm)     = " << totMoment << nl
            << endl;

    // Write state/results information
    setResult("Flow", volflux_);
    setResult("Density", rhoave);
    setResult("Gravity", grav);
    setResult("Omega", mag(omega));
    setResult("Head", dEmH.second());
    setResult("TOmega", TOmega);
    setResult("Efficiency", eff);
    setResult("Forces", totForce);
    setResult("Moments", totMoment);
}


void Foam::functionObjects::turboPerformance::executeFan
(
    const vector& totForce,
    const vector& totMoment,
    const vector& omega,
    const scalar& TOmega
)
{
    scalar dps(0);
    scalar dpt(0);
    scalar ns(0);
    scalar nt(0);

    const volVectorField& U = lookupObject<volVectorField>(UName_);
    const volScalarField& p = lookupObject<volScalarField>(pName_);
    const surfaceScalarField& phi = lookupObject<surfaceScalarField>(phiName_);
    tmp<volScalarField> trho = rho();
    const volVectorField::Boundary& Ub = U.boundaryField();
    const volScalarField::Boundary& pb = p.boundaryField();
    const surfaceScalarField::Boundary&  phib = phi.boundaryField();
    const volScalarField::Boundary&  rhob = trho().boundaryField();

    const surfaceScalarField::Boundary& magSfb = mesh_.magSf().boundaryField();

    scalar mflow1(VSMALL);
    scalar mflow2(VSMALL);
    scalar p1(0.0);
    scalar p2(0.0);
    scalar pt1(0.0);
    scalar pt2(0.0);
    scalar rho1(0.0);

    scalar inArea(0.0);
    scalar outArea(0.0);

    // Quantities at inlet
    forAllConstIter(labelHashSet, inletPatchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& patch = mesh_.boundary()[patchi];

        scalarField fpt1(0.5*rhob[patchi]*magSqr(Ub[patchi]));

        scalarField fp1(fpt1.size(), Zero);
        if (p.dimensions() == dimPressure)
        {
            fp1 = pb[patchi];
        }
        else
        {
            fp1 = pb[patchi]*rhob[patchi];
        }

        fpt1 += fp1;

        scalarField patchMagSf(patch.size(), Zero);
        forAll(patch, pfI)
        {
            const label polyFacei = mesh_.polyFacesBf()[patchi][pfI];
            const label bFacei = polyFacei - mesh_.nInternalFaces();

            // FV indices
            const labelUList patches = mesh_.polyBFacePatches()[bFacei];
            const labelUList patchFaces = mesh_.polyBFacePatchFaces()[bFacei];

            forAll(patches, i)
            {
                patchMagSf[pfI] += magSfb[patches[i]][patchFaces[i]];

                if (p.dimensions() == dimPressure)
                {
                    mflow1 += phib[patches[i]][patchFaces[i]];
                }
                else
                {
                    mflow1 +=
                        rhob[patches[i]][patchFaces[i]]
                       *phib[patches[i]][patchFaces[i]];
                }
            }
        }

        rho1 += sum(rhob[patchi]*patchMagSf);

        fpt1 *= patchMagSf;
        fp1 *= patchMagSf;

        pt1 += sum(fpt1);
        p1 += sum(fp1);
        inArea += sum(patchMagSf);
    }

    // Quantities at outlet
    forAllConstIter(labelHashSet, outletPatchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& patch = mesh_.boundary()[patchi];

        scalarField fpt2(0.5*magSqr(Ub[patchi])*rhob[patchi]);

        scalarField fp2(fpt2.size(), Zero);
        if (p.dimensions() == dimPressure)
        {
            fp2 = pb[patchi];
        }
        else
        {
            fp2 = pb[patchi]*rhob[patchi];
        }

        fpt2 += fp2;

        scalarField patchMagSf(patch.size(), Zero);
        forAll(patch, pfI)
        {
            const label polyFacei = mesh_.polyFacesBf()[patchi][pfI];
            const label bFacei = polyFacei - mesh_.nInternalFaces();

            // FV indices
            const labelUList patches = mesh_.polyBFacePatches()[bFacei];
            const labelUList patchFaces = mesh_.polyBFacePatchFaces()[bFacei];

            forAll(patches, i)
            {
                patchMagSf[pfI] += magSfb[patches[i]][patchFaces[i]];

                if (p.dimensions() == dimPressure)
                {
                    mflow2 += phib[patches[i]][patchFaces[i]];
                }
                else
                {
                    mflow2 +=
                        rhob[patches[i]][patchFaces[i]]
                       *phib[patches[i]][patchFaces[i]];
                }
            }
        }

        fpt2 *= patchMagSf;
        fp2 *= patchMagSf;

        pt2 += sum(fpt2);
        p2 += sum(fp2);
        outArea += sum(patchMagSf);
    }

    reduce(
        std::tie
        (
            pt1, pt2,
            p1, p2,
            inArea, outArea,
            mflow1, mflow2, rho1
        ),
        ParallelOp
        <
            sumOp<scalar>, sumOp<scalar>, sumOp<scalar>, sumOp<scalar>,
            sumOp<scalar>, sumOp<scalar>, sumOp<scalar>, sumOp<scalar>,
            sumOp<scalar>
        >{}
    );


    //- total massflow measured at the inlet
    volflux_ = Foam::mag(scalar(inletSectors_*mflow1));

    if ((inArea==0) || (outArea == 0))
    {
        FatalErrorInFunction
            << "In or outlet patch faces of the fan are zero. "
            << abort(FatalError);
    }
    p1 /= inArea;
    p2 /= outArea;
    pt1 /= inArea;
    pt2 /= outArea;
    rho1 /= inArea;

    if (debugInfo_)
    {
        Info<< endl;
        Info<< "//-- Debug turboPerformance  --//" << endl;
        Info<< "pIn:          " << p1 << endl;
        Info<< "pOut:         " << p2 << endl;
        Info<< "ptIn:         " << pt1 << endl;
        Info<< "ptOut:        " << pt2 << endl;
        Info<< "mass/Volflux: " << volflux_ << endl;
        Info<< "areaIn:       " << inArea << endl;
        Info<< "areaOut:      " << outArea << endl;
        Info<< "//-----------------------------//" << endl;
    }

    dps = p2 - p1;
    dpt = pt2 - pt1;
    ns = (volflux_*dps)/TOmega*100;
    nt = (volflux_*dpt)/TOmega*100;

    const scalar magOmega(mag(omega));

    file() << obr_.time().timeName() << tab
           << volflux_ << tab
           << rho1 << tab
           << magOmega << tab
           << dps << tab
           << dpt << tab
           << TOmega << tab
           << ns << tab
           << nt << tab
           << totForce[0] << tab << totForce[1] << tab << totForce[2] << tab
           << totMoment[0] << tab << totMoment[1] << tab << totMoment[2]
           << endl;


    if (p.dimensions() == dimPressure)
    {
        Log << "    Flow (kg/s)      = " << volflux_ << nl;
    }
    else
    {
        Log << "    Flow (m^3/s)     = " << volflux_ << nl;
    }


        Log << "    Density (kg/m^3) = " << rho1 << nl
            << "    Omega (rad/s)    = " << magOmega << nl
            << "    DPs (Pa)         = " << dps << nl
            << "    DPt (Pa)         = " << dpt << nl
            << "    Shaft Power (W)  = " << TOmega << nl
            << "    Eff_s (%)        = " << ns << nl
            << "    Eff_t (%)        = " << nt << nl
            << "    Forces (N)       = " << totForce << nl
            << "    Moments (Nm)     = " << totMoment << nl
            << endl;

    // Write state/results information
    setResult("Flow", volflux_);
    setResult("Density", rho1);
    setResult("Omega", magOmega);
    setResult("DPs", dps);
    setResult("DPt", dpt);
    setResult("TOmega", TOmega);
    setResult("EfficiencyStatic", ns);
    setResult("EfficiencyTotal", nt);
    setResult("Forces", totForce);
    setResult("Moments", totMoment);
}


void Foam::functionObjects::turboPerformance::executeCompGasTurbine
(
    const vector& totForce,
    const vector& totMoment,
    const vector& omega,
    const scalar& TOmega
)
{
    scalar rpts(0);
    scalar rptt(0);
    scalar nts(0);
    scalar ntt(0);

    const volVectorField& U = lookupObject<volVectorField>(UName_);
    const volScalarField& p = lookupObject<volScalarField>(pName_);
    const surfaceScalarField& phi = lookupObject<surfaceScalarField>(phiName_);
    tmp<volScalarField> trho = rho();
    const volVectorField::Boundary& Ub = U.boundaryField();
    const volScalarField::Boundary& pb = p.boundaryField();
    const surfaceScalarField::Boundary&  phib = phi.boundaryField();
    //const volScalarField::Boundary&  rhob = trho().boundaryField();
    const fluidThermo& thermo =
        lookupObject<fluidThermo>(fluidThermo::dictName);

    const volScalarField& T = thermo.T();
    const volScalarField::Boundary& Tb = T.boundaryField();

    const volScalarField& psi = thermo.psi();
    const volScalarField::Boundary& psib = psi.boundaryField();

    tmp<volScalarField> tGamma(thermo.gamma());

    const volScalarField& Cp = thermo.Cp();
    const volScalarField& gamma = tGamma();
    const volScalarField::Boundary& Cpb = Cp.boundaryField();
    const volScalarField::Boundary& gammab = gamma.boundaryField();
    const volScalarField::Boundary&  rhob = trho().boundaryField();

    const surfaceScalarField::Boundary& magSfb = mesh_.magSf().boundaryField();

    scalar mflow1(VSMALL);
    scalar mflow2(VSMALL);
    scalar p1(0.0);
    scalar p2(0.0);
    scalar pt1(0.0);
    scalar Tt1(0.0);
    scalar pt2(0.0);

    scalar gamma1(0.0);
    scalar cp1(0.0);
    scalar rho1(0.0);

    scalar w1(0.0);
    scalar w2(0.0);

    // Quantities at inlet
    forAllConstIter(labelHashSet, inletPatchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& patch = mesh_.boundary()[patchi];

        scalarField pw(patch.size(), Zero);
        forAll(patch, pfI)
        {
            const label polyFacei = mesh_.polyFacesBf()[patchi][pfI];
            const label bFacei = polyFacei - mesh_.nInternalFaces();

            // FV indices
            const labelUList patches = mesh_.polyBFacePatches()[bFacei];
            const labelUList patchFaces = mesh_.polyBFacePatchFaces()[bFacei];

            forAll(patches, i)
            {
                if (weightingType_ == area)
                {
                    pw[pfI] += magSfb[patches[i]][patchFaces[i]];
                }
                else if (weightingType_ == flux)
                {
                    pw[pfI] += mag(phib[patches[i]][patchFaces[i]]);
                }

                mflow1 += phib[patches[i]][patchFaces[i]];
            }
        }

        scalarField fpt1(patch.size(), 0.0);
        scalarField ftt1(patch.size(), 0.0);

        const fvPatchField<scalar>& gammap = gammab[patchi];
        const fvPatchField<scalar>& psip = psib[patchi];
        const fvPatchField<vector>& Up = Ub[patchi];
        forAll(fpt1, fI)
        {
            scalar gM1 = gammap[fI]-1;
            //- Sanity Check g > 0
            if (gammap[fI]>0)
            {
                scalar gM1ByGI = gM1/gammap[fI];
                if (gM1 == 0)
                {
                    //- g=1 transonic
                    fpt1[fI] = (1+0.5*psip[fI]*magSqr(Up[fI]));
                }
                else if (gM1 > 0)
                {
                    //- g>1 supersonic
                    fpt1[fI] = pow
                        (
                            (1+0.5*psip[fI]*gM1ByGI*magSqr(Up[fI])),
                            1/gM1ByGI
                        );
                }
                //- g>1 supersonic
                ftt1[fI] = (1+0.5*psip[fI]*gM1ByGI*magSqr(Up[fI]));
            }
        }
        scalarField fp1(pb[patchi]);
        fpt1 *= fp1;
        ftt1 *= Tb[patchi];

        gamma1 += sum(gammab[patchi]*pw);
        cp1 += sum(Cpb[patchi]*pw);

        fpt1 *= pw;
        ftt1 *= pw;
        fp1 *= pw;
        rho1 += sum(rhob[patchi]*pw);

        pt1 += sum(fpt1);
        p1 += sum(fp1);
        Tt1 += sum(ftt1);
        w1 += sum(pw);
    }

    // Quantities at outlet
    forAllConstIter(labelHashSet, outletPatchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& patch = mesh_.boundary()[patchi];

        scalarField pw(patch.size(), Zero);
        forAll(patch, pfI)
        {
            const label polyFacei = mesh_.polyFacesBf()[patchi][pfI];
            const label bFacei = polyFacei - mesh_.nInternalFaces();

            // FV indices
            const labelUList patches = mesh_.polyBFacePatches()[bFacei];
            const labelUList patchFaces = mesh_.polyBFacePatchFaces()[bFacei];

            forAll(patches, i)
            {
                if (weightingType_ == area)
                {
                    pw[pfI] += magSfb[patches[i]][patchFaces[i]];
                }
                else if (weightingType_ == flux)
                {
                    pw[pfI] += mag(phib[patches[i]][patchFaces[i]]);
                }

                mflow2 += phib[patches[i]][patchFaces[i]];
            }
        }

        scalarField fpt2(patch.size(), 0.0);

        const fvPatchField<scalar>& gammap = gammab[patchi];
        const fvPatchField<scalar>& psip = psib[patchi];
        const fvPatchField<vector>& Up = Ub[patchi];
        forAll(fpt2, fI)
        {
            scalar gM1 = gammap[fI]-1;
            if (gM1 == 0)
            {
                //- g=1 transonic
                fpt2[fI] = (1+0.5*psip[fI]*magSqr(Up[fI]));
            }
            else if (gM1 > 0)
            {
                //- g>1 supersonic
                scalar gM1ByGI = gM1/gammap[fI];
                fpt2[fI] = pow
                    (
                        (1+0.5*psip[fI]*gM1ByGI*magSqr(Up[fI])),
                        1/gM1ByGI
                    );
            }
        }
        scalarField fp2(pb[patchi]);
        fpt2 *= fp2;

        fpt2 *= pw;
        fp2 *= pw;

        pt2 += sum(fpt2);
        p2 += sum(fp2);
        w2 += sum(pw);
    }

    reduce(
        std::tie
        (
            pt1, pt2,
            p1, p2,
            w1, w2,
            mflow1, mflow2,
            Tt1, gamma1,
            cp1, rho1
        ),
        ParallelOp
        <
            sumOp<scalar>, sumOp<scalar>, sumOp<scalar>, sumOp<scalar>,
            sumOp<scalar>, sumOp<scalar>, sumOp<scalar>, sumOp<scalar>,
            sumOp<scalar>, sumOp<scalar>, sumOp<scalar>, sumOp<scalar>
        >{}
    );

    //- total massflow measured at the inlet
    volflux_ = Foam::mag(scalar(inletSectors_*mflow1));

    if (weightingType_ == area)
    {
        if ((w1<SMALL) || (w2<SMALL))
        {
            FatalErrorInFunction
                << "In or outlet patch faces of the machine are zero. "
                << abort(FatalError);
        }
    }
    else if (weightingType_ == flux)
    {
        if (w1<SMALL)
        {
            Info<< "Zero massflow weights at inlet patch" << endl;
            w1 = SMALL;
        }
        if (w2<SMALL)
        {
            Info<< "Zero massflow weights at outlet patch" << endl;
            w2 = SMALL;
        }
    }

    p1 /= w1;
    p2 /= w2;
    pt1 /= w1;
    Tt1 /= w1;
    pt2 /= w2;
    gamma1 /= w1;
    cp1 /= w1;
    rho1 /= w1;

    //- add reference pressure from material library since ratios require
    //  absolute pressures
    scalar pRef = thermo.pRef().value();
    if (mag(pRef) > SMALL)
    {
        p1 += pRef;
        pt1 += pRef;
        pt2 += pRef;
        p2 += pRef;
    }

    if (debugInfo_)
    {
        Info<< endl;
        Info<< "//-- Debug turboPerformance  --//" << endl;
        Info<< "pIn:      " << p1 << endl;
        Info<< "pOut:     " << p2 << endl;
        Info<< "ptIn:     " << pt1 << endl;
        Info<< "ptOut:    " << pt2 << endl;
        Info<< "cpIn:     " << cp1 << endl;
        Info<< "gammaIn:  " << gamma1 << endl;
        Info<< "ttIn:     " << Tt1 << endl;
        Info<< "massflux: " << volflux_ << endl;
        if (weightingType_ == area)
        {
            Info<< "areaIn:   " << w1 << endl;
            Info<< "areaOut:  " << w2 << endl;
        }
        else if (weightingType_ == flux)
        {
            Info<< "massIn:   " << w1 << endl;
            Info<< "massOut:  " << w2 << endl;
        }
        Info<< "//-----------------------------//" << endl;
    }

    if (TOmega > SMALL)
    {
        if (machineType_ == compressor)
        {
            scalar gm1d1 = (gamma1 - 1)/gamma1;

            if (mag(pt1) < SMALL)
            {
                WarningInFunction
                    << "Compressor with pt1 zero " << endl;

                rpts = 0.0;
                rptt = 0.0;
            }
            else
            {
                rpts = p2/pt1;
                rptt = pt2/pt1;
            }

            scalar powers = volflux_*cp1*Tt1*(pow(rpts, gm1d1) - 1);
            scalar powert = volflux_*cp1*Tt1*(pow(rptt, gm1d1) - 1);

            nts = powers/TOmega*100;
            ntt = powert/TOmega*100;
        }
        else if (machineType_ == gasTurbine)
        {
            scalar gdgm1 = (1 - gamma1)/gamma1;

            if (mag(p2) < SMALL)
            {
                WarningInFunction
                    << "gasTurbine with pOut zero " << endl;

                rpts = 0.0;
            }
            else
            {
                rpts = pt1/p2;
            }

            if (mag(pt2) < SMALL)
            {
                WarningInFunction
                    << "gasTurbine with pt2 zero " << endl;

                rptt = 0.0;
            }
            else
            {
                rptt = pt1/pt2;
            }

            scalar powers = volflux_*cp1*Tt1*(1 - pow(rpts, gdgm1));
            scalar powert = volflux_*cp1*Tt1*(1 - pow(rptt, gdgm1));

            if (powers > SMALL) nts = TOmega/powers*100;
            if (powert > SMALL) ntt = TOmega/powert*100;
        }
    }

    const scalar magOmega(mag(omega));

    file() << obr_.time().timeName() << tab
           << volflux_ << tab
           << rho1 << tab
           << magOmega << tab
           << rpts << tab
           << rptt << tab
           << TOmega << tab
           << nts << tab
           << ntt << tab
           << totForce[0] << tab << totForce[1] << tab << totForce[2] << tab
           << totMoment[0] << tab << totMoment[1] << tab << totMoment[2]
           << endl;


    if (p.dimensions() == dimPressure)
    {
        Log << "    Flow (kg/s)      = " << volflux_ << nl;
    }
    else
    {
        Log << "    Flow (m^3/s)     = " << volflux_ << nl;
    }


        Log << "    Density (kg/m^3) = " << rho1 << nl
            << "    Omega (rad/s)    = " << magOmega << nl
            << "    pRatio_ts        = " << rpts << nl
            << "    pRatio_tt        = " << rptt << nl
            << "    Shaft Power (W)  = " << TOmega << nl
            << "    Eff_isen_ts (%)  = " << nts << nl
            << "    Eff_isen_tt (%)  = " << ntt << nl
            << "    Forces (N)       = " << totForce << nl
            << "    Moments (Nm)     = " << totMoment << nl
            << endl;

    // Write state/results information
    setResult("Flow", volflux_);
    setResult("Density", rho1);
    setResult("Omega", magOmega);
    setResult("rPts", rpts);
    setResult("rPtt", rptt);
    setResult("TOmega", TOmega);
    setResult("EfficiencyTS", nts);
    setResult("EfficiencyTT", ntt);
    setResult("Forces", totForce);
    setResult("Moments", totMoment);
}


void Foam::functionObjects::turboPerformance::executeCompGasTurbineRealGas
(
    const vector& totForce,
    const vector& totMoment,
    const vector& omega,
    const scalar& TOmega
)
{
    scalar rpts(0);
    scalar rptt(0);
    scalar nts(0);
    scalar ntt(0);

    const volVectorField& U = lookupObject<volVectorField>(UName_);
    const volScalarField& p = lookupObject<volScalarField>(pName_);
    const surfaceScalarField& phi = lookupObject<surfaceScalarField>(phiName_);
    tmp<volScalarField> trho = rho();
    const volVectorField::Boundary& Ub = U.boundaryField();
    const volScalarField::Boundary& pb = p.boundaryField();
    const surfaceScalarField::Boundary&  phib = phi.boundaryField();
    fluidThermo& thermo =
        lookupObjectRef<fluidThermo>(fluidThermo::dictName);

    const volScalarField& T = thermo.T();
    const volScalarField::Boundary& Tb = T.boundaryField();

    const volScalarField::Boundary&  rhob = trho().boundaryField();

    const volScalarField& he = thermo.he();

    const surfaceScalarField::Boundary& magSfb = mesh_.magSf().boundaryField();

    scalar mflow1(VSMALL);
    scalar mflow2(VSMALL);

    //- Inlet = 1   ---   Outlet = 2

    scalar p1(0.0);
    scalar p2(0.0);

    scalar T1(0.0);
    scalar T2(0.0);

    scalar h1(0.0);
    scalar h2(0.0);

    scalar ht1(0.0);
    scalar ht2(0.0);

    scalar rho1(0.0);

    scalar w1(0.0);
    scalar w2(0.0);

    // Quantities at inlet
    forAllConstIter(labelHashSet, inletPatchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& patch = mesh_.boundary()[patchi];

        scalarField pw(patch.size(), Zero);
        forAll(patch, pfI)
        {
            const label polyFacei = mesh_.polyFacesBf()[patchi][pfI];
            const label bFacei = polyFacei - mesh_.nInternalFaces();

            // FV indices
            const labelUList patches = mesh_.polyBFacePatches()[bFacei];
            const labelUList patchFaces = mesh_.polyBFacePatchFaces()[bFacei];

            forAll(patches, i)
            {
                if (weightingType_ == area)
                {
                    pw[pfI] += magSfb[patches[i]][patchFaces[i]];
                }
                else if (weightingType_ == flux)
                {
                    pw[pfI] += mag(phib[patches[i]][patchFaces[i]]);
                }

                mflow1 += phib[patches[i]][patchFaces[i]];
            }
        }

        const fvPatchField<vector>& Up = Ub[patchi];
        const fvPatchField<scalar>& hep = he.boundaryField()[patchi];

        scalarField fp1(pb[patchi]);
        scalarField ft1(Tb[patchi]);

        scalarField fh1(hep);
        scalarField fht1(hep);
        fht1 += 0.5*magSqr(Up);

        fp1 *= pw;
        ft1 *= pw;
        fh1 *= pw;
        fht1 *= pw;

        rho1 += sum(rhob[patchi]*pw);

        p1 += sum(fp1);
        T1 += sum(ft1);
        h1 += sum(fh1);
        ht1 += sum(fht1);
        w1 += sum(pw);
    }

    // Quantities at outlet
    forAllConstIter(labelHashSet, outletPatchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& patch = mesh_.boundary()[patchi];

        scalarField pw(patch.size(), Zero);
        forAll(patch, pfI)
        {
            const label polyFacei = mesh_.polyFacesBf()[patchi][pfI];
            const label bFacei = polyFacei - mesh_.nInternalFaces();

            // FV indices
            const labelUList patches = mesh_.polyBFacePatches()[bFacei];
            const labelUList patchFaces = mesh_.polyBFacePatchFaces()[bFacei];

            forAll(patches, i)
            {
                if (weightingType_ == area)
                {
                    pw[pfI] += magSfb[patches[i]][patchFaces[i]];
                }
                else if (weightingType_ == flux)
                {
                    pw[pfI] += mag(phib[patches[i]][patchFaces[i]]);
                }

                mflow2 += phib[patches[i]][patchFaces[i]];
            }
        }

        const fvPatchField<vector>& Up = Ub[patchi];
        const fvPatchField<scalar>& hep = he.boundaryField()[patchi];

        scalarField fp2(pb[patchi]);
        scalarField ft2(Tb[patchi]);

        scalarField fh2(hep);
        scalarField fht2(hep);
        fht2 += 0.5*magSqr(Up);

        fp2 *= pw;
        ft2 *= pw;
        fh2 *= pw;
        fht2 *= pw;

        p2 += sum(fp2);
        T2 += sum(ft2);
        h2 += sum(fh2);
        ht2 += sum(fht2);
        w2 += sum(pw);
    }

    reduce(
        std::tie
        (
            p1, p2,
            T1, T2,
            h1, h2,
            ht1, ht2,
            w1, w2, mflow1, mflow2, rho1
        ),
        ParallelOp
        <
            sumOp<scalar>, sumOp<scalar>, sumOp<scalar>, sumOp<scalar>,
            sumOp<scalar>, sumOp<scalar>, sumOp<scalar>, sumOp<scalar>,
            sumOp<scalar>, sumOp<scalar>, sumOp<scalar>, sumOp<scalar>,
            sumOp<scalar>
        >{}
    );


    //- total massflow measured at the inlet
    volflux_ = Foam::mag(scalar(inletSectors_*mflow1));

    if (weightingType_ == area)
    {
        if ((w1<SMALL) || (w2<SMALL))
        {
            FatalErrorInFunction
                << "In or outlet patch faces of the machine are zero. "
                << abort(FatalError);
        }
    }
    else if (weightingType_ == flux)
    {
        if (w1<SMALL)
        {
            Info<< "Zero massflow weights at inlet patch" << endl;
            w1 = SMALL;
        }
        if (w2<SMALL)
        {
            Info<< "Zero massflow weights at outlet patch" << endl;
            w2 = SMALL;
        }
    }

    p1 /= w1;
    T1 /= w1;
    h1 /= w1;
    ht1 /= w1;
    rho1 /= w1;

    p2 /= w2;
    T2 /= w2;
    h2 /= w2;
    ht2 /= w2;


    //- add reference pressure from material library since ratios require
    //  absolute pressures
    scalar pRef = thermo.pRef().value();
    if (mag(pRef) > SMALL)
    {
        p1 += pRef;
        p2 += pRef;
    }


    // Find p0 T0 such that h(p0, T0) = h_total and s(p0, T0) = s_static
    // 2x2 system
    scalar s1 = entropy(thermo, p1, T1);
    scalar s2 = entropy(thermo, p2, T2);
    scalar pt1 = p1;
    scalar Tt1 = T2;
    calculateStagnation(thermo, p1, T1, ht1, s1, pt1, Tt1);

    scalar pt2 = p1;
    scalar Tt2 = T2;
    calculateStagnation(thermo, p2, T2, ht2, s2, pt2, Tt2);

    //scalar st1 = entropy(thermo, pt1, Tt1);

    // Initial guess
    scalar T2ist = T1;
    scalar T2iss = T1;

    // Find isentropic temperature and enthalpy
    const label maxIter = 100;
    const scalar tolerance = 1e-6;

    calculateTis(thermo, pt2, s1, maxIter, tolerance, T2ist);
    calculateTis(thermo, p2, s1, maxIter, tolerance, T2iss);

    scalar h2ist = enthalpy(thermo, pt2, T2ist);
    scalar h2iss = enthalpy(thermo, p2, T2iss);

    if (debugInfo_)
    {
        Info<< endl;
        Info<< "//-- Debug turboPerformance  --//" << endl;
        Info<< "p1:    " << p1 << endl;
        Info<< "p2:    " << p2 << endl;
        Info<< "pt1:   " << pt1 << endl;
        Info<< "pt2:   " << pt2 << endl;
        Info<< "T1:    " << T1 << endl;
        Info<< "Tt1:   " << Tt1 << endl;
        Info<< "T2iss: " << T2iss << endl;
        Info<< "T2ist: " << T2ist << endl;
        Info<< "h1    " << h1 << endl;
        Info<< "ht1   " << ht1 << endl;
        Info<< "h2    " << h2 << endl;
        Info<< "ht2   " << ht2 << endl;
        Info<< "h2iss " << h2iss << endl;
        Info<< "h2ist " << h2ist << endl;
        Info<< "massflux: " << volflux_ << endl;
        if (weightingType_ == area)
        {
            Info<< "areaIn:   " << w1 << endl;
            Info<< "areaOut:  " << w2 << endl;
        }
        else if (weightingType_ == flux)
        {
            Info<< "massIn:   " << w1 << endl;
            Info<< "massOut:  " << w2 << endl;
        }
        Info<< "//-----------------------------//" << endl;
    }

    if (machineType_ == compressor)
    {
        if (mag(pt1) < SMALL)
        {
            WarningInFunction
                << "Compressor with pt1 zero " << endl;

            rpts = 0.0;
            rptt = 0.0;
        }
        else
        {
            rpts = p2/pt1;
            rptt = pt2/pt1;
        }

        scalar den_s = h2-ht1;
        scalar den_t = ht2-ht1;
        if (mag(den_s) > 0) nts = (h2iss-ht1)/den_s*100;
        if (mag(den_t) > 0) ntt = (h2ist-ht1)/den_t*100;
    }
    else if (machineType_ == gasTurbine)
    {
        if (mag(p2) < SMALL)
        {
            WarningInFunction
                << "gasTurbine with pOut zero " << endl;

            rpts = 0.0;
        }
        else
        {
            rpts = pt1/p2;
        }

        if (mag(pt2) < SMALL)
        {
            WarningInFunction
                << "gasTurbine with pt2 zero " << endl;

            rptt = 0.0;
        }
        else
        {
            rptt = pt1/pt2;
        }


        scalar den_s = ht1-h2iss;
        scalar den_t = ht1-h2ist;

        if (mag(den_s) > 0) nts = (ht1-h2)/den_s*100;
        if (mag(den_t) > 0) ntt = (ht1-ht2)/den_t*100;
    }
    const scalar magOmega(mag(omega));
    file() << obr_.time().timeName() << tab
           << volflux_ << tab
           << rho1 << tab
           << magOmega << tab
           << rpts << tab
           << rptt << tab
           << TOmega << tab
           << nts << tab
           << ntt << tab
           << totForce[0] << tab << totForce[1] << tab << totForce[2] << tab
           << totMoment[0] << tab << totMoment[1] << tab << totMoment[2]
           << endl;


    if (p.dimensions() == dimPressure)
    {
        Log << "    Flow (kg/s)      = " << volflux_ << nl;
    }
    else
    {
        Log << "    Flow (m^3/s)     = " << volflux_ << nl;
    }

        Log << "    Density (kg/m^3) = " << rho1 << nl
            << "    Omega (rad/s)    = " << magOmega << nl
            << "    pRatio_ts        = " << rpts << nl
            << "    pRatio_tt        = " << rptt << nl
            << "    Shaft Power (W)  = " << TOmega << nl
            << "    Eff_isen_ts (%)  = " << nts << nl
            << "    Eff_isen_tt (%)  = " << ntt << nl
            << "    Forces (N)       = " << totForce << nl
            << "    Moments (Nm)     = " << totMoment << nl
            << endl;
    // Write state/results information
    setResult("Flow", volflux_);
    setResult("Density", rho1);
    setResult("Omega", magOmega);
    setResult("rPts", rpts);
    setResult("rPtt", rptt);
    setResult("TOmega", TOmega);
    setResult("EfficiencyTS", nts);
    setResult("EfficiencyTT", ntt);
    setResult("Forces", totForce);
    setResult("Moments", totMoment);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::turboPerformance::turboPerformance
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fluidPower(name, runTime, dict, true),
    machineType_(pump),
    weightingType_(area),
    movingWallSectors_(dict.lookupOrDefault<label>("movingWallSectors", 1)),
    realGas_(dict.lookupOrDefault<Switch>("realGas", false)),
    debugInfo_(dict.lookupOrDefault<Switch>("debugInfo", false))
{
    turboPerformance::read(dict);
    resetFile(turboPerformance::typeName);
    turboPerformance::writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::turboPerformance::~turboPerformance()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::turboPerformance::writeFileHeader(Ostream& os) const
{
    //- common
    const volScalarField& p
        = lookupObject<volScalarField>(pName_);
    os << "# Time" << tab;
    if (p.dimensions() == dimPressure)
    {
        os << "Flow (kg/s)" << tab;
    }
    else
    {
        os << "Flow (m^3/s)" << tab;
    }
    //- Specific to the machine
    if ((machineType_ == pump) || (machineType_ == hydraulicReactionTurbine))
    {
        writeFileHeaderPumpOrTurbine(os);
    }
    else if (machineType_ == fan)
    {
        writeFileHeaderFan(os);
    }
    else
    {
        writeFileHeaderCompGasTurbine(os);
    }
}


bool Foam::functionObjects::turboPerformance::read(const dictionary& dict)
{
    Log << turboPerformance::type() << " " << name() <<  " read:" << nl;

    fluidPower::read(dict);

    if (!dict.found("machineType"))
    {
        WarningInFunction
            << "No machineType entry in dictionary " << dict << endl
            << "Assigning it as a pump" << endl;
    }
    else
    {
        machineType_ = machineTypeNames_.read(dict.lookup("machineType"));
    }

    if (dict.found("weighting"))
    {
        weightingType_ = weightingTypeNames_.read(dict.lookup("weighting"));
    }

    // For now omega (in rad/s) is the only additional info we need
    if (dict.found("omega"))
    {
        FatalErrorInFunction
            << "Turbo performance interface has changed. "
            << "Omega must be specified "
            << "via referenceFrame interface"
            << abort(FatalError);
    }
    if (dict.found("referenceFrame"))
    {
        coorFramePtr_ = coordinateFrame::lookupNew(mesh_, dict);
        definedInFrame_ =
            dict.lookupOrDefault<Switch>("definedInFrame", false);
        localSystem_ = true;
    }
    else
    {
        FatalErrorInFunction
            << "Turbo performance: define the reference frame for the motion. "
            << abort(FatalError);
    }

    // Rotating patch
    patchSet_ =
        mesh_.boundaryMesh().patchSet(wordReList(dict.lookup("patches")));

    return true;
}


bool Foam::functionObjects::turboPerformance::execute()
{
    Log << type() << " " << name() << " execute:" << nl;

    forces::calcForcesMoment();

    vector totForce = forces::forceEff();
    vector totMoment = forces::momentEff();

    // Add sectors quantities to 360
    transformSectorQuantitiesTo360(totForce, totMoment);

    vector omega = coorFramePtr_->Omega();

    // Shaft power (W)
    scalar TOmega = fabs( totMoment & omega);

    //------------------------------------------------------------//
    //- Specific to the machine
    if ((machineType_ == pump) || (machineType_ == hydraulicReactionTurbine))
    {
        executePumpOrTurbine(totForce, totMoment, omega, TOmega);
    }
    else if (machineType_ == fan)
    {
        executeFan(totForce, totMoment, omega, TOmega);
    }
    else
    {
        if (!realGas_)
        {
            executeCompGasTurbine(totForce, totMoment, omega, TOmega);
        }
        else
        {
            executeCompGasTurbineRealGas(totForce, totMoment, omega, TOmega);
        }
    }

    return true;
}


bool Foam::functionObjects::turboPerformance::write()
{
    return true;
}

// ************************************************************************* //
