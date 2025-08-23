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
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "liquidProperties/liquidProperties/liquidProperties.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "global/etcFiles/etcFiles.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liquidProperties, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquidProperties::liquidProperties(const dictionary& dict)
:
    thermophysicalProperties(dict),
    name_(extractFluidName(dict)),
    Tc_(dict.lookup<scalar>("Tc")),
    Pc_(dict.lookup<scalar>("Pc")),
    Vc_(dict.lookup<scalar>("Vc")),
    Zc_(dict.lookup<scalar>("Zc")),
    Tt_(dict.lookup<scalar>("Tt")),
    Pt_(dict.lookup<scalar>("Pt")),
    Tb_(dict.lookup<scalar>("Tb")),
    dipm_(dict.lookup<scalar>("dipm")),
    omega_(dict.lookup<scalar>("omega")),
    delta_(dict.lookup<scalar>("delta")),
    pRef_(dict.lookupOrDefault<scalar>("pRef", Pstd)),
    rho_
    (
        Function1<scalar>::New
        (
            "rho",
            dict.optionalSubDict("equationOfStateCoeffs")
        )
    ),
    pv_(Function1<scalar>::New("pv", dict.optionalSubDict("pvModelCoeffs"))),
    hl_(Function1<scalar>::New("hl", dict.optionalSubDict("hlModelCoeffs"))),
    Cp_
    (
        Function1<scalar>::New
        (
            "Cp",
            dict.optionalSubDict("thermodynamicsCoeffs")
        )
    ),
    h_(Function1<scalar>::New("h", dict.optionalSubDict("hModelCoeffs"))),
    Cpg_(Function1<scalar>::New("Cpg", dict.optionalSubDict("CpgModelCoeffs"))),
    B_(Function1<scalar>::New("B", dict.optionalSubDict("BModelCoeffs"))),
    mu_(Function1<scalar>::New("mu", dict.optionalSubDict("muModelCoeffs"))),
    mug_(Function1<scalar>::New("mug", dict.optionalSubDict("mugModelCoeffs"))),
    kappa_
    (
        Function1<scalar>::New
        (
            "kappa",
            dict.optionalSubDict("kappaModelCoeffs")
        )
    ),
    kappag_
    (
        Function1<scalar>::New
        (
            "kappag",
            dict.optionalSubDict("kappagModelCoeffs")
        )
    ),
    sigma_
    (
        Function1<scalar>::New
        (
            "sigma",
            dict.optionalSubDict("sigmaModelCoeffs")
        )
    ),
    D_(dict.found("D") ? new Function2s::APIdiffCoef("D", dict) : nullptr),
    Dfunc1_
    (
        dict.found("Dfunc1")
      ? Function1<scalar>::New("Dfunc1", dict)
      : nullptr
    ),
    hf_(h_->value(Tstd))
{}


Foam::liquidProperties::liquidProperties(const liquidProperties& l)
:
    thermophysicalProperties(l.W()),
    name_(l.name_),
    Tc_(l.Tc_),
    Pc_(l.Pc_),
    Vc_(l.Vc_),
    Zc_(l.Zc_),
    Tt_(l.Tt_),
    Pt_(l.Pt_),
    Tb_(l.Tb_),
    dipm_(l.dipm_),
    omega_(l.omega_),
    delta_(l.delta_),
    rho_(l.rho_.clone()),
    pv_(l.pv_.clone()),
    hl_(l.hl_.clone()),
    Cp_(l.Cp_.clone()),
    h_(l.h_.clone()),
    Cpg_(l.Cpg_.clone()),
    B_(l.B_.clone()),
    mu_(l.mu_.clone()),
    mug_(l.mug_.clone()),
    kappa_(l.kappa_.clone()),
    kappag_(l.kappag_.clone()),
    sigma_(l.sigma_.clone()),
    D_(l.D_.valid() ? new Function2s::APIdiffCoef(*l.D_) : nullptr),
    Dfunc1_(l.Dfunc1_.valid() ? l.Dfunc1_.clone() : nullptr)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::liquidProperties> Foam::liquidProperties::New
(
    const word& name
)
{
    if (debug)
    {
        InfoInFunction << "Constructing liquidProperties" << endl;
    }

    IFstream streamRead(findEtcFile(fileName("thermoData/liquidProperties")));
    const dictionary liquidDict(streamRead);
    if (!liquidDict.isDict(name))
    {
        FatalErrorInFunction
            << "liquidProperties type " << name << " not found" << nl
            << Foam::exit(FatalError);
    }
    return autoPtr<liquidProperties>
    (
        new liquidProperties(liquidDict.subDict(name))
    );
}


void Foam::liquidProperties::setDictValue
(
    dictionary& liquidDict,
    const dictionary& dict,
    const word& key,
    const word& subDictName
)
{
    if (dict.found(key) && subDictName == word::null)
    {
        liquidDict.set(key, dict.lookup<scalar>(key));
    }
    else if (subDictName != word::null)
    {
        if (dict.optionalSubDict(subDictName).found(key))
        {
            liquidDict.subDict(subDictName).set
            (
                key,
                dict.optionalSubDict(subDictName).lookup(key)
            );
        }
    }
}


Foam::dictionary Foam::liquidProperties::mergeDicts
(
    const dictionary& dict,
    const word& name
)
{
    // Merging the default liquid properties with the user specified properties
    IFstream streamRead(findEtcFile(fileName("thermoData/liquidProperties")));
    const dictionary liquidsDict(streamRead);
    dictionary liquidDict
    (
        !liquidsDict.isDict(name)
       ? dict
       : liquidsDict.subDict(name)
    );
    liquidDict.changeKeyword("molWeight", "W");
    wordList keys
    ({
        "W",
        "Tc",
        "Pc",
        "Vc",
        "Zc",
        "Tt",
        "Pt",
        "Tb",
        "dipm",
        "omega",
        "delta",
        "D",
        "Dfunc1"
    });
    wordList subDicts(keys.size(), word::null);

    // Append functional dicts
    subDicts.append
    (
        wordList
        ({
            "equationOfStateCoeffs",
            "pvModelCoeffs",
            "hlModelCoeffs",
            "thermodynamicsCoeffs",
            "hModelCoeffs",
            "CpgModelCoeffs",
            "BModelCoeffs",
            "muModelCoeffs",
            "mugModelCoeffs",
            "kappaModelCoeffs",
            "kappagModelCoeffs",
            "sigmaModelCoeffs"
        })
    );
    keys.append
    (
        wordList
        ({
            "rho",
            "pv",
            "hl",
            "Cp",
            "h",
            "Cpg",
            "B",
            "mu",
            "mug",
            "kappa",
            "kappag",
            "sigma"
        })
    );
    forAll(keys, i)
    {
        setDictValue(liquidDict, dict, keys[i], subDicts[i]);
    }
    return liquidDict;
}


Foam::autoPtr<Foam::liquidProperties> Foam::liquidProperties::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing liquidProperties" << endl;
    }

    // If the type is not specified use the entry name as the liquid type name
    const word name
    (
        dict.found("type") ? dict.lookup("type") : dict.dictName()
    );

    return autoPtr<liquidProperties>
    (
        new liquidProperties(mergeDicts(dict, name))
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word Foam::liquidProperties::extractFluidName
(
    const dictionary& dict
) const
{
    string liquidName(fileName::name(dict.dictName()));
    liquidName.removeStart("liquidProperties.");
    return liquidName;
}

const Foam::word& Foam::liquidProperties::name() const
{
    return name_;
}


Foam::scalar Foam::liquidProperties::s(scalar p, scalar T) const
{
    NotImplemented;
}


Foam::scalar Foam::liquidProperties::pvInvert(scalar p) const
{
    // Check for critical and solid phase conditions
    if (p >= Pc_)
    {
        return Tc_;
    }
    else if (p < Pt_)
    {
        if (debug)
        {
            WarningInFunction
                << "Pressure below triple point pressure: "
                << "p = " << p << " < Pt = " << Pt_ <<  nl << endl;
        }
        return -1;
    }

    // Set initial upper and lower bounds
    scalar Thi = Tc_;
    scalar Tlo = Tt_;

    // Initialise T as boiling temperature under normal conditions
    scalar T = Tb_;

    while ((Thi - Tlo) > 1.0e-4)
    {
        if ((pv(p, T) - p) <= 0)
        {
            Tlo = T;
        }
        else
        {
            Thi = T;
        }

        T = (Thi + Tlo)*0.5;
    }

    return T;
}


void Foam::liquidProperties::readIfPresent(const dictionary& dict)
{
    thermophysicalProperties::readIfPresent(dict);
    dict.readIfPresent("Tc", Tc_);
    dict.readIfPresent("Pc", Pc_);
    dict.readIfPresent("Vc", Vc_);
    dict.readIfPresent("Zc", Zc_);
    dict.readIfPresent("Tt", Tt_);
    dict.readIfPresent("Pt", Pt_);
    dict.readIfPresent("Tb", Tb_);
    dict.readIfPresent("dipm", dipm_);
    dict.readIfPresent("omega", omega_);
    dict.readIfPresent("delta", delta_);
}


void Foam::liquidProperties::write(Ostream& os) const
{
    thermophysicalProperties::write(os);

    os.writeEntry("Tc", Tc_);
    os.writeEntry("Pc", Pc_);
    os.writeEntry("Vc", Vc_);
    os.writeEntry("Zc", Zc_);
    os.writeEntry("Tt", Tt_);
    os.writeEntry("Pt", Pt_);
    os.writeEntry("Tb", Tb_);
    os.writeEntry("dipm", dipm_);
    os.writeEntry("omega", omega_);
    os.writeEntry("delta", delta_);

    os << nl;
    rho_->writeData(os); os << nl;
    pv_->writeData(os); os << nl;
    hl_->writeData(os); os << nl;
    Cp_->writeData(os); os << nl;
    h_->writeData(os); os << nl;
    Cpg_->writeData(os); os << nl;
    B_->writeData(os); os << nl;
    mu_->writeData(os); os << nl;
    mug_->writeData(os); os << nl;
    kappa_->writeData(os); os << nl;
    kappag_->writeData(os); os << nl;
    sigma_->writeData(os); os << nl;
    D_->writeData(os); os << endl;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const liquidProperties& l)
{
    l.write(os);
    return os;
}


// ************************************************************************* //
