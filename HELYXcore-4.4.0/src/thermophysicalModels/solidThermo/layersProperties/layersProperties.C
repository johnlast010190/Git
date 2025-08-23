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
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/
#include "layersProperties/layersProperties.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::layersProperties::calculateTotalResistance() const
{
    rtot_ = 0.0;

    if (layerNames_.size()>0)
    {
        forAll(layerNames_, li)
        {
            // sum layer and contact resistances
            rtot_ += (t_[li]/max(VSMALL, kappa_[li]) + rcontact_[li]);
        }
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::layersProperties::layersProperties
(
    const dictionary& dict
)
:
    layerNames_(),
    t_(),
    kappa_(),
    rcontact_(),
    rtot_()
{
    read(dict);

    calculateTotalResistance();
}

Foam::layersProperties::layersProperties
(
    const layersProperties& lp
)
:
    layerNames_(lp.layerNames_),
    t_(lp.t_),
    kappa_(lp.kappa_),
    rcontact_(lp.rcontact_),
    rtot_(lp.rtot_)
{
    // copy constructor
}

// * * * * * * * * * * * * Member Functions  * * * * * * * * * * * //

bool Foam::layersProperties::read(const dictionary& dict)
{

    if (dict.found("layers"))
    {
        const dictionary& layerDicts(dict.subDict("layers"));

        layerNames_.setSize(layerDicts.size());

        t_.setSize(layerDicts.size());
        kappa_.setSize(layerDicts.size());
        rcontact_.setSize(layerDicts.size());

        label nLayers = 0;

        forAllConstIter(dictionary, layerDicts, iter)
        {

            if (!iter().isDict())
            {
                continue;
            }

            const word& key = iter().keyword();
            const dictionary& ldict = iter().dict();

            layerNames_[nLayers] = key;

            t_[nLayers] = ldict.lookup<scalar>("thickness");
            kappa_[nLayers] = ldict.lookup<scalar>("kappa");
            rcontact_[nLayers] = ldict.lookup<scalar>("rcontact");

            nLayers++;
        }
    }

    return true;

}

void Foam::layersProperties::write(Ostream& os) const
{
    os.beginBlock("layers");
    forAll(t_, li)
    {
        os.beginBlock(layerNames_[li]);
        os.writeEntry("thickness", t_[li]);
        os.writeEntry("kappa", kappa_[li]);
        os.writeEntry("rcontact", rcontact_[li]);
        os.endBlock();
    }
    os.endBlock();

}


// ************************************************************************* //
