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
    (c) 2010-2012 Engys Ltd.
    (c) 2011-2013 OpenFOAM Foundation

Class
    volumetricNurbs

Description
    Creates a Nurbs based morphing method

\*---------------------------------------------------------------------------*/

#ifndef volumetricNurbs_H
#define volumetricNurbs_H

#include "db/dictionary/dictionary.H"
#include "vNurbsDeformation/controlBox/controlBox.H"
#include "vNurbsDeformation/splineBasis/splineBasis.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class volumetricNurbs Declaration
\*---------------------------------------------------------------------------*/

class volumetricNurbs
{
    // Private data

    //- controlBox + basis functions

    controlBox box;

    splineBasis uBasis, vBasis, wBasis;

        scalarField weights;

        //- Rational bools

        bool checkRational, isRational;

        void checkWeights();

    // Private member functions

        //- Rational bases handling

    point getRationalPoint
    (
        const scalar u,
        const scalar v,
        const scalar w
    ) const;

    //- returns the spatial derivative along the u isos
    vector getRationalUDer
    (
        const scalar u,
        const scalar v,
        const scalar w
    ) const;

    //- returns the spatial derivative along the v isos
    vector getRationalVDer
    (
        const scalar u,
        const scalar v,
        const scalar w
    ) const;

    //- returns the spatial derivative along the w isos
    vector getRationalWDer
    (
        const scalar u,
        const scalar v,
        const scalar w
    ) const;

    //- returns potential displacement of point at (u, v, w)
    vector getRationalDisplacement
    (
        const scalar u,
        const scalar v,
        const scalar w
    ) const;

    scalarField getRationalBasesList
    (
        const scalar u,
        const scalar v,
        const scalar w
    ) const;

    inline bool contains(const point& pt) const;

public:

    // Constructors

        //- Construct from dictionary
        volumetricNurbs
        (
            const dictionary& dict
        );

    volumetricNurbs
    (
        point pMin,
        point pMax,
        label nXCPs,
        label XDegree,
        label nYCPs,
        label YDegree,
        label nZCPs,
        label ZDegree
    );

        //- Construct as copy
        volumetricNurbs(const volumetricNurbs& vNurbs);

        // Destructor

        ~volumetricNurbs()
    {};


    // Member Functions


    //- returns point coors at parameters u, v, w
    point getPoint
    (
        const scalar u,
        const scalar v,
        const scalar w
    );

    //- returns the spatial derivative along the u isos
    vector getUDer
    (
        const scalar u,
        const scalar v,
        const scalar w
    );

    //- returns the spatial derivative along the v isos
    vector getVDer
    (
        const scalar u,
        const scalar v,
        const scalar w
    );

    //- returns the spatial derivative along the w isos
    vector getWDer
    (
        const scalar u,
        const scalar v,
        const scalar w
    );

    //- returns potential displacement of point at (u, v, w)
    vector getDisplacement
    (
        const scalar u,
        const scalar v,
        const scalar w
    );

        const scalarField& getWeights() const;

        scalarField getBasesList
    (
        const scalar u,
        const scalar v,
        const scalar w
    );

        void setWeight
        (
            const label& i,
            const label& j,
            const label& k,
            const scalar& w
        );

    //- returns a point res with u = res.x(), v = res,y() and w = res.z()
    //- if p is out of the box then ( -1 -1 -1 ) is returned

    point invert
    (
        const point& p
    );

    controlBox& getBox();

    const controlBox& getBox() const;

        void assignUknots(const scalarField& uKnots);
        void assignVknots(const scalarField& vKnots);
        void assignWknots(const scalarField& wKnots);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
