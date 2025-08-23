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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::fv::limitTemperature

Group
    grpFvOptionsCorrections

Description
    Correction for temperature to apply limits between minimum and maximum
    values.

Usage
    Example usage:
    \verbatim
    limitT
    {
        type            limitTemperature;
        active          yes;
        field           T;
        selectionMode   all;
        min             200;
        max             500;
    }
    \endverbatim

SourceFiles
    limitTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef limitTemperature_H
#define limitTemperature_H

#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class limitTemperature Declaration
\*---------------------------------------------------------------------------*/

class limitTemperature
:
    public cellSetOption
{
protected:

    // Protected data

        //- Energy eqn. field name, default = h
        word TName_;

        //- Minimum temperature limit [K]
        scalar Tmin_;

        //- Maximum temperature limit [K]
        scalar Tmax_;


public:

    //- Runtime type information
    TypeName("limitTemperature");


    // Constructors

        //- Construct from components
        limitTemperature
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        limitTemperature(const limitTemperature&) = delete;


    //- Destructor
    virtual ~limitTemperature()
    {}


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames) override;

        //- Read dictionary
        virtual bool read(const dictionary& dict) override;

        //- Correct the energy field
        virtual void correct(volScalarField& he) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const limitTemperature&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
