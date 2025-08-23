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
    (c) 2016-2023 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::moleFractions

Description
    This function object calculates mole-fraction fields from the mass-fraction
    fields of the multicomponent thermo. The names of the mole-fraction fields
    are obtained from the corresponding mass-fraction fields prepended by "X_".

    Example of function object specification:
    \verbatim
    moleFractions
    {
        type moleFractions;
    }
    \endverbatim

    Optionally, the name of the phase can be specified for multiphase cases.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    moleFractions.C

\*---------------------------------------------------------------------------*/

#ifndef moleFractions_H
#define moleFractions_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fields/volFields/volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class moleFractions Declaration
\*---------------------------------------------------------------------------*/

class moleFractions
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Optional phase name
        word phaseName_;

        //- Species mole fractions
        PtrList<volScalarField> X_;


public:

    //- Runtime type information
    TypeName("moleFractions");


    // Constructors

        //- Construct from Time and dictionary
        moleFractions(const word& name, const Time& t, const dictionary& dict);

        //- Disallow default bitwise copy construction
        moleFractions(const moleFractions&) = delete;


    //- Destructor
    virtual ~moleFractions();


    // Member Functions

        //- Read the moleFractions data
        virtual bool read(const dictionary&);

        //- Calculate the mole-fraction fields
        virtual bool execute();

        //- The mole-fraction fields auto-write
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const moleFractions&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
