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

Class
    Foam::functionObjects::objectiveFunctionObject::
    uniformityIndexFunctionObject

Description
    Calculates uniformity index on patch.

SourceFiles
    uniformityIndexFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef uniformityIndexFunctionObject_H
#define uniformityIndexFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class uniformityIndexFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class uniformityIndexFunctionObject
:
    public objectiveFunctionObject
{

protected:

    // Protected data

        //- total objective patch area
        scalar objPatchA_;


private:

    // Private Member Functions

        //- calculate total area of objective patches
        scalar objectivePatchArea() const;

        //- calculate mean normal velocity magnitude
        scalar targetVelocity() const;

        //- Disallow default bitwise copy construct and assignment
        uniformityIndexFunctionObject
        (
            const uniformityIndexFunctionObject&
        );

        void operator=(const uniformityIndexFunctionObject&);


public:

    //- Runtime type information
    TypeName("uniformityIndex");


    // Constructors

        //- Construct from time, dictionary, adjoint solver and use old-style
        //  adjoint file-format
        uniformityIndexFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from time, dictionary and adjoint solver
        uniformityIndexFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~uniformityIndexFunctionObject();


    // Member Functions

        //- read and update from dictionary
        virtual bool read(const dictionary&) override;

        //- Calculate cost function value
        virtual bool execute() override;

        //- Track and write
        virtual bool write() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
