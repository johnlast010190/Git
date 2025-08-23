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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2022-2024 Engys Ltd.

Class
    Foam::fv::forceActuationDiskSource

Description
    Force actuation disk (momentum) source with force input based on
    local coordinate system.

SourceFiles
    forceActuationDiskSource.C
    forceActuationDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef forceActuationDiskSource_H
#define forceActuationDiskSource_H

#include "cellSetOption/cellSetOption.H"
#include "coordinate/systems/coordinateSystem.H"
#include "coordinate/systems/cartesianCS.H"
#include "referenceFrames/coordinateFrame.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class forceActuationDiskSource Declaration
\*---------------------------------------------------------------------------*/

class forceActuationDiskSource
:
    public cellSetOption
{
public:

    // Public data

        //- Enumeration for volume types
        enum volumeModeType
        {
            vmAbsolute,
            vmSpecific
        };

        //- Word list of volume mode type names
        static const NamedEnum<volumeModeType, 2> volumeModeTypeNames_;


protected:

    // Protected data

        //- Volume mode
        volumeModeType volumeMode_;

        //- Volume normalisation
        scalar VDash_;

        //- Force
        autoPtr<Function1<vector>> force_;

        //- Local coordinate System
        autoPtr<coordinateSystem> coorSysPtr_;

        //- Reference frame ptr
        coordinateFrame* coorFramePtr_;

        //- Convert force to radial distribution
        bool radialDistribution_;

        //- Access to local coordinate system
        inline const coordinateSystem& csys() const;


private:

    // Private Member Functions

        //- Create coordinate system
        void createCS();

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addForceActuationDiskResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;


public:

    //- Runtime type information
    TypeName("forceActuationDiskSource");


    // Constructors

        //- Construct from components
        forceActuationDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        forceActuationDiskSource
        (
            const forceActuationDiskSource&
        ) = delete;


    //- Destructor
    virtual ~forceActuationDiskSource()
    {}


    // Member Functions

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);

        //- Source term to fvMatrix<vector>
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldI
        );
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldI
        );

        //- Add contribution to block Momentum equation
        virtual void addSup
        (
            fvBlockMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add contribution to block Momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvBlockMatrix<vector>& eqn,
            const label fieldi
        );

        //- Write data
        virtual void writeData(Ostream&) const;

        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const forceActuationDiskSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "forceActuationDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
