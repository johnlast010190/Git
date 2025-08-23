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
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::aerodynamicPower

Group
    grpForcesFunctionObjects

Description
    Calculates the aerodynamicPower of the bins, similar to forces.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    aerodynamicPower.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_aerodynamicPower_H
#define functionObjects_aerodynamicPower_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "coordinate/systems/coordinateSystem.H"
#include "referenceFrames/coordinateFrame.H"
#include "fields/volFields/volFieldsFwd.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "primitives/Tuple2/Tuple2.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "primitives/bools/Switch/Switch.H"
#include "meshes/polyMesh/zones/cellZone/cellZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class aerodynamicPower Declaration
\*---------------------------------------------------------------------------*/

class aerodynamicPower
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Pressure, viscous and porous force per bin
        List<Field<vector>> force_;

        // File streams

            //- Forces
            autoPtr<OFstream> powerFilePtr_;

            //- Force bins
            autoPtr<OFstream> powerBinFilePtr_;

            autoPtr<OFstream> velBinFilePtr_;

            autoPtr<OFstream> powerAvgBinFilePtr_;

        // Read from dictionary

            //- Patches to integrate aerodynamicPower over
            labelHashSet patchSet_;

            //- Name of pressure field
            word pName_;

            //- Name of velocity field
            word UName_;

            //- Name of density field (optional)
            word rhoName_;

            //- Is the force density being supplied directly?
            Switch directForceDensity_;

            //- The name of the force density (fD) field
            word fDName_;

            //- Whether to cache/use cached boundary devRhoReff object for
            //  sharing with other forces FOs
            Switch devRhoReffCaching_;

            //- Controls caching is activated, to avoid using it when we are not
            //  being used as a function object
            bool devRhoReffCachingActive_;

            //- Reference density needed for incompressible calculations
            scalar rhoRef_;

            //- Reference pressure
            scalar pRef_;

            //- Coordinate system used when evaluting aerodynamicPower
            coordinateSystem coordSys_;

            //- Reference Frame Ptr
            coordinateFrame* coorFramePtr_;

            //- Are variables defined in the frame of reference?
            Switch definedInFrame_;

            //- Flag to indicate whether we are using a local co-ordinate sys
            bool localSystem_;

            //- Flag to include porosity effects
            bool porosity_;

            //- Names of zones that will be taken into account (if present)
            wordList porosityZones_;


            // Bin information

                //- Number of bins
                label nBin_;

                //- Direction used to determine bin orientation
                vector binDir_;

                //- Distance between bin divisions
                scalar binDx_;

                //- Minimum bin bounds
                scalar binMin_;

                //- Bin positions along binDir
                List<point> binPoints_;

                //- Should bin data be cumulative?
                bool binCumulative_;

                scalarField aeroPower_;

                scalarField aeroPowerAvg_;

                scalar power_;

                scalar timePower_;

                scalar totalTime_;

                label totalIter_;

                scalar window_;


            //- Write fields flag
            bool writeFields_;

            //- Initialised flag
            bool initialised_;

            //- Write aerodynamicPower in old single file format
            bool combineFiles_;

            scalar timePowerPatch_;
            bool areaWeighted_;

            const word bDevRhoReffName() const
            {
                return typeName+":bDevRhoReff";
            }

    // Protected Member Functions

        //- Local coordinate system
        const coordinateSystem& csys() const;

        //- Frame bin direction
        virtual vector binDir() const;

        //- Create a field name
        word fieldName(const word& name) const;

        //- Create the output files
        virtual void createFiles();

        //- Write header for old style aerodynamicPower.dat
        void writeIntegratedHeader();

        // //- Write header for integrated data
        // virtual void writeIntegratedHeader
        // (
        //     const word& header,
        //     Ostream& os
        // ) const;

        //- Write header for binned data (old format)
        void writeBinHeader();

        //- Write header for binned data
        // virtual void writeBinHeader(const word& header, Ostream& os) const;

        //- Initialise the fields
        void initialise();

        //- Initialise the collection bins
        void initialiseBins();

        //- Reset the fields prior to accumulation of force/moments
        void resetFields();

        //- Reset the fields prior to accumulation of force/moments
        bool checkZone(const cellZone&);

        //- Return effective boundary viscous stress (laminar + turbulent)
        virtual tmp<volSymmTensorField::Boundary> bDevRhoReff() const;

        //- Calculate or lookup cached 'compressible' boundary viscous stress
        tmp<volSymmTensorField::Boundary> cachedBDevRhoReff() const;

        //- Dynamic viscosity field
        tmp<volScalarField> mu() const;

        //- Return rho if specified otherwise rhoRef
        tmp<volScalarField> rho() const;

        //- Return rhoRef if the pressure field is dynamic, i.e. p/rho
        //  otherwise return 1
        scalar rho(const volScalarField& p) const;

        //- return p
        virtual const tmp<volScalarField> getp();

        //- Accumulate bin data
        void applyBins
        (
            const vectorField& fN,
            const vectorField& fT,
            const vectorField& fP,
            const vectorField& d
        );

        void integrateOnPatch
        (
            const vectorField& fN,
            const vectorField& fT,
            const vectorField& fP,
            label patchID
        );

        //- Add patch contributions to force and moment fields
        void addToFields
        (
            const label patchi,
            const vectorField& fN,
            const vectorField& fT,
            const vectorField& fP
        );

        //- Add cell contributions to force and moment fields
        void addToFields
        (
            const labelList& cellIDs,
            const vectorField& fN,
            const vectorField& fT,
            const vectorField& fP
        );

        //- Helper function to write integrated forces and moments (old format)
        void writeIntegratedAeroPower();

        //- Write force data
        void writeAeroPower();

        //- Helper function to write binned forces and moments (old format)
        void writeBinnedAeroPower();

        autoPtr<vectorField> solidVelocityBin();

        autoPtr<vectorField> solidVelocity(label patchID);

        autoPtr<scalarField> calculateMeanField
        (
            const scalarField& baseField
        );

        scalar calculateBinsIntegral(const scalarField& baseField);

        scalar aeroPower() const
        {
            return power_;
        };

public:

    //- Runtime type information
    TypeName("aerodynamicPower");


    // Constructors

        //- Construct from Time and dictionary
        aerodynamicPower
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool baseClass = false
        );

        //- Construct from objectRegistry and dictionary
        aerodynamicPower
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool baseClass = false
        );

        //- Disallow default bitwise copy construction
        aerodynamicPower(const aerodynamicPower&) = delete;


    //- Destructor
    virtual ~aerodynamicPower();


    // Member Functions

        //- Read the aerodynamicPower data
        virtual bool read(const dictionary&);

        //- Calculate the forces and moments
        virtual void calcForcesMoment();

        //- Return the total force
        virtual vector forceEff() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the aerodynamicPower
        virtual bool write();

        //- Called when time step is incremented
        virtual void timeIncremented();


    // Access Functions

        //- Reference pressure
        const scalar& pRef()
        {
            return pRef_;
        }

        //- Reference density
        const scalar& rhoRef()
        {
            return rhoRef_;
        }

        //- Name of pressure field
        const word& pName()
        {
            return pName_;
        }

        //- Name of velocity field
        const word& UName()
        {
            return UName_;
        }

        //- Name of density field
        const word& rhoName()
        {
            return rhoName_;
        }

        //- Patch set
        labelHashSet& patchSet()
        {
            return patchSet_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const aerodynamicPower&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
