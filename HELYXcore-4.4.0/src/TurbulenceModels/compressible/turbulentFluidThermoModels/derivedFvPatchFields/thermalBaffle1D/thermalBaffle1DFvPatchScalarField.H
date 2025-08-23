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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd

Class
    Foam::compressible::thermalBaffle1DFvPatchScalarField

Description
    This BC solves a steady 1D thermal baffle.

    The solid properties are specify as dictionary. Boundary sources such as
    radiation will be incorporated into the balance. Some under-relaxation might
    be needed on the boundary sources.  Baffle and solid properties need to be
    specified on the master side of the baffle.

Usage
    Example of the boundary condition specification using constant
    solid thermo :

    \verbatim
    <masterPatchName>
    {
        type            thermalBaffle1D;
        samplePatch     <slavePatchName>;

        thickness       uniform 0.005;  // thickness [m]
        qs              uniform 100;    // heat flux [W/m2]

        relaxation      1;
        kappa           boundaryKappa; // Name of the kappa model
        value           uniform 300;
    }

    <slavePatchName>
    {
        type            thermalBaffle1D;
        samplePatch     <masterPatchName>;
        relaxation      1;
    }
    \endverbatim

SourceFiles
    thermalBaffle1DFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef thermalBaffle1DFvPatchScalarField_H
#define thermalBaffle1DFvPatchScalarField_H

#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"
#include "memory/autoPtr/autoPtr.H"
#include "mappedPatches/mappedPolyPatch/mappedPatchBase.H"
#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "materialModels/materialTables/materialTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
            Class thermalBaffle1DFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class thermalBaffle1DFvPatchScalarField
:
    public mappedPatchBase,
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of the temperature field
        word TName_;

        //- Baffle is activated
        bool baffleActivated_;

        //- Baffle thickness [m]
        mutable scalarField thickness_;

        //- Superficial heat source [W/m2]
        mutable scalarField qs_;

        //- Solid dictionary
        dictionary solidDict_;

        //- Cache qr for relaxation
        scalarField QPrevious_;

        //- Relaxation for qr
        scalar QRelaxation_;


    // Private members

        //- Neighbour patch access
        const thermalBaffle1DFvPatchScalarField& neighbour() const;

        //- Return const solid thermo
        const baseModels<scalar>& kappaModel() const;

        //- Return qs from master
        tmp<scalarField> qs() const;

        //- Return thickness from master
        tmp<scalarField> baffleThickness() const;

        //- Is Owner
        bool owner() const;


public:

    //- Runtime type information
    TypeName("thermalBaffle1D");


    // Constructors

        //- Construct from patch and internal field
        thermalBaffle1DFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thermalBaffle1DFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  thermalBaffle1DFvPatchScalarField onto a new patch
        thermalBaffle1DFvPatchScalarField
        (
            const thermalBaffle1DFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        thermalBaffle1DFvPatchScalarField
        (
            const thermalBaffle1DFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new thermalBaffle1DFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        thermalBaffle1DFvPatchScalarField
        (
            const thermalBaffle1DFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new thermalBaffle1DFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Retun kappa (backward compatibility)
        scalar kappa() const;

        //- Retun kappa (backward compatibility)
        bool isLegacy() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
