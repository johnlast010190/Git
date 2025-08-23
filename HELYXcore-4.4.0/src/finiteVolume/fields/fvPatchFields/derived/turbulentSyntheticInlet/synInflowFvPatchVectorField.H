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
    (c) 2017 OpenFOAM Foundation
    (c) 2018-2020 CWE@LAMC (University of Bologna)
    (c) 2024 Engys Ltd.

Class
    Foam::synInflowFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    This boundary condition describes a synthetic inflow condition
    corrected by applying the VBIC procedure.

Usage
    Run inflowGeneration.py and inflowCorrection.py to generate
    inflowGenerationDict and inflowCorrectionDict inside the system directory.
    Example of the boundary condition specification:

    \verbatim
    <patchName>
    {
        type            synInflow;
        origin          (0 0 0);
        zDir            (0 0 1);
        sf_dictName     inflowGenerationDict;
        vc_dictName     inflowCorrectionDict;
        value           uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    synInflowFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef synInflowFvPatchVectorField_H
#define synInflowFvPatchVectorField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
         Class synInflowFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class synInflowFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
        /* ------------------------ BC interface ---------------------------*/
        vector sf_origin;
        vector sf_zDir;
        const word sf_dictName;
        const word vc_dictName;

        /* ------------------------ Inflow generation ---------------------------*/
        scalar sf_cutOff;
        scalarField sf_seeds;
        vectorField sf_tProfile;
        vectorField sf_uProfile;
        vectorField sf_sProfile;
        vectorField sf_lProfile;
        vectorField sf_p;
        vectorField sf_q;
        vectorField sf_k;

        /* ------------------------ Inflow correction ---------------------------*/
        scalar vc_corr;

        vector vc_t1;
        vector vc_t2;

        scalarField vc_LdOp;
        vectorField vc_LoOp;
        scalarField vc_UdOp;
        vectorField vc_UoOp;

        scalarField vc_LUpc;
        scalarField vc_LUpr;

        vectorField vc_nodes;
        scalarField vc_faces;
        vectorField vc_faceCentres;
        vectorField vc_grad1Op;
        vectorField vc_grad2Op;
        vectorField vc_bcEdges;
        vectorField vc_bcNorm;


        // Initializing local variables
        vector      sf_xDir;
        vector      sf_yDir;

        scalar      sf_uSmall;
        scalar      sf_lSmall;

        bool sf_calculateCorrections;
        vectorField sf_faceC;
        scalarField sf_nodalArea;
        vectorField sf_nodesTrans;
        vectorField sf_nodesProps;
        vectorField sf_seedsWeights;
        scalarField sf_faceAreas;
        scalarField sf_facesIds;
        vectorField sf_seedsProps;
        vectorField sf_uTot_c;
        scalar      sf_dxScale;
        label       sf_curTimeIndex;

        // For LU solver
        scalarField ls_Ld;
        vectorField ls_Lo;
        scalarField ls_Ud;
        vectorField ls_Uo;
        scalarField ls_LUpc ;
        scalarField ls_LUpr;

public:

   //- Runtime type information
   TypeName("synInflow");

   // Constructors

        //- Construct from patch and internal field
        synInflowFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        synInflowFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping
        synInflowFvPatchVectorField
        (
            const synInflowFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        synInflowFvPatchVectorField
        (
            const synInflowFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new synInflowFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        synInflowFvPatchVectorField
        (
            const synInflowFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new synInflowFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

private:

        void initializeQuantities
        (
            const IOdictionary& dict1,
            const IOdictionary& dict2
        );

        void fluctuatingField
        (
            scalar intOff,
            label iN,
            scalar time,
            vector& uvw_flu
        );

        void LUSolver
        (
            const scalarField& ls_b,
            label& id1,
            label& id2,
            scalarField& ls_xxn
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// ************************************************************************* //

