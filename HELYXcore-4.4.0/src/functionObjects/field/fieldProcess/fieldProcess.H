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
    (c) 2010-2024 Engys Ltd.
    (c) 2008 Icon-CG Ltd.

Class
    Foam::functionObjects::fieldProcess

Group
    grpFieldFunctionObjects

Description
    Performs field operations

SourceFiles
    fieldProcess.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldProcess_H
#define functionObjects_fieldProcess_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/Time/Time.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "cfdTools/general/include/fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class fieldProcess Declaration
\*---------------------------------------------------------------------------*/

class fieldProcess
:
    public fvMeshFunctionObject
{
    // Private data

public:

        //- Calculation type options
        enum calcTypes
        {
            ctnone,
            ctptot,
            ctCp,
            ctCpComponents,
            ctpGauge,
            ctCg,
            ctlambda2,
            ctQ,
            ctGcriterion,
            ctvorticity,
            ctenstrophy,
            ctwallShear,
            cty,
            ctyPlus,
            ctFrictionVelocity,
            ctmag,
            ctgrad,
            ctmaggrad,
            ctplus,
            ctminus,
            ctmultiply,
            ctdivide,
            ctmax,
            ctmin,
            ctbound,
            ctpow,
            ctaverage,
            ctdot,
            ctcross,
            ctCo,
            ctdiv,
            ctPe,
            ctUnw,
            ctHtc,
            ctNwf,
            ctDdt,
            ctPatchWrite,
            ctStrainRate,
            ctUrel,
            ctUrtz,
            ctUCartesianToCylindrical,
            ctUCylindricalToCartesian,
            ctHeatFlux,
            ctMach,
            ctSkinFriction,
            ctCylindricalTransform,
            ctProudmanSources
        };

        //- Field type options
        enum opFieldTypes
        {
            ftnone,
            ftscalar,
            ftvector,
            fttensor,
            ftsymmTensor,
            ftsphTensor
        };

private:

            static const NamedEnum<calcTypes, 47> calcTypeNames_;
            static const NamedEnum<opFieldTypes, 6> opFieldTypeNames_;

        // Read from dictonary

            //- input
            PtrList<dictionary> operations_;

        //- Result values
            //actual fields are stored in objectRegistry
            DynamicList<word> resultFieldNames_;
            DynamicList<bool> writeStatus_;

            //Performs checks to see if duplicate fields
            bool checkFields_;


    // Private Member Functions

        //- write named vol field to disk
        void writeField(word);

        //- clear named vol field
        void clearField(word);

        //- clear all fields associated with this function object
        void clearFields();

        //- search for a vol field of specified name and return its type
        opFieldTypes volFieldType(word);

        //- create zero field with boundary field = patchInternalField
        template<class Type>
        tmp<VolField<Type>>
        createNearWallField
        (
            const VolField<Type>& inField
        );

        //- set boundaries to near-cell values
        template<class geoField>
        void setWallsToNearCell(geoField&);

        template<class geoField>
        void setWallsToNearCell(autoPtr<geoField>&);

        //- generate field with zero values exept in the specified patches
        template<class Type>
        tmp<VolField<Type>>
        createPatchOnlyField
        (
            const VolField<Type>& inField,
            const List<wordRe>& patchesToWrite
        );


        //- Read transportProperties for rho
        dimensionedScalar readConstRho(word rhoName);

        //- Unsupported operator warning
        void notImplementedYet(word);

        //- retrieve muEff
        const tmp<volScalarField> getMuEff();

        //- retrieve muEff
        const tmp<volScalarField> getK();

        //- retrieve muEff
        const tmp<volScalarField> getEpsilon();

        //- Calculate shear stress from turbulence models
        const tmp<volVectorField> wallShear();

        //- manipulate and store result fields in registry
        template<class geoField>
        void setField(word, tmp<geoField>, const dictionary&);
        template<class geoField>
        void setField(word, geoField*, const dictionary&);
        template<class geoField>
        void setField(word, geoField&, const dictionary&);

        scalar referenceValue
        (
            const dictionary& dict,
            const word refName,
            const word fieldName,
            const scalar defaultRefValue = -1
        );

        template<class Type, class Type2>
        Type closestReferenceValueParallel
        (
            const word fieldName,
            const dictionary& dict
        );

        //- operations
        void none(const dictionary& dict);
        void ptot(const dictionary& dict);
        void Cp(const dictionary& dict);
        void CpComponents(const dictionary& dict);
        void pGauge(const dictionary& dict);
        void Cg(const dictionary& dict);
        void lambda2(const dictionary& dict);
        void Q(const dictionary& dict);
        void Gcriterion(const dictionary& dict);
        void vorticity(const dictionary& dict);
        void enstrophy(const dictionary& dict);
        void Co(const dictionary& dict);
        void divOp(const dictionary& dict);
        void magOp(const dictionary& dict);
        void gradOp(const dictionary& dict);
        void magGradOp(const dictionary& dict);
        void multiplyOp(const dictionary& dict);
        void maxOp(const dictionary& dict);
        void minOp(const dictionary& dict);
        void powOp(const dictionary& dict);
        void dotOp(const dictionary& dict);
        void crossOp(const dictionary& dict);
        void wallShear(const dictionary& dict);
        void y(const dictionary& dict);
        void yPlus(const dictionary& dict);
        void frictionVelocity(const dictionary& dict);
        void Pe(const dictionary& dict);
        void Unw(const dictionary& dict);
        void wallHTC(const dictionary& dict);
        void nearWallValue(const dictionary& dict);
        void ddtOp(const dictionary& dict);
        void patchWrite(const dictionary& dict);
        void strainRate(const dictionary& dict);
        void Urel(const dictionary& dict);
        void Urtz(const dictionary& dict);
        void UCartesianToCylindrical(const dictionary& dict);
        void UCylindricalToCartesian(const dictionary& dict);
        void wallHeatFlux(const dictionary& dict);
        void Mach(const dictionary& dict);
        void skinFriction(const dictionary& dict);
        void cylindricalTransform(const dictionary& dict);
        void proudmanSources(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("fieldProcess");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        fieldProcess
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        fieldProcess(const fieldProcess&) = delete;


    //- Destructor
    virtual ~fieldProcess();


    // Member Functions

        //- write results
        void writeFields();

        //- Update solver
        void calculate();

        //- Execute, currently does nothing
        virtual bool execute() override;

        //- Sample and write
        virtual bool write() override;

        //- Clean up
        bool end() override;

        //- Read the parameters
        virtual bool read(const dictionary&) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fieldProcess&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
