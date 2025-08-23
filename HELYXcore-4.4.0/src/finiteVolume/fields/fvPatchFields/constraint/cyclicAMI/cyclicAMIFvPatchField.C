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
    (c) 2017-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/solutionControl/helyxCoupledControl/helyxCoupledControl.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


template<class Type>
bool Foam::cyclicAMIFvPatchField<Type>::implicitCoeffs() const
{
    bool doImplicitCoeffs = false;

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    if (mesh.foundObject<helyxCoupledControl>(solutionControl::typeName))
    {
        Switch implicitAMI
        (
            debug::optimisationSwitches().lookupOrDefault<Switch>
            ("implicitNonCoveredAMI", true)
        );
        const word& fieldName = this->internalField().name();

        if
        (
            implicitAMI &&
            (
                (fieldName == ("U")) || fieldName == ("p")
            )
        )
        {
            doImplicitCoeffs = true;
        }
    }
    return doImplicitCoeffs;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::interpPatchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();
    const labelUList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().nbrPatch().faceCells();

    Field<Type> pnf(iField, nbrFaceCells);

    tmp<Field<Type>> tpnf;

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        tmp<vectorField> normalFt = this->patch().nf();
        const vectorField& normalF = normalFt();
        Field<tensor> trans(this->size(), Zero);
        forAll(trans, cI)
        {
            trans[cI] = I - 2.0*sqr(normalF[cI]);
        }
        Field<Type> mirrorField
        (
            Foam::transform(trans, this->patchInternalField()())
        );

        tpnf = cyclicAMIPatch_.interpolate(pnf, mirrorField);
    }
    else
    {
        tpnf = cyclicAMIPatch_.interpolate(pnf);
    }

    transform().transform(tpnf.ref(), tpnf());

    return tpnf;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict, dict.found("value")),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{
    if (!isA<cyclicAMIFvPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    if (!dict.found("value") && this->coupled())
    {
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{
    if (!isA<cyclicAMIFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_)
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicAMIFvPatchField<Type>::coupled() const
{
    return cyclicAMIPatch_.coupled();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::patchNeighbourField
(
    const Pstream::commsTypes
) const
{
    if (this->calculatedMode())
    {
        // The boundary values were stored as the patch neighbour field values
        // in evaluate(). We use those stored values instead of the
        // implementation of patchNeighbourField below because a derived type
        // may have overridden it, but will be substituted by this constraint
        // type when used in field calculations
        return *this;
    }
    else
    {
        return interpPatchNeighbourField();
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::patchInternalField() const
{
    const Field<Type>& iField = this->internalField();
    const labelUList& faceCells =
        cyclicAMIPatch_.cyclicAMIPatch().faceCells();

    tmp<Field<Type>> tpnf(new Field<Type>(iField, faceCells));
    Field<Type>& pnf = tpnf.ref();

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        const scalarField& weights = cyclicAMIPatch_.cyclicAMIPatch().
            lowCorrAMIWeights();
        tmp<vectorField> normalFt = this->patch().nf();
        const vectorField& normalF = normalFt();
        Field<tensor> trans(this->size(), Zero);
        forAll(trans, cI)
        {
            trans[cI] = I - 2.0*sqr(normalF[cI]);
        }
        Field<Type> mirrorField(Foam::transform(trans, pnf));

        pnf = (1-weights)*mirrorField + weights*pnf;
    }

    return tpnf;
}


template<class Type>
const Foam::cyclicAMIFvPatchField<Type>&
Foam::cyclicAMIFvPatchField<Type>::nbrPatchField() const
{
    const VolField<Type>& fld =
        static_cast<const VolField<Type>&>
        (
            this->primitiveField()
        );

    return refCast<const cyclicAMIFvPatchField<Type>>
    (
        fld.boundaryField()[cyclicAMIPatch_.nbrPatchID()]
    );
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // We co-opt the patch field to store the neighbour field to be used
    // in field calculations. We can do this since it isn't used for snGrad
    // (which uses the coupled value), and other operators that use the face
    // value will use the coupled interpolated value, rather than the
    // patch value.
    if (this->calculatedMode())
    {
        // This is an exception: If evaluate is called on this base class in
        // calculated mode (which shouldn't really happen, but is allowed
        // for convenience), we avoid a circular assignment
        Field<Type>::operator=(interpPatchNeighbourField());
    }
    else
    {
        Field<Type>::operator=
        (
            this->patchNeighbourField(commsType)
        );
    }

    fvPatchField<Type>::evaluate();
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().nbrPatch().faceCells();

    scalarField pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    // patch Internal field
    Field<scalar> pf(psiInternal,cyclicAMIPatch_.faceCells());

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        scalarField pif(psiInternal, cyclicAMIPatch_.faceCells());
        tmp<vectorField> normalFt = this->patch().nf();
        const vectorField& normalF = normalFt();
        Field<tensor> trans(this->size(), Zero);
        forAll(trans, cI)
        {
            trans[cI] = I - 2.0*sqr(normalF[cI]);
        }
        Field<scalar> mirrorField(Foam::transform(trans, pf));

        pnf = cyclicAMIPatch_.interpolate(pnf, mirrorField);
    }
    else
    {
        pnf = cyclicAMIPatch_.interpolate(pnf);
    }

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, coeffs, pnf);
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().nbrPatch().faceCells();

    Field<Type> pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        Field<Type> pif(psiInternal, cyclicAMIPatch_.faceCells());
        pnf = cyclicAMIPatch_.interpolate(pnf, pif);
    }
    else
    {
        pnf = cyclicAMIPatch_.interpolate(pnf);
    }

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, coeffs, pnf);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    tmp<Field<Type>> tvIC =
        this->coupledFvPatchField<Type>::valueInternalCoeffs(w);
    Field<Type>& vIC = tvIC.ref();

    if
    (
        cyclicAMIPatch_.applyLowWeightCorrection()
      && this->cyclicAMIPatch().cyclicAMIPatch().implicitNonCovered()
    )
    {
        if (implicitCoeffs())
        {
            //- Mimic implicit zero Neumann or Dirichlet
            // with cyclicAMI at non-overlaping faces
            // A bit hacky implementation,
            // Implementation needs revision
            // Only for coupled solver currently
            Type defValue = pTraits<Type>::zero;
            if (this->internalField().name() == "p")
            {
                defValue = pTraits<Type>::one;
            }

            if (cyclicAMIPatch_.cyclicAMIPatch().owner())
            {
                const scalarList& srcWSum =
                    cyclicAMIPatch_.AMI().srcWeightsSum();
                forAll(vIC, fI)
                {
                    scalar w2 = 1-srcWSum[fI];
                    vIC[fI] = srcWSum[fI]*vIC[fI]+w2*defValue;
                }
            }
            else
            {
                const scalarList& tgtWSum =
                    cyclicAMIPatch_.nbrPatch().AMI().tgtWeightsSum();

                forAll(vIC, fI)
                {
                    scalar w2 = 1-tgtWSum[fI];
                    vIC[fI] = tgtWSum[fI]*vIC[fI]+w2*defValue;
                }
            }
        }
    }
    return tvIC;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    tmp<Field<Type>> tvBC =
        this->coupledFvPatchField<Type>::valueBoundaryCoeffs(w);
    Field<Type>& vBC = tvBC.ref();

    if
    (
        cyclicAMIPatch_.applyLowWeightCorrection()
      && this->cyclicAMIPatch().cyclicAMIPatch().implicitNonCovered()
    )
    {
        if (implicitCoeffs())
        {
            Type defValue = pTraits<Type>::zero;
            if (cyclicAMIPatch_.cyclicAMIPatch().owner())
            {
                const labelList& nsrcFaces =
                    cyclicAMIPatch_.AMI().nonOverlapSourceFaces();

                forAll(nsrcFaces, fI)
                {
                    vBC[nsrcFaces[fI]] = defValue;
                }
            }
            else
            {
                const labelList& ntgtFaces = cyclicAMIPatch_.
                    nbrPatch().AMI().nonOverlapTargetFaces();

                forAll(ntgtFaces, fI)
                {
                    vBC[ntgtFaces[fI]] = defValue;
                }
            }
        }
    }
    return tvBC;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::gradientInternalCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    tmp<Field<Type>> tvIC =
        this->coupledFvPatchField<Type>::gradientInternalCoeffs(deltaCoeffs);
    Field<Type>& vIC = tvIC.ref();

    if
    (
        cyclicAMIPatch_.applyLowWeightCorrection()
      && this->cyclicAMIPatch().cyclicAMIPatch().implicitNonCovered()
    )
    {
        if (implicitCoeffs())
        {
            Type defValue = pTraits<Type>::zero;
            if (cyclicAMIPatch_.cyclicAMIPatch().owner())
            {
                const labelList& nsrcFaces =
                    cyclicAMIPatch_.AMI().nonOverlapSourceFaces();
                forAll(nsrcFaces, fI)
                {
                    vIC[nsrcFaces[fI]] = defValue;
                }
            }
            else
            {
                const labelList& ntgtFaces =
                    cyclicAMIPatch_.nbrPatch().AMI().
                    nonOverlapTargetFaces();
                forAll(ntgtFaces, fI)
                {
                    vIC[ntgtFaces[fI]] = defValue;
                }
            }
        }
    }
    return tvIC;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    tmp<Field<Type>> tvBC =
        this->coupledFvPatchField<Type>::gradientBoundaryCoeffs(deltaCoeffs);
    Field<Type>& vBC = tvBC.ref();

    if
    (
        cyclicAMIPatch_.applyLowWeightCorrection()
      && this->cyclicAMIPatch().cyclicAMIPatch().implicitNonCovered()
    )
    {
        if (implicitCoeffs())
        {
            Type defValue = pTraits<Type>::zero;
            if (cyclicAMIPatch_.cyclicAMIPatch().owner())
            {
                const labelList& nsrcFaces =
                    cyclicAMIPatch_.AMI().nonOverlapSourceFaces();
                forAll(nsrcFaces, fI)
                {
                    vBC[nsrcFaces[fI]] = defValue;
                }
            }
            else
            {
                const labelList& ntgtFaces =
                    cyclicAMIPatch_.nbrPatch().AMI().
                    nonOverlapTargetFaces();
                forAll(ntgtFaces, fI)
                {
                    vBC[ntgtFaces[fI]] = defValue;
                }
            }
        }
    }
    return tvBC;
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
