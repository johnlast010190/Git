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

\*---------------------------------------------------------------------------*/

#include "fields/Fields/transformField/transformField.H"
#include "fields/Fields/Field/FieldM.H"
#include "primitives/DiagTensor/diagTensor/diagTensor.H"

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

void Foam::transformPoints
(
    vectorField& rtf,
    const spatialTransform& tr,
    const vectorField& tf
)
{
    forAll(rtf, i)
    {
        rtf[i] = tr.transformPoint(tf[i]);
    }
}


Foam::tmp<Foam::vectorField> Foam::transformPoints
(
    const spatialTransform& tr,
    const vectorField& tf
)
{
    tmp<vectorField > tranf(new vectorField(tf.size()));
    transformPoints(tranf.ref(), tr, tf);
    return tranf;
}


Foam::tmp<Foam::vectorField> Foam::transformPoints
(
    const spatialTransform& tr,
    const tmp<vectorField>& ttf
)
{
    tmp<vectorField > tranf = New(ttf);
    transformPoints(tranf.ref(), tr, ttf());
    ttf.clear();
    return tranf;
}


void Foam::transform
(
    vectorField& rtf,
    const quaternion& q,
    const vectorField& tf
)
{
    tensor t = q.R();
    TFOR_ALL_F_OP_FUNC_S_F(vector, rtf, =, transform, tensor, t, vector, tf)
}


Foam::tmp<Foam::vectorField> Foam::transform
(
    const quaternion& q,
    const vectorField& tf
)
{
    tmp<vectorField > tranf(new vectorField(tf.size()));
    transform(tranf.ref(), q, tf);
    return tranf;
}


Foam::tmp<Foam::vectorField> Foam::transform
(
    const quaternion& q,
    const tmp<vectorField>& ttf
)
{
    tmp<vectorField > tranf = New(ttf);
    transform(tranf.ref(), q, ttf());
    ttf.clear();
    return tranf;
}


void Foam::transformPoints
(
    vectorField& rtf,
    const septernion& tr,
    const vectorField& tf
)
{
    vector T = tr.t();

    // Check if any translation
    if (mag(T) > VSMALL)
    {
        TFOR_ALL_F_OP_F_OP_S(vector, rtf, =, vector, tf, -, vector, T);
    }
    else
    {
        rtf = tf;
    }

    // Check if any rotation
    if (mag(tr.r().R() - I) > SMALL)
    {
        transform(rtf, tr.r(), rtf);
    }
}


Foam::tmp<Foam::vectorField> Foam::transformPoints
(
    const septernion& tr,
    const vectorField& tf
)
{
    tmp<vectorField > tranf(new vectorField(tf.size()));
    transformPoints(tranf.ref(), tr, tf);
    return tranf;
}


Foam::tmp<Foam::vectorField> Foam::transformPoints
(
    const septernion& tr,
    const tmp<vectorField>& ttf
)
{
    tmp<vectorField > tranf = New(ttf);
    transformPoints(tranf.ref(), tr, ttf());
    ttf.clear();
    return tranf;
}


// ************************************************************************* //
