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
    (c) 2019-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "helyxMap.H"
#include "global/argList/argList.H"
#include "db/IOobjectList/IOobjectList.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::helyxMap::trimLeft(word& str, const word& del)
{
    label i0 = str.find_first_not_of(del);
    if (i0 > 0)
    {
        label sz = str.size() - i0;
        str = str.substr(i0, sz);
    }

    return i0;
}


Foam::label Foam::helyxMap::trimRight(word& str, const word& del)
{
    label i0 = str.find_last_not_of(del);
    label sz = str.size() - 1;
    if (i0 > 0 && i0 < sz)
    {
        str = str.substr(0, i0 + 1);
    }

    return i0;
}


void Foam::helyxMap::trim(word& str, const word& del)
{
    trimLeft(str, del);
    trimRight(str, del);
}


void Foam::helyxMap::replaceString(word& s1)
{
    word nstr;
    word s2;
    char c0 = ' ';

    // Remove \t, \r, \n
    bool tostart = false;
    label i2 = s1.find_last_not_of(' ');

    for (label i = 0; i <= i2; i++)
    {
        char c = s1[i];
        if (!tostart)
        {
            if (c != ' ' && c != '\r' && c != '\t' && c != '\n')
            {
                tostart = true;

                if (c0 == ',' && c == ',')
                {
                    s2.push_back('-');
                }

                c0 = c;
                s2.push_back(c);
            }
        }
        else
        {
            if (c == '\t' || c == '\r' || c == '\n')
            {
                continue;
            }

            if (c == ' ' && c0 == ' ')
            {
                continue;
            }

            if (c0 == ',' && c == ',')
            {
                s2.push_back('-');
            }

            s2.push_back(c);
            c0 = c;
        }
    }

    s1 = s2;
}


void Foam::helyxMap::stringTokenize
(
    wordList& container,
    const word& in0,
    const word& del
)
{
    word in = in0;
    if (del == ",")
    {
        replaceString(in);
    }

    label len = in.length();
    label i = 0;

    if (container.size() >= 1)
    {
        container.resize(0);
    }

    while (i < len)
    {
        i = in.find_first_not_of(del, i);
        if (i < 0)
        {
            return;
        }
        label j = in.find_first_of(del, i);

        if (j < 0)
        {
            word ss = in.substr(i);
            trim(ss, " ");
            container.append(ss);

            break;
        }
        else
        {
            word ss1 = in.substr(i, j - i);
            trim(ss1, " ");
            container.append(ss1);
        }

        i = j + 1;
    }

    for (label k = 0; k < container.size(); k++)
    {
        word sk = container[k];
        trim(sk, " ");
        container[k] = sk;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::helyxMap::helyxMap()
:
    norm_(false),
    parModel_(2),    // memory-efficient mode
    extCoef_(0.02),
    bboxType_("smartBoundBox"),    // or 'boundBox'
    boxCellNum_(40000),
    reportMem_(false),
    scaleSource_(false),
    caseType_("vehicle"),
    nProcs_(1),
    myId_(0),
    masterId_(0),
    deltaBox_(0.1),
    wDistMap_(false),
    alphaMax_(0.5),
    errorBound_(1.0e-7),
    rhoRefSource_(1.205),
    rhoRefTarget_(1.205),
    UrotDegreeFromSource_(0),
    interp_(false),
    mapBoundary_(true),
    mapFixedBCs_(wordList(0)),
    function_("mapping"),    // or 'validate'
    mapTimeName_("latestTime"),
    tgtTimeName_("latestTime"),
    sourceCase_(word::null),
    nwdist_(100),
    sourceRegions_(wordList(1, fvMesh::defaultRegion)),
    targetRegions_(wordList(1, fvMesh::defaultRegion)),
    targetRegionDirs_(wordList(1, word::null)),
    allRegions_(false)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::helyxMap::~helyxMap()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::helyxMap::mapError(const word& name)
{
    const word unknown("unknown");
    const word& type = fieldTypes_.lookup(name, unknown);

    if (type == "unknown")
    {
        Info<< "Error validating field mapping for " << name
            << " unknown field." << endl;

        return 0;
    }

    scalar perror = 0;
    scalar small = 1.0e-25;
    const fvMesh* mesh = target().mesh_;
    scalar coef = target().scale(name);

    forAll(mesh->C(), i)
    {
        scalar x = mesh->C()[i].x();
        scalar y = mesh->C()[i].y();
        scalar z = mesh->C()[i].z();

        x = xbar(x);
        y = ybar(y);
        z = zbar(z);

        scalar alphai = target().wDists_[i];
        std::vector<label> nb;
        std::vector<scalar>dist;
        source().searchNearest(nb, dist, alphai, x, y, z);

        if (nb.empty())
        {
            perror += 1.0;    // totally missed
            continue;
        }

        label k = nb[0];
        label ic = source().cellIndex(alphai, k);

        scalar value = source().scalarFields_[name][ic];
        scalar pi = value*coef;
        scalar prefi = target().getScalarField(name, i);

        scalar perrori = fabs(pi - prefi)/(fabs(prefi) + small);

        target().getScalarField(name).primitiveFieldRef()[i] = pi;
        perror += perrori;
    }

    Info<< "Total relative mapping error for scalar field: " << name
        << " is " << perror << endl;

    return perror;
}


Foam::scalar Foam::helyxMap::UmapError(const word& name)
{
    const word unknown("unknown");
    const word& type = fieldTypes_.lookup(name, unknown);

    if (type == "unknown")
    {
        Info<< "Error validating field mapping for " << name
            << " unknown field." << endl;

        return 0;
    }

    scalar pError = 0;
    scalar small = 1.0e-25;
    const fvMesh* mesh = target().mesh_;
    scalar coef = target().scale(name);

    forAll(mesh->C(), i)
    {
        scalar x = mesh->C()[i].x();
        scalar y = mesh->C()[i].y();
        scalar z = mesh->C()[i].z();

        x = xbar(x);
        y = ybar(y);
        z = zbar(z);

        scalar alphai = target().wDists_[i];
        std::vector<label> nb;
        std::vector<scalar>dist;
        source().searchNearest(nb, dist, alphai, x, y, z);

        if (nb.empty())
        {
            pError += 1.0;    // totally missed
            continue;
        }

        label k = nb[0];
        label ic = source().cellIndex(alphai, k);

        point value = source().vectorFields_[name][ic];
        point pi = value*coef;
        point prefi = target().getVectorField(name, i);

        scalar pErrori = mag(pi - prefi)/(mag(prefi) + small);
        target().getVectorField(name).primitiveFieldRef()[i] = pi;
        pError += pErrori;
    }

    Info<< "Total relative mapping error for vector field: " << name
        << " is " << pError << endl;

    return pError;
}


Foam::label Foam::helyxMap::searchError()
{
    // Search error from target field
    label cnt = 0;
    const fvMesh* mesh = target().mesh_;
    label  ncells = mesh->C().size();
    List<label> ref(ncells);

    for (label i = 0; i < ncells; i++)
    {
        ref[i] = i;
    }

    forAll(mesh->C(), i)
    {
        scalar x = mesh->C()[i].x();
        scalar y = mesh->C()[i].y();
        scalar z = mesh->C()[i].z();

        x = xbar(x);
        y = ybar(y);
        z = zbar(z);
        scalar alphai = target().wDists_[i];

        std::vector<label> nb;
        std::vector<scalar> dist;

        source().searchNearest(nb, dist, alphai, x, y, z);

        if (nb.empty())
        {
            cnt++;
            continue;
        }

        label k = nb[0];
        label ic = source().cellIndex(alphai, k);

        if (ic != ref[i])
        {
            cnt++;
        }
    }

    Info<< "Number of mis-matched cells: " << cnt << endl;

    return cnt;
}


void Foam::helyxMap::mapFields()
{
    // Scalar fields
    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& scName = iter();

        if (parRun())
        {
            Info<< "Mapping scalar field in parallel: " << scName << endl;
            mapScalarFieldPar(scName);
        }
        else
        {
            mapScalarField(scName);
        }
    }

    wait();

    // Vector fields
    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& vecName = iter();

        if (parRun())
        {
            Info<< "Mapping vector field in parallel: " << vecName << endl;
            mapVectorFieldPar(vecName);
        }
        else
        {
            mapVectorField(vecName);
        }
    }
}


void Foam::helyxMap::mapScalarFieldPar(const word& scName)
{
    const fvMesh* mesh = target().mesh_;

    scalar coef = target().scale(scName);

    // Internal field map
    forAll(mesh->C(), i)
    {
        scalar x = mesh->C()[i].x();
        scalar y = mesh->C()[i].y();
        scalar z = mesh->C()[i].z();

        x = xbar(x);
        y = ybar(y);
        z = zbar(z);

        scalar alphai = wDistMap_ ? target().wDists_[i] : GREAT;

        std::vector<label> nb;
        std::vector<scalar> dist;

        searchNearest(nb, dist, alphai, x, y, z);

        if (nb.empty())
        {
            continue;
        }

        label k = nb[0];

        if (k < 0 || k >= sourceScalarFields_[scName].size())
        {
            continue;
        }

        // Properly considered alpha > alphaMax case
        label ic = cellIndex(alphai, k);

        if (ic < 0 || ic >= sourceScalarFields_[scName].size())
        {
            continue;
        }

        scalar value = sourceScalarFields_[scName][ic];

        scalar pi = value*coef;
        target().getScalarField(scName).primitiveFieldRef()[i] = pi;
    }

    if (mapBoundary_)
    {
        // Boundary field
        Info<< "Mapping boundary field" << endl;

        forAll(mesh->boundaryMesh(), r)
        {
            word type = mesh->boundaryMesh()[r].type();
            if
            (
                type == "empty"
             || type == "symmetryPlane"
             || type == "symmetry"
            )
            {
                continue;
            }

            if
            (
                mapFixedBCs_.find(mesh->boundaryMesh()[r].name()) == -1
             && isA<fixedValueFvPatchScalarField>
                (
                    target().getScalarField(scName).boundaryField()[r]
                )
            )
            {
                Info<< "Not mapping " << scName << " on patch "
                    << mesh->boundaryMesh()[r].name()
                    << " with boundary condition "
                    << target().getScalarField(scName).boundaryField()[r].type()
                    << endl;

                continue;
            }

            forAll(mesh->boundaryMesh()[r], j)
            {
                point ptf = mesh->Cf().boundaryField()[r][j];

                scalar x = ptf.x();
                scalar y = ptf.y();
                scalar z = ptf.z();
                x = xbar(x);
                y = ybar(y);
                z = zbar(z);

                std::vector<label> nb;
                std::vector<scalar> dist;

                searchNearest(nb, dist, x, y, z);

                if (nb.empty())
                {
                    continue;
                }

                label k = nb[0];
                label ic = cellMap_[k];

                scalar pi = sourceScalarFields_[scName][ic];
                pi *= coef;
                target().getScalarField(scName).boundaryFieldRef()[r][j] = pi;
            }
        }
    }

    Info<< "Scalar field " << scName << " mapped." << nl << endl;
}


void Foam::helyxMap::getFieldTypes(const fvMesh& mesh, const word& timeName)
{
    IOobjectList objects(mesh, timeName);

    IOobjectList fields  = objects.lookupClass(volScalarField::typeName);
    IOobjectList vfields = objects.lookupClass(volVectorField::typeName);

    mapScalarFields_ = HashSet<word>();
    mapVectorFields_ = HashSet<word>();

    if (mapFieldNames_.size() == 0)
    {
        forAllIter(IOobjectList, fields, fieldIter)
        {
            mapScalarFields_.insert(fieldIter()->name());
        }
        forAllIter(IOobjectList, vfields, fieldIter)
        {
            mapVectorFields_.insert(fieldIter()->name());
        }
    }
    else
    {
        forAll(mapFieldNames_, i)
        {
            const word& name = mapFieldNames_[i];
            if (fields.lookup(name) != nullptr)
            {
                mapScalarFields_.insert(name);
            }
            else if (vfields.lookup(name) != nullptr)
            {
                mapVectorFields_.insert(name);
            }
        }
    }
}


void Foam::helyxMap::scalePressure
(
    const dimensionSet& sourceDim,
    const dimensionSet& targetDim
)
{
    // Unset parallel flag to avoid problems with rho field missing
    // on some processors
    const bool oldParRun = Pstream::parRun();
    Pstream::parRun() = false;

    if
    (
        sourceDim == dimensionSet(1, -1, -2, 0, 0, 0, 0)
     && targetDim == dimensionSet(0,  2, -2, 0, 0, 0, 0)
    )
    {
        Info<< "Dividing source pressure by rho field..." << endl;

        volScalarField rhoField
        (
            IOobject
            (
                "rho",
                source().runTime_->timeName(),
                *(source().mesh_),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            *(source().mesh_),
            dimensionedScalar("one", dimMass/dimVolume, 1)
        );

        Info<< "rho dimensions: " << rhoField.dimensions() << endl;

        for (label i = 0; i < rhoField.primitiveField().size(); i++)
        {
            source().scalarFields_["p"][i] /= rhoField.primitiveField()[i];
        }

        // Boundaries
        if (mapBoundary_)
        {
            label cnt = rhoField.primitiveField().size();

            forAll(source().mesh_->boundaryMesh(), r)
            {
                word type = source().mesh_->boundaryMesh()[r].type();
                if (type == "empty" || type == "symmetryPlane")
                {
                    continue;
                }

                forAll(source().mesh_->boundaryMesh()[r], j)
                {
                    scalar sc = rhoField.boundaryField()[r][j];
                    source().scalarFields_["p"][cnt] /= sc;
                    cnt++;
                }
            }
        }
    }
    else if
    (
        targetDim == dimensionSet(1, -1, -2, 0, 0, 0, 0)
     && sourceDim == dimensionSet(0,  2, -2, 0, 0, 0, 0)
    )
    {
        Info<< "Multiplying source pressure by rho (rhoRefSource = "
            << rhoRefSource_ << ")" << endl;

        for (label i = 0; i < source().scalarFields_["p"].size(); i++)
        {
            source().scalarFields_["p"][i] *= rhoRefSource_;
        }
    }

    Pstream::parRun() = oldParRun;
    Info<< endl;
}


void Foam::helyxMap::checkFieldDimensions()
{
    Info<< "Checking field dimensions..." << endl;

    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& name = iter();

        const volScalarField& sourceField = source().getScalarField(name);
        const dimensionSet& sourceDim = sourceField.dimensions();

        const volScalarField& targetField = target().getScalarField(name);
        const dimensionSet& targetDim = targetField.dimensions();

        if (sourceDim != targetDim)
        {
            if (name == "p")
            {
                Info<< "Source field " << name << ", dimensions: "
                    << sourceDim << endl;
                Info<< "Target field " << name << ", dimensions: "
                    << targetDim << endl;

                scalePressure(sourceDim, targetDim);
            }
            else
            {
                FatalErrorInFunction
                    << "Field " << name << " dimensions do not match." << nl
                    << "Source dimensions: " << sourceDim << nl
                    << "Target dimensions: " << targetDim << nl
                    << exit(FatalError);
            }
        }
    }

    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& name = iter();

        const volVectorField& sourceField = source().getVectorField(name);
        const dimensionSet& sourceDim = sourceField.dimensions();

        const volVectorField& targetField = target().getVectorField(name);
        const dimensionSet& targetDim = targetField.dimensions();

        if (sourceDim != targetDim)
        {
            FatalErrorInFunction
                << "Field " << name << " dimensions do not match." << nl
                << "Source dimensions: " << sourceDim << nl
                << "Target dimensions: " << targetDim << nl
                << exit(FatalError);
        }
    }
}


void Foam::helyxMap::mapScalarField(const word& scName)
{
    const fvMesh* mesh = target().mesh_;

    scalar coef = target().scale(scName);

    // Internal field map
    label cnt = 0;
    forAll(mesh->C(), i)
    {
        scalar x = mesh->C()[i].x();
        scalar y = mesh->C()[i].y();
        scalar z = mesh->C()[i].z();
        x = xbar(x);
        y = ybar(y);
        z = zbar(z);

        scalar alphai = wDistMap_ ? target().wDists_[i] : GREAT;

        std::vector<label> nb;
        std::vector<scalar> dist;
        nb.resize(2);
        dist.resize(2);

        source().searchNearest(nb, dist, alphai, x, y, z);

        if (nb.empty())
        {
            cnt += 1;
            continue;
        }

        label k = nb[0];

        label ic = source().cellIndex(alphai, k);

        if (ic != i)
        {
            cnt++;
        }

        scalar value = source().scalarFields_[scName][ic];
        scalar pi = value*coef;
        target().getScalarField(scName).primitiveFieldRef()[i] = pi;
    }

    if (mapBoundary_)
    {
        // Boundary field
        Info<< "Mapping boundary field..." << endl;

        forAll(mesh->boundaryMesh(), r)
        {
            word type = mesh->boundaryMesh()[r].type();
            if
            (
                type == "empty"
             || type == "symmetryPlane"
             || type == "symmetry"
            )
            {
                continue;
            }

            if
            (
                mapFixedBCs_.find(mesh->boundaryMesh()[r].name()) == -1
             && isA<fixedValueFvPatchScalarField>
                (
                    target().getScalarField(scName).boundaryField()[r]
                )
            )
            {
                Info<< "Not mapping " << scName << " on patch "
                    << mesh->boundaryMesh()[r].name()
                    << " with boundary condition "
                    << target().getScalarField(scName).boundaryField()[r].type()
                    << endl;

                continue;
            }

            forAll(mesh->boundaryMesh()[r], j)
            {
                point ptf = mesh->Cf().boundaryField()[r][j];
                scalar x = ptf.x();
                scalar y = ptf.y();
                scalar z = ptf.z();
                x = xbar(x);
                y = ybar(y);
                z = zbar(z);

                std::vector<label> nb;
                std::vector<scalar> dist;

                source().searchNearestBFace(nb, dist, x, y, z);

                if (nb.empty())
                {
                    continue;
                }

                label k = nb[0];
                label ic = source().bfaceIndex(k);

                scalar pi = source().scalarFields_[scName][ic];
                pi *= coef;
                target().getScalarField(scName).boundaryFieldRef()[r][j] = pi;
            }
        }
    }

    Info<< "Scalar field " << scName << " mapped." << " " << cnt << nl << endl;
}


void Foam::helyxMap::mapVectorFieldPar(const word& vecName)
{
    const fvMesh* mesh = target().mesh_;

    scalar coef = target().scale(vecName);

    // Internal field map
    forAll(mesh->C(), i)
    {
        scalar x = mesh->C()[i].x();
        scalar y = mesh->C()[i].y();
        scalar z = mesh->C()[i].z();
        x = xbar(x);
        y = ybar(y);
        z = zbar(z);

        scalar alphai = wDistMap_ ? target().wDists_[i] : GREAT;

        std::vector<label> nb;
        std::vector<scalar> dist;
        searchNearest(nb, dist, alphai, x, y, z);

        if (nb.empty())
        {
            continue;
        }

        label k = nb[0];

        if (k < 0 || k >= sourceVectorFields_[vecName].size())
        {
            continue;
        }

        label ic = cellIndex(alphai, k);
        point value = sourceVectorFields_[vecName][ic];

        point pi = value*coef;
        target().getVectorField(vecName).primitiveFieldRef()[i] = pi;
    }

    if (mapBoundary_)
    {
        // Boundary field
        Info<< "Mapping boundary field..." << endl;

        forAll(mesh->boundaryMesh(), r)
        {
            word type = mesh->boundaryMesh()[r].type();
            if
            (
                type == "empty"
             || type == "symmetryPlane"
             || type == "symmetry"
            )
            {
                continue;
            }

            if
            (
                mapFixedBCs_.find(mesh->boundaryMesh()[r].name()) == -1
             && isA<fixedValueFvPatchVectorField>
                (
                    target().getVectorField(vecName).boundaryField()[r]
                )
            )
            {
                Info<< "Not mapping " << vecName << " on patch "
                    << mesh->boundaryMesh()[r].name()
                    << " with boundary condition "
                    << target().getVectorField(vecName).boundaryField()[r].type()
                    << endl;

                continue;
            }

            forAll(mesh->boundaryMesh()[r], j)
            {
                point ptf = mesh->Cf().boundaryField()[r][j];
                scalar x = ptf.x();
                scalar y = ptf.y();
                scalar z = ptf.z();
                x = xbar(x);
                y = ybar(y);
                z = zbar(z);

                std::vector<label> nb;
                std::vector<scalar> dist;

                searchNearest(nb, dist, x, y, z);

                if (nb.empty())
                {
                    continue;
                }

                label k = nb[0];
                label ic = cellMap_[k];
                point pi = sourceVectorFields_[vecName][ic];
                pi *= coef;
                target().getVectorField(vecName).boundaryFieldRef()[r][j] = pi;
            }
        }
    }

    wait();

    Info<< "Vector field " << vecName << " mapped." << nl << endl;
}


void Foam::helyxMap::mapVectorField(const word& vecName)
{
    const fvMesh* mesh = target().mesh_;

    scalar coef = target().scale(vecName);
    scalar beta = UrotDegreeFromSource_*3.1415926/180.0;

    // Internal field map
    forAll(mesh->C(), i)
    {
        scalar x = mesh->C()[i].x();
        scalar y = mesh->C()[i].y();
        scalar z = mesh->C()[i].z();
        x = xbar(x);
        y = ybar(y);
        z = zbar(z);

        scalar alphai = wDistMap_ ? target().wDists_[i] : GREAT;

        std::vector<label> nb;
        std::vector<scalar> dist;

        source().searchNearest(nb, dist, alphai, x, y, z);

        if (nb.empty())
        {
            continue;
        }

        label k = nb[0];
        label ic = source().cellIndex(alphai, k);

        point value = source().vectorFields_[vecName][ic];
        scalar pmag = mag(value);

        if (fabs(beta) > 1.0e-4 && pmag > 1.0e-3)
        {
            scalar alpha0 = asin(value.y()/pmag);
            scalar vx = pmag*cos(alpha0 - beta);
            scalar vy = pmag*sin(alpha0 - beta);
            scalar vz = value.z();
            value = point(vx, vy, vz);
        }
        point pi = value*coef;

        target().getVectorField(vecName).primitiveFieldRef()[i] = pi;
    }

    if (mapBoundary_)
    {
        // Boundary field
        forAll(mesh->boundaryMesh(), r)
        {
            word type = mesh->boundaryMesh()[r].type();
            if
            (
                type == "empty"
             || type == "symmetryPlane"
             || type == "symmetry"
            )
            {
                continue;
            }

            if
            (
                mapFixedBCs_.find(mesh->boundaryMesh()[r].name()) == -1
             && isA<fixedValueFvPatchVectorField>
                (
                    target().getVectorField(vecName).boundaryField()[r]
                )
            )
            {
                Info<< "Not mapping " << vecName << " on patch "
                    << mesh->boundaryMesh()[r].name()
                    << " with boundary condition "
                    << target().getVectorField(vecName).boundaryField()[r].type()
                    << endl;

                continue;
            }

            forAll(mesh->boundaryMesh()[r], j)
            {
                point ptf = mesh->Cf().boundaryField()[r][j];
                scalar x = ptf.x();
                scalar y = ptf.y();
                scalar z = ptf.z();
                x = xbar(x);
                y = ybar(y);
                z = zbar(z);

                std::vector<label> nb;
                std::vector<scalar> dist;

                source().searchNearestBFace(nb, dist, x, y, z);

                if (nb.empty())
                {
                    continue;
                }

                label k = nb[0];
                label ic = source().bfaceIndex(k);

                point pi = source().vectorFields_[vecName][ic];
                scalar pmag = mag(pi);

                if (fabs(beta) > 1.0e-4 && pmag > 1.0e-3)
                {
                    scalar alpha0 = asin(pi.y()/pmag);
                    scalar vx = pmag*cos(alpha0 - beta);
                    scalar vy = pmag*sin(alpha0 - beta);
                    scalar vz = pi.z();
                    pi = point(vx, vy, vz);
                }

                pi *= coef;
                target().getVectorField(vecName).boundaryFieldRef()[r][j] = pi;
            }
        }
    }

    Info<< "Vector field " << vecName << " mapped." << nl << endl;
}


void Foam::helyxMap::setInput(const dictionary& dict)
{
    if (sourceCase_ == word::null)
    {
        // helyxSample initializes sourceCase to "./"
        sourceCase_ = fileName(dict.lookup("sourceCase"));
    }

    if (dict.found("fields"))
    {
        mapFieldNames_ = wordList(dict.lookup("fields"));
    }

    if (dict.found("mapScalarFields"))
    {
        mapScalarFields_ = HashSet<word>(dict.lookup("mapScalarFields"));
    }

    if (dict.found("mapVectorFields"))
    {
        mapVectorFields_ = HashSet<word>(dict.lookup("mapVectorFields"));
    }

    nwdist_ = dict.lookupOrDefault<label>("nwdist", 80);

    rhoRefSource_ = dict.lookupOrDefault<scalar>("rhoRefSource", 1.205);

    rhoRefTarget_ = dict.lookupOrDefault<scalar>("rhoRefTarget", 1.205);

    vector uref(30, 0, 0);
    UrefSource_ = dict.lookupOrDefault<vector>("UrefSource", uref);
    UrefTarget_ = dict.lookupOrDefault<vector>("UrefTarget", uref);

    interp_ = dict.lookupOrDefault<bool>("interpolation", false);

    UrotDegreeFromSource_ =
        dict.lookupOrDefault<scalar>("UrotDegreeFromSource", 0);

    alphaMax_ = dict.lookupOrDefault<scalar>("alphaMax", 0.5);

    errorBound_ = dict.lookupOrDefault<scalar>("errorBound", 1.0e-7);

    if (dict.found("sourceTime"))
    {
        ITstream sourceTimeEntry = dict.lookup("sourceTime");

        if (sourceTimeEntry[0] == word("latestTime"))
        {
            mapTimeName_ = word("latestTime");
        }
        else
        {
            scalar sourceTimeScalar = dict.lookup<scalar>("sourceTime");
            mapTimeName_ = name(sourceTimeScalar);
        }
    }

    if (dict.found("targetTime"))
    {
        ITstream targetTimeEntry = dict.lookup("targetTime");

        if (targetTimeEntry[0] == word("latestTime"))
        {
            tgtTimeName_ = word("latestTime");
        }
        else
        {
            scalar targetTimeScalar = dict.lookup<scalar>("targetTime");
            tgtTimeName_ = name(targetTimeScalar);
        }
    }

    if (dict.found("sourceRegion"))
    {
        sourceRegions_ = wordList(1, word(dict.lookup("sourceRegion")));
    }

    if (dict.found("targetRegion"))
    {
        targetRegions_ = wordList(1, word(dict.lookup("targetRegion")));
    }

    if (dict.found("allRegions"))
    {
        allRegions_ = dict.lookup<bool>("allRegions");
    }

    if (dict.found("regionMaps"))
    {
        List<Pair<word>> regionMap;
        regionMap = List<Pair<word>>(dict.lookup("regionMaps"));

        sourceRegions_.resize(0);
        targetRegions_.resize(0);

        forAll(regionMap, regioni)
        {
            sourceRegions_.append(regionMap[regioni].first());
            targetRegions_.append(regionMap[regioni].second());
        }

        targetRegionDirs_ = targetRegions_;
    }

    mapBoundary_ = dict.lookupOrDefault<bool>("mapBoundary", true);

    mapFixedBCs_ = dict.lookupOrDefault<wordList>("mapFixedBCs", wordList(0));
    if (dict.found("mapFixedBC"))
    {
        WarningInFunction
            << "mapFixedBC has been deprecated. "
            << "Please, use mapFixedBCs (List<word>) instead!" << endl;
    }

    caseType_ = dict.lookupOrDefault<word>("caseType", "vehicle");
    wDistMap_ = dict.lookupOrDefault<bool>("wdist_map", false);
    extCoef_ = dict.lookupOrDefault<scalar>("extensionCoef", 0.1);
    parModel_ = dict.lookupOrDefault<label>("parallelModel", 2);
    function_ = dict.lookupOrDefault<word>("function", "mapping");
    norm_ = dict.lookupOrDefault<bool>("coordsNormalise", false);

    // Memory usage info
    reportMem_ = dict.lookupOrDefault<bool>("reportMemoryUsage", false);

    // Debug flag
    bool debug = dict.lookupOrDefault<bool>("debug", false);
    if (debug)
    {
        reportMem_ = true;
    }

    if (dict.found("fieldTypes"))
    {
        fieldTypes_ = HashTable<word>(dict.lookup("fieldTypes"));
    }

    // Wall patch ID's for calculating wall distances
    if (dict.found("patches"))
    {
        wallPatchNames_ = wordReList(dict.lookup("patches"));
    }

    // Smart boundBox for HPC run
    bboxType_ = dict.lookupOrDefault<word>("boundBoxType", "smartBoundBox");

    // Backward-compatibility for legacy bound box keyword
    if (bboxType_ == "smBoundBox")
    {
        DeprecationWarningInFunction
        (
            "smBoundBox",
            "bound box type name",
            40400,
            "Please use the 'smartBoundBox' type name instead."
        );

        bboxType_ = "smartBoundBox";
    }

    boxCellNum_ = dict.lookupOrDefault<label>("gridCells", 1000000);

    scaleSource_ = dict.lookupOrDefault<bool>("scaleSource", false);
}


void Foam::helyxMap::setInput(const dictionary& dict, const argList& args)
{
    if (sourceCase_ == word::null && !args.optionFound("sourceCase"))
    {
        // helyxSample initializes sourceCase to "./"
        sourceCase_ = fileName(dict.lookup("sourceCase"));
    }

    if (dict.found("fields"))
    {
        mapFieldNames_ = wordList(dict.lookup("fields"));
    }

    if (dict.found("mapScalarFields"))
    {
        mapScalarFields_ = HashSet<word>(dict.lookup("mapScalarFields"));
    }

    if (dict.found("mapVectorFields"))
    {
        mapVectorFields_ = HashSet<word>(dict.lookup("mapVectorFields"));
    }

    nwdist_ = dict.lookupOrDefault<label>("nwdist", 80);

    rhoRefSource_ = dict.lookupOrDefault<scalar>("rhoRefSource", 1.205);

    rhoRefTarget_ = dict.lookupOrDefault<scalar>("rhoRefTarget", 1.205);

    vector uref(30, 0, 0);
    UrefSource_ = dict.lookupOrDefault<vector>("UrefSource", uref);
    UrefTarget_ = dict.lookupOrDefault<vector>("UrefTarget", uref);

    interp_ = dict.lookupOrDefault<bool>("interpolation", false);

    UrotDegreeFromSource_ =
        dict.lookupOrDefault<scalar>("UrotDegreeFromSource", 0);

    alphaMax_ = dict.lookupOrDefault<scalar>("alphaMax", 0.5);

    errorBound_ = dict.lookupOrDefault<scalar>("errorBound", 1.0e-7);

    if (dict.found("sourceTime"))
    {
        ITstream sourceTimeEntry = dict.lookup("sourceTime");

        if (sourceTimeEntry[0] == word("latestTime"))
        {
            mapTimeName_ = word("latestTime");
        }
        else
        {
            scalar sourceTimeScalar = dict.lookup<scalar>("sourceTime");
            mapTimeName_ = name(sourceTimeScalar);
        }
    }

    if (dict.found("targetTime"))
    {
        ITstream targetTimeEntry = dict.lookup("targetTime");

        if (targetTimeEntry[0] == word("latestTime"))
        {
            tgtTimeName_ = word("latestTime");
        }
        else
        {
            scalar targetTimeScalar = dict.lookup<scalar>("targetTime");
            tgtTimeName_ = name(targetTimeScalar);
        }
    }

    if (dict.found("sourceRegion"))
    {
        sourceRegions_ = wordList(1, word(dict.lookup("sourceRegion")));
    }

    if (dict.found("targetRegion"))
    {
        targetRegions_ = wordList(1, word(dict.lookup("targetRegion")));
    }

    if (dict.found("allRegions"))
    {
        allRegions_ = dict.lookup<bool>("allRegions");
    }

    if (dict.found("regionMaps"))
    {
        List<Pair<word>> regionMap;
        regionMap = List<Pair<word>>(dict.lookup("regionMaps"));

        sourceRegions_.resize(0);
        targetRegions_.resize(0);

        forAll(regionMap, regioni)
        {
            sourceRegions_.append(regionMap[regioni].first());
            targetRegions_.append(regionMap[regioni].second());
        }

        targetRegionDirs_ = targetRegions_;
    }

    mapBoundary_ = dict.lookupOrDefault<bool>("mapBoundary", true);

    mapFixedBCs_ = dict.lookupOrDefault<wordList>("mapFixedBCs", wordList(0));
    if (dict.found("mapFixedBC"))
    {
        WarningInFunction
            << "mapFixedBC has been deprecated. "
            << "Please, use mapFixedBCs (List<word>) instead!" << endl;
    }

    caseType_ = dict.lookupOrDefault<word>("caseType", "vehicle");
    wDistMap_ = dict.lookupOrDefault<bool>("wdist_map", false);
    extCoef_ = dict.lookupOrDefault<scalar>("extensionCoef", 0.1);
    parModel_ = dict.lookupOrDefault<label>("parallelModel", 2);
    function_ = dict.lookupOrDefault<word>("function", "mapping");
    norm_ = dict.lookupOrDefault<bool>("coordsNormalise", false);

    // Memory usage info
    reportMem_ = dict.lookupOrDefault<bool>("reportMemoryUsage", false);

    // Debug
    bool debug = dict.lookupOrDefault<bool>("debug", false);
    if (debug)
    {
        reportMem_ = true;
    }

    if (dict.found("fieldTypes"))
    {
        fieldTypes_ = HashTable<word>(dict.lookup("fieldTypes"));
    }

    // Wall patch ID's for calculating wall distances
    if (dict.found("patches"))
    {
        wallPatchNames_ = wordReList(dict.lookup("patches"));
    }

    // Smart boundBox for HPC run
    bboxType_ = dict.lookupOrDefault<word>("boundBoxType", "smartBoundBox");

    // Backward-compatibility for legacy bound box keyword
    if (bboxType_ == "smBoundBox")
    {
        DeprecationWarningInFunction
        (
            "smBoundBox",
            "bound box type name",
            40400,
            "Please use the 'smartBoundBox' type name instead."
        );

        bboxType_ = "smartBoundBox";
    }

    boxCellNum_ = dict.lookupOrDefault<label>("gridCells", 1000000);

    scaleSource_ = dict.lookupOrDefault<bool> ("scaleSource", false);
}


void Foam::helyxMap::getWallDistPatchs(const fvMesh& mesh)
{
    Info<< "\nCalculating wall distance for ";

    // Wall patch ID's for calculating wall distances
    if (wallPatchNames_.size() == 0)
    {
        Info<< "all patches." <<endl;
        return;
    }

    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    wallDistPatchs_ = bMesh.patchSet(wallPatchNames_);

    const wordList allPatchNames(bMesh.names());

    Info<< wallDistPatchs_.size() << " patch"
        << (wallDistPatchs_.size() > 1 ? "es:" : ":")
        << endl;

    forAllConstIter(labelHashSet, wallDistPatchs_, it)
    {
        Info<< token::TAB << allPatchNames[*it] <<endl;
    }

    Info<< nl;
}


void Foam::helyxMap::setOptions(const argList& args)
{
    args.optionReadIfPresent("sourceCase", sourceCase_);

    if (args.optionFound("mapScalarFields"))
    {
        args.optionLookup("mapScalarFields")() >> mapScalarFields_;
    }

    if (args.optionFound("mapVectorFields"))
    {
        args.optionLookup("mapVectorFields")() >> mapVectorFields_;
    }

    if (args.optionFound("mapSurfaceScalarFields"))
    {
        args.optionLookup("mapSurfaceScalarFields")()
            >> mapSurfaceScalarFields_;
    }

    args.optionReadIfPresent("alphaMax", alphaMax_);

    args.optionReadIfPresent("errorBound", errorBound_);

    if (args.optionFound("mapBoundary"))
    {
        mapBoundary_ = true;
    }

    if (args.optionFound("mapFixedBCs"))
    {
        mapFixedBCs_ = wordList(1, args["mapFixedBCs"]);

        mapFixedBCs_.resize(0);
        List<word> fixedBCs;
        args.optionLookup("mapFixedBCs")() >> fixedBCs;

        forAll(fixedBCs, bci)
        {
            mapFixedBCs_.append(fixedBCs[bci]);
        }
    }

    args.optionReadIfPresent("sourceTime", mapTimeName_);

    args.optionReadIfPresent("targetTime", tgtTimeName_);

    if (args.optionFound("sourceRegion"))
    {
        sourceRegions_ = wordList(1, args["sourceRegion"]);
    }

    if (args.optionFound("targetRegion"))
    {
        targetRegions_ = wordList(1, args["targetRegion"]);
        targetRegionDirs_ = targetRegions_;
    }

    if (args.optionFound("allRegions"))
    {
        allRegions_ = true;
    }

    if (args.optionFound("regionMaps"))
    {
        sourceRegions_.resize(0);
        targetRegions_.resize(0);

        List<Pair<word>> regionMap;
        args.optionLookup("regionMaps")() >> regionMap;

        forAll(regionMap, regioni)
        {
            sourceRegions_.append(regionMap[regioni].first());
            targetRegions_.append(regionMap[regioni].second());
        }

        targetRegionDirs_ = targetRegions_;
    }

    args.optionReadIfPresent("function", function_);

    label nndist = 0;
    args.optionReadIfPresent("nwdist", nndist);
    if (nndist > 2 && nndist <= 1001)
    {
        nwdist_ = nndist;
    }

    args.optionReadIfPresent("rhoRefSource", rhoRefSource_);

    args.optionReadIfPresent("rhoRefTarget", rhoRefTarget_);

    args.optionReadIfPresent("UrefSource", UrefSource_);

    args.optionReadIfPresent("UrefTarget", UrefTarget_);

    args.optionReadIfPresent("UrotDegreeFromSource", UrotDegreeFromSource_);

    if (args.optionFound("interpolation"))
    {
        interp_ = true;
    }

    args.optionReadIfPresent("fieldTypes", fieldTypes_);
}


void Foam::helyxMap::setMapTime(Time& runTime)
{
    scalar itime = atof(runTime.timeName().c_str());
    scalar maptime = atof(mapTimeName_.c_str());

    runTime += (maptime - itime);
}


void Foam::helyxMap::buildKdTreesPar()
{
    Info<< "Building parallel search trees..." << endl;

    sourceAlphaMax_ = alphaMax_;

    label treeNum = nwdist_ + 1;
    deltaY_ = 1.0/float(nwdist_);

    kdTrees_.resize(treeNum);
    alphaMap_.resize(treeNum);

    const point& pmin = target().gBox().min();
    const point& pmax = target().gBox().max();

    xmin_ = pmin.x();
    xmax_ = pmax.x();
    ymin_ = pmin.y();
    ymax_ = pmax.y();
    zmin_ = pmin.z();
    zmax_ = pmax.z();

    // Put the cell into groups
    Array2d<point> pxyz;
    pxyz.resize(treeNum);

    label n2 = sourceWDists_.size();

    forAll(sourceXyz_, i)
    {
        const point& pt = sourceXyz_[i];
        if (i > n2 - 1 || sourceWDists_[i] > sourceAlphaMax_) continue;

        label ic = inGroup(i);
        if (ic < 0) ic = 0;
        if (ic > treeNum - 1) ic = treeNum - 1;

        scalar xb = xbar(pt.x());
        scalar yb = ybar(pt.y());
        scalar zb = zbar(pt.z());
        point pt1 (xb, yb, zb);
        pxyz[ic].push_back(pt1);
        alphaMap_[ic].push_back(i);
    }

    label inactiveKdTrees = 0;
    if (wDistMap_)
    {
        // Use relative coordinates to build kdTrees
        for (unsigned ic = 0; ic < pxyz.size(); ic++)
        {
            label npt = pxyz[ic].size();
            if (npt < 1)
            {
                kdTrees_[ic].deactivate();
                inactiveKdTrees++;

                continue;
            }

            farray2d trainData(npt);
            for (label i = 0; i < npt; i++)
            {
                trainData[i].resize(3);

                point pt = pxyz[ic][i];
                trainData[i][0] = pt.x();
                trainData[i][1] = pt.y();
                trainData[i][2] = pt.z();
            }

            if (!interp_)
            {
                kdTrees_[ic].k = 1;
            }

            kdTrees_[ic].eps = errorBound_;

            kdTrees_[ic].build(trainData);
        }
    }

    if (wDistMap_ && inactiveKdTrees < label(kdTrees_.size()))
    {
        scalar coef = 0.8;

        constructInternalKnn(coef*alphaMax_);
    }
    else
    {
        wDistMap_ = false;

        if (!interp_)
        {
            internalMap_.k = 1;
        }

        constructInternalKnn(-0.1);
    }

    wait();

    Info<< "Build of parallel search trees completed." << endl;
}


void Foam::helyxMap::createSource(const fvMesh* mesh, const Time* runTime)
{
    source_ = new readFields(mesh, runTime, this, "source");

    source().setNorm(norm_);
    source().constructBoundBox(0.0, bboxType_, boxCellNum_);
}


void Foam::helyxMap::createTarget(const fvMesh* mesh, const Time* runTime)
{
    target_ = new readFields(mesh, runTime, this, "target");

    target().setNorm(norm_);
    target().constructBoundBox(extCoef_, bboxType_, boxCellNum_);
}


void Foam::helyxMap::clearSourceTarget()
{
    internalMap_.deAllocate();

    if (wDistMap_)
    {
        for (size_t i = 0; i < kdTrees_.size(); i++)
        {
            if (kdTrees_[i].active)
            {
                kdTrees_[i].deAllocate();
            }
        }
    }

    sourceXyz_.clear();
    sourceScalarFields_.clear();
    sourceVectorFields_.clear();
    sourceWDists_.clear();

    source_.clear();
    target_.clear();
}


void Foam::helyxMap::constructInternalKnn(scalar cutoff)
{
    label cnt = 0;

    if (!wDistMap_)
    {
        cutoff = -0.1;
    }

    farray2d trainData;

    // For gridMap (cutoff < 0), no need to do alpha-mapping.
    // Internal map contains all cells.
    label ncells = sourceXyz_.size();
    if (wDistMap_)
    {
        ncells = min(sourceXyz_.size(), sourceWDists_.size());
    }

    for (label i = 0; i < ncells; i++)
    {
        if (wDistMap_)
        {
            scalar alphai = sourceWDists_[i];
            if (alphai < cutoff)
            {
                continue;
            }
        }

        point pt = sourceXyz_[i];
        scalar xb = xbar(pt.x());
        scalar yb = ybar(pt.y());
        scalar zb = zbar(pt.z());
        std::vector<scalar> pts(3);
        pts[0] = xb;
        pts[1] = yb;
        pts[2] = zb;

        trainData.push_back(pts);

        cellMap_[cnt] = i;
        cnt++;
    }

    internalMap_.eps = errorBound_;

    internalMap_.build(trainData);
}


void Foam::helyxMap::parInit(label np)
{
    // Initialization for parallel run
    if (!parRun())
    {
        return;
    }

    procs_.resize(np);
    for (label i = 0; i < np; i++)
    {
        procs_[i] = i;
    }
}


bool Foam::helyxMap::inList(label k, const DynamicList<label>& list)
{
    if (list.size() == 0)
    {
        return false;
    }

    for (label i = 0; i < list.size(); i++)
    {
        if (list[i] == k)
        {
            return true;
        }
    }

    return false;
}


bool Foam::helyxMap::inList(const word& key, const wordList& list)
{
    if (list.size() == 0)
    {
        return false;
    }

    for (label i = 0; i < list.size(); i++)
    {
        if (list[i] == key)
        {
            return true;
        }
    }

    return false;
}


template<class T>
void Foam::interProcTrans
(
    DynamicList<T>& field,
    const DynamicList<point>& verts,
    const List<boundBox>& tboxs,
    const List<boundBox>& sboxs,
    const DynamicList<T>& source,
    const word& name
)
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    label myId = Pstream::myProcNo();

    for (label ib = 0; ib < Pstream::nProcs(); ib++)
    {
        const boundBox& tgtbox = tboxs[ib];
        const boundBox& mySourceBox = sboxs[myId];

        // Skip if mySourceBox does not overlap target box
        DynamicList<T> cellField;
        bool boxInDomain = mySourceBox.overlaps(tgtbox);
        if (boxInDomain)
        {
            for (label i = 0; i < source.size(); i++)
            {
                const point& pt = verts[i];

                if (tgtbox.contains(pt))
                {
                    cellField.append(source[i]);
                }
            }
        }

        UOPstream toProc(ib, pBufs);
        toProc << cellField;
    }

    pBufs.finishedSends();

    field.resize(0);
    for (label ib = 0; ib < Pstream::nProcs(); ib++)
    {
        DynamicList<T> posGlob;

        UIPstream fromProc(ib, pBufs);
        fromProc >> posGlob ;
        for (label i = 0; i < posGlob.size(); i++)
        {
            field.append(posGlob[i]);
        }
    }

    pBufs.finishedSends();

    Info<< "Field transferred for: " << name << ", size: "
        << field.size() << endl;
}


void Foam::helyxMap::reportMemory
(
    memInfo& mem,
    const word& desc,
    scalar& maxUsage
)
{
    scalar memory = scalar(mem.update().size())*1.0e-3;

    if (memory > maxUsage)
    {
        scalar maxUsage = memory;

        Pout<< "Maximum memory usage now: " << maxUsage
            << " MB, " << desc << endl;
    }
}


void Foam::helyxMap::getDomainFields()
{
    if (parModel() == 1)
    {
        return getDomainFields0();
    }

    const smartBoundBox& tBox = target().bBox();

    List<boundBox> tBoxs(Pstream::nProcs());
    label myId = Pstream::myProcNo();
    tBoxs[myId] = tBox;

    Pstream::allGatherList(tBoxs);

    const smartBoundBox& sBox = source().bBox();

    List<boundBox> sBoxs(Pstream::nProcs());
    sBoxs[myId] = sBox;
    Pstream::allGatherList(sBoxs);

    scalar maxUsage = 0;
    word desc = "getDomainFields: before get sourceXYZ";
    memInfo mem;

    if (reportMemUsage())
    {
        reportMemory(mem, desc, maxUsage);
    }

    DynamicList<point> sourceXyz;

    interProcTrans
    (
        sourceXyz,
        source().xyz(),
        tBoxs,
        sBoxs,
        source().xyz(),
        "cellVerts"
    );

    List<bool> status(sourceXyz.size(), false);

    sourceXyz_.resize(0);

    label cellCount = 0;
    forAll(sourceXyz, i)
    {
        const point& pt = sourceXyz[i];
        if (tBox.pointInBox(pt))
        {
            status[i] = true;
            sourceXyz_.append(pt);
            cellCount++;
        }
    }

    if (reportMemUsage())
    {
        Pout<< "Number of cells in this processor: "
            << source().xyz().size() << endl;

        Pout<< "Points from all processors: " << sourceXyz.size() << endl;
    }

    sourceXyz.clear();

    if (reportMemUsage())
    {
        desc = "getDomainFields: after get sourceXYZ";

        reportMemory(mem, desc, maxUsage);

        Pout<< "sourceXyz_: " << sourceXyz_.size() << endl;
        Pout<< "Number of active sourceXyz_ count: " << cellCount << endl;
    }

    wait();

    if (wDistMap_)
    {
        DynamicList<scalar> sourceWdists;

        // Wall distance
        interProcTrans
        (
            sourceWdists,
            source().xyz(),
            tBoxs,
            sBoxs,
            source().wDists_,
            "wallDistance"
        );

        forAll(sourceWdists, i)
        {
            if (i < status.size() && status[i])
            {
                sourceWDists_.append(sourceWdists[i]);
            }
        }

        sourceWdists.clear();

        if (reportMemUsage())
        {
            desc = "getDomainFields: after get wallDists";

            reportMemory(mem, desc, maxUsage);
        }

        wait();
    }

    // Scalar fields
    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& name = iter();

        const DynamicList<scalar>& scfld = source().scalarFields_[name];

        DynamicList<scalar> sourceScalarFields;
        interProcTrans
        (
            sourceScalarFields,
            source().xyz(),
            tBoxs,
            sBoxs,
            scfld,
            name
        );

        sourceScalarFields_[name].resize(0);
        forAll(sourceScalarFields, i)
        {
            if (i < status.size() && status[i])
            {
                sourceScalarFields_[name].append(sourceScalarFields[i]);
            }
        }

        sourceScalarFields.clear();

        if (reportMemUsage())
        {
            desc = "getDomainFields: after get scalar field: ";
            desc += name;

            reportMemory(mem, desc, maxUsage);
        }
    }

    wait();

    // Vector fields
    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& name = iter();

        const DynamicList<vector>& vecfld = source().vectorFields_[name];

        DynamicList<vector> sourceVectorFields;
        interProcTrans
        (
            sourceVectorFields,
            source().xyz(),
            tBoxs,
            sBoxs,
            vecfld,
            name
        );

        sourceVectorFields_[name].resize(0);
        forAll(sourceVectorFields, i)
        {
            if (i < status.size() && status[i])
            {
                sourceVectorFields_[name].append(sourceVectorFields[i]);
            }
        }

        sourceVectorFields.clear();
    }

    if (reportMemUsage())
    {
        desc = "getDomainFields: after get vector fields";

        reportMemory(mem, desc, maxUsage);
    }

    wait();
}


void Foam::helyxMap::getDomainFields0()
{
    if (!Pstream::parRun())
    {
        return;
    }

    Info<< "Single-block parallel model..." << endl;

    List<vectorField> posLoc(Pstream::nProcs());
    label myId = Pstream::myProcNo();
    posLoc[myId] = source().xyz();

    Pstream::allGatherList(posLoc);

    vectorField posGlob =
        ListListOps::combine<vectorField>(posLoc, accessOp<vectorField>());

    cellInDomain_.resize(posGlob.size());
    cellInDomain_ = false;

    sourceXyz_.resize(0);
    forAll(posGlob, i)
    {
        point pt = posGlob[i];
        bool indm = false;

        if (target().inDomain(pt))
        {
            sourceXyz_.append(pt);
            indm = true;
        }

        cellInDomain_[i] = indm;
    }

    posGlob.resize(0);
    posLoc.clear();

    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& name = iter();

        const DynamicList<scalar>& fld = source().scalarFields_[name];
        label sz = fld.size();

        scalarField scfldLoc(sz);
        for (label i = 0; i < sz; i++)
        {
            scfldLoc[i] = fld[i];
        }

        List<scalarField> sField(Pstream::nProcs());
        sField[myId] = scfldLoc;
        Pstream::allGatherList(sField);

        scalarField fldGlobal =
            ListListOps::combine<scalarField>(sField, accessOp<scalarField>());

        label sz1 = fldGlobal.size();
        if (sz1 > cellInDomain_.size())
        {
            sz1 = cellInDomain_.size();
        }

        sourceScalarFields_[name].resize(0);

        for (label i = 0; i < sz1; i++)
        {
            if (cellInDomain_[i])
            {
                sourceScalarFields_[name].append(fldGlobal[i]);
            }
        }

        wait();

        fldGlobal.clear();
        sField.clear();
    }

    // Wall distance

    label sz = source().wDists_.size();

    scalarField scfldLoc(sz);
    for (label i = 0; i < sz; i++)
    {
        // Recover absolute value
        scfldLoc[i] = source().wDists_[i]*source().maxWallDist_;
    }

    List<scalarField> wDists(Pstream::nProcs());
    wDists[myId] = scfldLoc;
    Pstream::allGatherList(wDists);

    scalarField fldGlobal =
        ListListOps::combine<scalarField>(wDists, accessOp<scalarField>());
    scalar maxwdist = max(fldGlobal);

    fldGlobal /= maxwdist;

    label sz2 = fldGlobal.size();
    if (sz2 > cellInDomain_.size())
    {
        sz2 = cellInDomain_.size();
    }

    sourceWDists_.resize(0);
    for (label i = 0; i < sz2; i++)
    {
        if (cellInDomain_[i])
        {
            sourceWDists_.append(fldGlobal[i]);
        }
    }

    wait();

    fldGlobal.clear();
    wDists.clear();

    // Velocity fields
    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& name = iter();

        const DynamicList<point>& fld = source().vectorFields_[name];
        label sz = fld.size();

        vectorField vectLoc(sz);
        for (label i = 0; i < sz; i++)
        {
            vectLoc[i] = fld[i];
        }

        List<vectorField> vField(Pstream::nProcs());
        vField[myId] = vectLoc;
        Pstream::allGatherList(vField);

        vectorField fldGlobal =
            ListListOps::combine<vectorField>(vField, accessOp<vectorField>());

        label sz1 = fldGlobal.size();
        if (sz1 > cellInDomain_.size())
        {
            sz1 = cellInDomain_.size();
        }

        sourceVectorFields_[name].resize(0);

        for (label i = 0; i < sz1; i++)
        {
            if (cellInDomain_[i])
            {
                sourceVectorFields_[name].append(fldGlobal[i]);
            }
        }

        wait();

        vField.clear();
        fldGlobal.clear();
    }
}


// ************************************************************************* //
