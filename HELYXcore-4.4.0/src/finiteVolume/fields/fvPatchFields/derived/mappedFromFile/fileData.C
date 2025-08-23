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
    (c) 2009-2009 OpenCFD Ltd.
    (c) 2017-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/mappedFromFile/fileData.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * //

namespace Foam
{
    template<class Type>
    PtrList<Field<Type>> fileData<Type>::fields_(0);

    template<class Type>
    triSurface fileData<Type>::surface_ = triSurface();

    template<class Type>
    fileName fileData<Type>::filename_("");

    template<class Type>
    word fileData<Type>::fieldName_("");

    template<class Type>
    bool fileData<Type>::convertToMillimeters_(false);

    template<class Type>
    autoPtr<indexedOctree<treeDataPoint>> fileData<Type>::pointTreePtr_(nullptr);

    template<class Type>
    bool fileData<Type>::geometryLoaded_(false);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fileData<Type>::fileData()
:
    mapBack_()
{}


template<class Type>
Foam::fileData<Type>::fileData(const dictionary& dict)
:
    mapBack_(dict.lookupOrDefault<bool>("mapBack", false))
{
    if (filename_ == "")
    {
        filename_ = fileName(dict.lookup("file"));
    }
    if (fieldName_ == "")
    {
        fieldName_ = word(dict.lookup("field"));
    }
    convertToMillimeters_ =
        dict.lookupOrDefault<bool>("convertToMillimeters", false);
}


template<class Type>
Foam::fileData<Type>::fileData(const fileData& fd)
:
    mapBack_(fd.mapBack_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::fileData<Type>::~fileData()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::fileData<Type>::readFields()
{
    if (filename_.substr(filename_.size()-3,filename_.size()) == "nas")
    {
        readFieldsNAS();
    }
    else if (filename_.substr(filename_.size()-3,filename_.size()) == "ntl")
    {
        readFieldsNTL();
    }
}


template<typename Type>
void Foam::fileData<Type>::readFieldsNTL()
{
    // Nonsense, but preserves original insane semantics
    // without breaking the linker.
    if constexpr (!std::is_same_v<Type, scalar>)
    {
        return;
    }

    const Field<point>& pnts = surface_.points();
    PtrList<scalarField> pntFields(2);
    pntFields.set(0, new scalarField(pnts.size(), 0.0));
    pntFields.set(1, new scalarField(pnts.size(), 0.0));

    Info<< "Reading fields values from file " << filename_ << endl;
    IFstream is(filename_);

    label ctr(0);

    while (is.good())
    {
        size_t linei = 2;
        string line;
        is.getLine(line);

        // Check for a tag - space followed by a digit, or two digits,
        // followed by whitespace at the beginning of a line
        if
        (
            !(
                line.length() >= 3
             && (isdigit(line[0]) || line[0] == ' ')
             && isdigit(line[1])
             && isspace(line[2])
            )
        )
        {
            // If not a tag, just continue and read next line
            continue;
        }

        string tag = line.substr(0, 2);

        // Reading Node Temperatures
        if (tag == "10")
        {
            if (ctr == pnts.size())
            {
                // Warn and stop reading to avoid a buffer overrun
                WarningInFunction
                    << "File appears to be invalid: More nodal temperature"
                    << "values than nodes were found."
                    << endl;
                break;
            }

            is.getLine(line);
            linei = 0;

            // Get front and back temperature on nodes
            pntFields[0][ctr] = parseScalar(readToken(line, 16, linei));
            pntFields[1][ctr] = parseScalar(readToken(line, 16, linei));

            ctr++;
        }
        else if (tag == " 1")
        {
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Material Properties
        else if (tag == " 3")
        {
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Element Properties
        else if (tag == " 4")
        {
            is.getLine(line);
        }
        // Reading Coordinate Frames
        else if (tag == " 5")
        {
            is.getLine(line);
            is.getLine(line);
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Distributed Loads (Data Card number varies)
        else if (tag == " 6")
        {
            is.getLine(line);
        }
        // Reading Node Forces (Data Card number varies)
        else if (tag == " 7")
        {
            is.getLine(line);
        }
        // Reading Node Displacements
        else if (tag == " 8")
        {
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Element Temperatures
        else if (tag == "11")
        {
            is.getLine(line);
        }
        // Reading MPC Data
        else if (tag == "14")
        {
            is.getLine(line);
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Nodal Heat Source
        else if (tag == "15")
        {
            is.getLine(line);
        }
        // Reading Distributed Heat Source (Data Card number varies)
        else if (tag == "16")
        {
            is.getLine(line);
        }
        // Reading Convection Coefficients
        else if (tag == "17")
        {
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Emissivity Values
        else if (tag == "18")
        {
            is.getLine(line);
            is.getLine(line);
        }
        if (tag == "25")
        {
            is.getLine(line);
        }
        // Reading/Skipping Summary Data
        else if (tag == "26")
        {
            is.getLine(line);
        }
        else if (tag == "98")
        {
            is.getLine(line);
        }
        else if (tag == "99")
        {
            break;
        }
    }

    Info<< "Interpolating fields from face points to centres" << endl;

    fields_.setSize(2);
    fields_.set(0, new scalarField(surface_.Cf().size(), 0.0));
    fields_.set(1, new scalarField(surface_.Cf().size(), 0.0));

    forAll(surface_, fI)
    {
        forAll(surface_[fI],pI)
        {
            label pntI = surface_[fI][pI];
            fields_[0][fI] += pntFields[0][pntI]/3.;
            fields_[1][fI] += pntFields[1][pntI]/3.;
        }
    }
    Info<< "Done" << endl;
}


template<typename Type>
void Foam::fileData<Type>::readFieldsNAS()
{
    // Nonsense, but preserves original insane semantics without breaking the linker.
    if constexpr (!std::is_same_v<Type, scalar>)
    {
        return;
    }

    fields_.setSize(2);
    forAll(fields_, fldI)
    {
        const Field<point>& pnts = surface_.points();
        scalarField pntField(pnts.size(), 0.0);

        fileName fn = filename_;
        if (fldI == 1)
        {
            fn = fn.substr(0, fn.size() - 4) + "_back.nas";
        }

        Info<< "Reading fields values from file " << fn << endl;
        IFstream is(filename_);

        label ctr(0);

        while (is.good())
        {
            size_t linei = 0;
            string line;
            is.getLine(line);

            if (line.empty() || line[0] == '$')
            {
                // Skip empty or comment
                continue;
            }

            // Read first word
            word cmd(IStringStream(readToken(line, 8, linei))());

            if (cmd == fieldName_)
            {
                if (ctr == pntField.size())
                {
                    Info<< "exceeded list!!" << endl;
                }

                readToken(line, 8, linei);
                readToken(line, 8, linei);

                scalar fieldValue = parseScalar(readToken(line, 8, linei));

                pntField[ctr] = fieldValue;

                ctr++;
            }
        }

        Info<< "Interpolating fields from face points to centres" << endl;

        fields_.set(fldI, new scalarField(surface_.Cf().size(), 0.0));

        forAll(surface_, fI)
        {
            forAll(surface_[fI],pI)
            {
                label pntI = surface_[fI][pI];
                fields_[fldI][fI] += pntField[pntI]/3.;
            }
        }

        Info<< "Done" << endl;
    }
}


template<class Type>
Foam::Field<Type>* Foam::fileData<Type>::mapData
(
    const pointField& tarPoints,
    const vectorField& tarNormals
)
{
    static_assert
    (
        std::is_same_v<Type, scalar>,
        "Function only supported for scalars"
    );

    if (!geometryLoaded_)
    {
        Info<< "Reading thermal mesh" << endl;

        surface_ = triSurface(filename_);

        const vectorField& pnts = surface_.Cf();

        // Construct octree
        indexedOctree<treeDataPoint> pointTree
        (
            treeDataPoint(pnts),
            treeBoundBox(treeBoundBox(pnts).extend(1e-4)),  // overall search
            8,          // maxLevel
            10,         // leafsize
            3.0         // duplicity
        );

        pointTreePtr_ = pointTree.clone();

        geometryLoaded_ = true;
    }

    if (fields_.size() == 0)
    {
        readFields();
    }

    scalarField* mappedFieldPtr(nullptr);
    mappedFieldPtr = new scalarField(tarPoints.size(),0.);
    scalarField& mappedField = *mappedFieldPtr;

    // Loop over points, find nearest and assign mapped data
    scalar spanSqr = Foam::sqr(pointTreePtr_->bb().mag());

    forAll(tarPoints, pI)
    {
        const point& p = tarPoints[pI];

        pointIndexHit info = pointTreePtr_->findNearest(p, spanSqr);

        if (!info.hit())
        {
            info = pointTreePtr_->findNearest(p, GREAT*spanSqr);
        }

        label side = 0;

        if (mapBack_ && (tarNormals[pI] & surface_.Sf()[info.index()]) > 0)
        {
            side = 1;
        }

        mappedField[pI] = fields_[side][info.index()];
    }

    return mappedFieldPtr;
}


template<class Type>
Foam::scalar Foam::fileData<Type>::parseScalar(const string& s)
{
    size_t expSign = s.find_last_of("+-");

    if (expSign != string::npos && expSign > 0 && !isspace(s[expSign-1]))
    {
        scalar mantissa = readScalar(IStringStream(s.substr(0, expSign-1))());
        scalar exponent = readScalar(IStringStream(s.substr(expSign+1))());

        if (s[expSign] == '-')
        {
            exponent = -exponent;
        }
        return mantissa*pow(10, exponent);
    }
    else
    {
        return readScalar(IStringStream(s)());
    }
}


template<class Type>
std::string Foam::fileData<Type>::readToken
(
    const string& line,
    const size_t& width,
    size_t& index
)
{
    size_t indexStart, indexEnd;

    indexStart = index;

    indexEnd = line.find(',', indexStart);
    index = indexEnd + 1;

    if (indexEnd == std::string::npos)
    {
        indexEnd = indexStart + width;
        index = indexEnd;
    }

    return line.substr(indexStart, indexEnd - indexStart);
}

// ************************************************************************* //
