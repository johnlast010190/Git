/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH

Class
    Foam::groovyBCCommon

Description
    Foam::groovyBCCommon

SourceFiles
    groovyBCCommon.C

Contributors/Copyright:
    2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef groovyBCCommon_H
#define groovyBCCommon_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "db/dictionary/dictionary.H"
#include "helpers/exprString.H"

namespace Foam
{

    class pointPatch;
    class fvPatch;

/*---------------------------------------------------------------------------*\
                     Class groovyBCFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class groovyBCCommon
{

    bool evaluateDuringConstruction_;

protected:

    // Protected data

    bool debug_;
    bool hasGradient_;

    // the expressions
    exprString valueExpression_;
    exprString gradientExpression_;
    exprString fractionExpression_;

    bool evaluateDuringConstruction() const
        { return evaluateDuringConstruction_; }

public:

    static exprString nullValue();

    static const fvPatch &getFvPatch(const pointPatch &pp);

    // Constructors

    //- Construct from nothing
    groovyBCCommon
    (
        bool hasGradient,
        bool isPoint=false,
        exprString fractionExpression="1"
    );

    //- Construct from dictionary
    groovyBCCommon
    (
        const dictionary&,
        bool hasGradient,
        bool isPoint=false,
        exprString fractionExpression="1"
    );

    //- Copy constructor
    groovyBCCommon
    (
        const groovyBCCommon<Type> &
    );

    // Member functions

    //- Write
    void write(Ostream&) const;
};


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #ifdef NoRepository
    #   include "groovyBCCommon.C"
    #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #endif

// ************************************************************************* //
