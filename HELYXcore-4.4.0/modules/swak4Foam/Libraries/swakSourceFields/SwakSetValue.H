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
    Foam::SwakSetValue

Description
    Fixes field values according to an expression

SourceFiles
    SwakSetValue.C

Contributors/Copyright:
    2010, 2012-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef SwakSetValue_H
#define SwakSetValue_H

#include "SwakBasicSourceCommon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SwakSetValue Declaration
\*---------------------------------------------------------------------------*/

#ifdef FOAM_HAS_FVOPTIONS
    namespace fv {
#endif

template<class T>
class SwakSetValue
:
    public SwakBasicSourceCommon<T>
{
    // Private data

    SwakSetValue(const SwakSetValue&);

    bool useMaskExpression_;

    exprString maskExpression_;

    bool getMask(DynamicList<label> &,const word &psi);

public:

    //- Runtime type information
    TypeName("SwakSetValue");

    // Constructors

        //- Construct from a dictionary
        SwakSetValue
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

    // Destructor

        virtual ~SwakSetValue();

        // Evaluation

            //- Set value on vector field
            virtual void setValue(fvMatrix<T>& eqn, const label fieldI);

            //- in newer versions setValue was renamed
            virtual void constrain(fvMatrix<T>& eqn, const label fieldI) {
                setValue(eqn,fieldI);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef FOAM_HAS_FVOPTIONS
    }
#endif

} // End namespace Foam

#ifdef NoRepository
#   include "SwakSetValue.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
