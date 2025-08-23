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
    Foam::ReaderParticle

Description

SourceFiles

Contributors/Copyright:
    2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef ReaderParticle_H
#define ReaderParticle_H

#include "include/swak.H"

#include "include/swakCloudTypes.H"


#ifdef FOAM_PARTICLE_CLASS_IS_TEMPLATE
#include "Particle.H"
#else
#include "particle/particle.H"
#endif

#include "db/IOstreams/IOstreams/IOstream.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ReaderParticle Declaration
\*---------------------------------------------------------------------------*/

    class ReaderParticleCloud;

#ifndef FOAM_PARTICLE_CLASS_IS_TEMPLATE
    typedef particle ReaderParticleBaseClass;
#else
    class ReaderParticle;

    typedef Particle<ReaderParticle> ReaderParticleBaseClass;
#endif

class ReaderParticle
:
    public ReaderParticleBaseClass
{

public:

#ifndef FOAM_PARTICLE_CLASS_IS_TEMPLATE
    TypeName("ReaderParticle");
#endif

    // Constructors

#ifndef FOAM_PARTICLE_CLASS_IS_TEMPLATE
        //- Construct from components
        ReaderParticle
        (
            const polyMesh& mesh,
#ifdef FOAM_BARYCENTRIC_PARTICLES
            const barycentric &position,
#else
            const vector& position,
#endif
            const label celli,
            const label tetFaceI,
            const label tetPtI

        )
        :
            ReaderParticleBaseClass(mesh, position, celli, tetFaceI, tetPtI)
        {}

        ReaderParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            bool doCellFacePt = true
        )
        :
            ReaderParticleBaseClass(
                mesh, position, cellI
#ifndef FOAM_BARYCENTRIC_PARTICLES
                , doCellFacePt
#endif
            )
        {}

        //- Construct from Istream
        ReaderParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        )
        :
            ReaderParticleBaseClass(mesh, is, readFields)
        {}
#else
        //- Construct from components
        ReaderParticle
        (
            const Cloud<ReaderParticle>& c,
            const vector& position,
            const label celli
        )
        :
            Particle<ReaderParticle>(c, position, celli)
        {}

        //- Construct from Istream
        ReaderParticle
        (
            const Cloud<ReaderParticle>& c,
            Istream& is,
            bool readFields = true
        )
        :
            Particle<ReaderParticle>(c, is, readFields)
        {}
#endif

        //- Construct as copy
        ReaderParticle(const ReaderParticle& p)
        :
            ReaderParticleBaseClass(p)
        {}

        //- Construct and return a clone
        virtual autoPtr<ReaderParticleBaseClass> clone() const
        {
            return autoPtr<ReaderParticleBaseClass>(new ReaderParticle(*this));
        }

        static void readFields (ReaderParticleCloud &c);

        static void writeFields (const Cloud<ReaderParticle> &c);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
