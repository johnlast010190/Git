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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "include/swakCloudTypes.H"

#ifndef FOAM_BARYCENTRIC_PARTICLES

#include "CorrectParticleCell.H"
#include "db/IOstreams/Pstreams/Pstream.H"
#include "containers/Lists/ListListOps/ListListOps.H"
#include "db/IOobjects/IOPtrList/IOPtrList.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CorrectParticleCell<CloudType>::CorrectParticleCell
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    out_(
        this->outputDir(),
        this->owner().time()
    ),
    search_(
        this->owner().mesh()
    ),
    logCorrected_(
        dict.lookupOrDefault<bool>("logCorrected",false)
    )
{
    out_.addSpec(
        "correctedCell.*",
        "corrected outside"
    );
    out_.addSpec(
        "logCorrectedParticles.*",
        CloudType::parcelType::propertyList()
    );
    out_.addSpec(
        "logOutsideParticles.*",
        CloudType::parcelType::propertyList()
    );
}


template<class CloudType>
Foam::CorrectParticleCell<CloudType>::CorrectParticleCell
(
    const CorrectParticleCell<CloudType>& ppm
)
:
    CloudFunctionObject<CloudType>(ppm),
    out_(
        this->outputDir(),
        this->owner().time()
    ),
    search_(
        this->owner().mesh()
    ),
    logCorrected_(
        ppm.logCorrected_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CorrectParticleCell<CloudType>::~CorrectParticleCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CorrectParticleCell<CloudType>::write()
{
}

template<class CloudType>
void Foam::CorrectParticleCell<CloudType>::preEvolve()
{
    if (
        this->owner().mesh().changing()
    ) {
        Info<< this->modelName() << ":" << this->owner().name()
            << ":" << this->modelType()
            << ": Mesh moving" << endl;

        search_.correct();
    }

    label cnt=0;
    label outCnt=0;

    forAllIter(typename CloudType,this->owner(),iter) {
        parcelType &p=iter();
        label oldCellI=p.cell();
        label cellI=search_.findCell(
            p.position(),
            oldCellI
        );
        if (
            cellI<0
            // ||
            // (cellI % 4)==0
        ) {
            cnt++;
            if (logCorrected_) {
                out_("logOutsideParticles")
                    << p << endl;
            }
            //            Info<< "Not in Mesh" << endl;
            // label tetC=-1,tetP=-1,newCell=-1;
            // this->owner().mesh().findCellFacePt(
            //     p.position(),
            //     newCell,
            //     tetC,
            //     tetP
            // );
            // Info<< p.position() << " " << newCell << " " << tetC << " " << tetP << endl;
            // Info<< "Old: " << oldCellI << " " << p.tetFace() << " " << p.tetPt() << endl;
        } else if (cellI!=oldCellI) {
            if (logCorrected_) {
                out_("logCorrectedParticles")
                    << p << endl;
            }
            // Info<< "Cell: " << cellI << " old: " << oldCellI << endl;
            label tetC=-1,tetP=-1,newCell=-1;
            this->owner().mesh().findCellFacePt(
                p.position(),
                newCell,
                tetC,
                tetP
            );
            // Info<< "Corrected: " << p.position() << " " << newCell << " " << tetC << " " << tetP << endl;
            // Info<< "Old: " << p.cell() << " "<< p.tetFace() << " " << p.tetPt() << endl;

            outCnt++;
            p.cell()=newCell;
            p.tetFace()=tetC;
            p.tetPt()=tetP;
            p.cell()=cellI;
            p.initCellFacePt();
        }
    }
    if (outCnt>0) {
        Pout<< outCnt << " particles not in the right cell" << endl;
    }

    if (Pstream::parRun()) {
        out_["correctedCellProc"+Foam::name(Pstream::myProcNo())]
            << outCnt << tab << cnt << endl;
    }
    reduce(
        std::tie(cnt, outCnt),
        UniformParallelOp<plusOp<label>, 2>{}
    );
    if (Pstream::master()) {
        out_["correctedCellTotal"]
            << outCnt << tab << cnt << endl;
    }
    if (outCnt>0) {
        Info<< this->modelName() << ":" << this->owner().name()
            << ":" << this->modelType()
            << "Corrected " << outCnt << " particles" << endl;
    }
    if (cnt>0) {
        Info<< this->modelName() << ":" << this->owner().name()
            << ":" << this->modelType()
            << "Not in mesh " << cnt << " particles" << endl;
    }
}

#endif

// ************************************************************************* //
