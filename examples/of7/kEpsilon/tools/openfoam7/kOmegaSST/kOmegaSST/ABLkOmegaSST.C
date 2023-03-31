/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ABLkOmegaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
ABLkOmegaSST<BasicTurbulenceModel>::ABLkOmegaSST
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    Foam::ABLkOmegaSST
    <
        eddyViscosity<RASModel<BasicTurbulenceModel>>,
        BasicTurbulenceModel
    >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )

   /* betaStarMax_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "BetaStarMax",
            this->coeffDict_,
            0.15
        )
    ),

    z0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "z0",
            this->coeffDict_,
            dimensionSet(0, 1, 0, 0, 0, 0, 0),
            0.01
        )
    ),
    
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.41
        )
    ),
    
    ustar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ustar",
            this->coeffDict_,
            dimensionSet(0, 1, -1, 0, 0, 0, 0),
            0.641
        )
    ),
    
    zDir_
    (
        dimensioned<vector>::lookupOrAddToDict
        (
            "zDir",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            vector(0, 0, 1)
        )
    ),

    flowDir_
    (
        dimensioned<vector>::lookupOrAddToDict
        (
            "flowDir",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            vector(1, 0, 0)
        )
    ),
    
    centres_
    (
        IOobject
        (
             "centres",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag(zDir_ & this->mesh_.C())

    ),

    uref_
    (
        IOobject
        (
            "uref",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("uref", dimensionSet(0,1,-1,0,0,0,0),0.0)
    ),

    kref_
    (
        IOobject
        (
            "kref",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("kref", dimensionSet(0,2,-2,0,0,0,0),0.0)
    ),

    oref_
    (
        IOobject
        (
            "oref",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("oref", dimensionSet(0,0,-1,0,0,0,0),0.0)
    ),

    Uerr_
    (
        IOobject
        (
            "Uerr",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Uerr", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    
    Kerr_
    (
        IOobject
        (
            "Kerr",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Kerr", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    
    Oerr_
    (
        IOobject
        (
            "Oerr",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Eerr", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    
    Hybriderr_
    (
        IOobject
        (
            "Hybriderr",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Hybriderr", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    
    BlendingFun_
    (
        IOobject
        (
            "BlendingFun",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("BlendingFun", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),

    A_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A",
            this->coeffDict_,
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.1
        )
    ),
    
    B_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "B",
            this->coeffDict_,
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.1
        )
    ),
    
    C_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C",
            this->coeffDict_,
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.1
        )
    ),
    
    D_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "D",
            this->coeffDict_,
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            1.3696
        )
    ),
    
    E_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "E",
            this->coeffDict_,
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.1
        )
    ),

    blendCoeff_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "blendCoeff",
            this->coeffDict_,
            0
        )
    ),
    
    thresholdU_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "thresholdU",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.05
        )
    ),
    
    thresholdK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "thresholdK",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.3
        )
    ),
    
    thresholdO_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "thresholdO",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.1
        )
    ),
    
    constantU_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "constantU",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            1
        )
    ),
    
    constantK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "constantK",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            3
        )
    ),
    
    constantO_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "constantO",
            this->coeffDict_,
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            10
        )
    ),

    msel_
     (
         dimensioned<vector>::lookupOrAddToDict
         (
             "msel",
             this->coeffDict_,
             dimensionSet(0, 0, 0, 0, 0, 0, 0),
             vector(1, 0, 0)
         )
     ),

    usmall("usmall",dimensionSet(0,1,-1,0,0,0,0),scalar(1e-10)),
    ksmall("ksmall",dimensionSet(0,2,-2,0,0,0,0),scalar(1e-10)),
    osmall("esmall",dimensionSet(0,0,-1,0,0,0,0),scalar(1e-10)),
    
    NLbetaStar_
    (
        IOobject
        (
            "NLBetaStar",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("NLBetaStar", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    
    hombetaStar_
    (
        IOobject
        (
            "homBetaStar",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("homBetaStar", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    
    blendbetaStar_
    (
        IOobject
        (
            "blendBetaStar",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("blendBetaStar", dimensionSet(0,0,0,0,0,0,0),0.0)
    ),*/

{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
