/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
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

#include "mainABL.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::mainABL::kappaDefault_ = 0.41;

const Foam::scalar Foam::mainABL::betaStarDefault_ = 0.09;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mainABL::init()
{
    if (mag(flowDir_) < small || mag(zDir_) < small)
    {
        FatalErrorInFunction
            << "magnitude of n or z must be greater than zero"
            << abort(FatalError);
    }

    // Ensure direction vectors are normalized
    flowDir_ /= mag(flowDir_);
    zDir_ /= mag(zDir_);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mainABL::mainABL()
:
    flowDir_(Zero),
    zDir_(Zero),
    kappa_(0.41),
    betaStar_(0.09),
    Uref_(0),
    Zref_(0),
    A_(0),
    B_(0),
    C_(0),
    D_(0.1),
    E_(0),
    betaStarMax_(0.15),
    z0_(0),
    ustar_(0),
    BW_(0),
    hABL_(100)
{}


Foam::mainABL::mainABL
(
    const vector& flowDir,
    const vector& zDir,
    const scalar Uref,
    const scalar Zref,
    const scalar A,
    const scalar B,
    const scalar C,
    const scalar D,
    const scalar E,
    const scalar betaStarMax,
    const scalar z0,
    const scalar ustar,
    const scalar BW,
    const scalar hABL,
    const scalar kappa,
    const scalar betaStar
)
:
    flowDir_(flowDir),
    zDir_(zDir),
    kappa_(kappa),
    betaStar_(betaStar),
    Uref_(Uref),
    Zref_(Zref),
    A_(A),
    B_(B),
    C_(C),
    D_(D),
    E_(E),
    betaStarMax_(betaStarMax),
    z0_(z0),
    ustar_(ustar),
    BW_(BW),
    hABL_(hABL)

{
    init();
}


Foam::mainABL::mainABL
(
    const vectorField& p,
    const dictionary& dict
)
:
    flowDir_(dict.lookup("flowDir")),
    zDir_(dict.lookup("zDir")),
    kappa_(dict.lookupOrDefault<scalar>("kappa", kappaDefault_)),
    betaStar_(dict.lookupOrDefault<scalar>("betaStar", betaStarDefault_)),
    Uref_(readScalar(dict.lookup("Uref"))),
    Zref_(readScalar(dict.lookup("Zref"))),
    A_(readScalar(dict.lookup("A"))),
    B_(readScalar(dict.lookup("B"))),
    C_(readScalar(dict.lookup("C"))),
    D_(readScalar(dict.lookup("D"))),
    E_(readScalar(dict.lookup("E"))),
    betaStarMax_(readScalar(dict.lookup("betaStarMax"))),
    z0_(readScalar(dict.lookup("z0"))),
    ustar_(readScalar(dict.lookup("ustar"))),
    BW_(readScalar(dict.lookup("BW"))),
    hABL_(readScalar(dict.lookup("hABL")))
{
    init();
}


Foam::mainABL::mainABL
(
    const mainABL& abl,
    const fvPatchFieldMapper& mapper
)
:
    flowDir_(abl.flowDir_),
    zDir_(abl.zDir_),
    kappa_(abl.kappa_),
    betaStar_(abl.betaStar_),
    Uref_(abl.Uref_),
    Zref_(abl.Zref_),
    A_(abl.A_),
    B_(abl.B_),
    C_(abl.C_),
    D_(abl.D_),
    E_(abl.E_),
    betaStarMax_(abl.betaStarMax_),
    z0_(abl.z0_),
    ustar_(abl.ustar_),
    BW_(abl.BW_),
    hABL_(abl.hABL_)
{}


Foam::mainABL::mainABL(const mainABL& abl)
:
    flowDir_(abl.flowDir_),
    zDir_(abl.zDir_),
    kappa_(abl.kappa_),
    betaStar_(abl.betaStar_),
    Uref_(abl.Uref_),
    Zref_(abl.Zref_),
    A_(abl.A_),
    B_(abl.B_),
    C_(abl.C_),
    D_(abl.D_),
    E_(abl.E_),
    betaStarMax_(abl.betaStarMax_),
    z0_(abl.z0_),
    ustar_(abl.ustar_),
    BW_(abl.BW_),
    hABL_(abl.hABL_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mainABL::autoMap(const fvPatchFieldMapper& m)
{
    //m(z0_, z0_);
    //m(Ustar_, Ustar_);
}


void Foam::mainABL::rmap
(
    const mainABL& blptf,
    const labelList& addr
)
{
    //z0_.rmap(blptf.z0_, addr);
    //Ustar_.rmap(blptf.Ustar_, addr);
}


Foam::tmp<Foam::vectorField> Foam::mainABL::U
(
    const vectorField& p
) const
{
    const scalarField Un
    (
        (ustar_/kappa_)
       *log(((zDir_ & p) + z0_)/z0_)
    );
        return flowDir_*Un;

}


Foam::tmp<Foam::scalarField> Foam::mainABL::k
(
    const vectorField& p
) const
{
    if(BW_ == 0)
    {tmp<scalarField> tk
    (
     A_*log(((zDir_ & p) + z0_)/z0_) + B_*sqr(((zDir_ & p) + z0_)/z0_) + C_*(((zDir_ & p) + z0_)/z0_)+ D_ +E_*log((zDir_ & p) + z0_)
    );
        return tk;
    }
    else
    {
    tmp<scalarField> tk
    (
    (pow(ustar_,2)/2)*(8.7-6*(zDir_ & p)/hABL_)
    );
        return tk;
    }

}


Foam::tmp<Foam::scalarField> Foam::mainABL::omega
(
    const vectorField& p
) const
{
 
    scalarField arg = (zDir_ & p) + z0_;
    scalarField arg2 = ((zDir_ & p) + z0_)/z0_;
    scalarField kref_ = A_*log(arg2) + B_*sqr(arg2) + C_*(arg2)+ D_ +E_*log(arg);
    
    scalarField hombetaStar = min(pow(ustar_,4)/sqr(kref_),betaStarMax_);

    
    tmp<scalarField> tomega
    (
        (ustar_/sqrt(hombetaStar))/(kappa_*arg)
    );

    return tomega;
}


void Foam::mainABL::write(Ostream& os) const
{
    writeEntry(os, "flowDir", flowDir_);
    writeEntry(os, "zDir", zDir_);
    writeEntry(os, "kappa", kappa_);
    writeEntry(os, "betaStar", betaStar_);
    writeEntry(os, "Uref", Uref_);
    writeEntry(os, "Zref", Zref_);
    writeEntry(os, "A", A_);
    writeEntry(os, "B", B_);
    writeEntry(os, "C", C_);
    writeEntry(os, "D", D_);
    writeEntry(os, "E", E_);
    writeEntry(os, "betaStarMax", betaStarMax_);
    writeEntry(os, "z0", z0_) ;
    writeEntry(os, "ustar", ustar_) ;
    writeEntry(os, "BW", BW_) ;
    writeEntry(os, "hABL", hABL_) ;
}


// ************************************************************************* //
