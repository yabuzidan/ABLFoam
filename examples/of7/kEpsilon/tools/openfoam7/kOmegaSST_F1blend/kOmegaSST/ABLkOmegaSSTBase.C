/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "ABLkOmegaSSTBase.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //


template<class TurbulenceModel, class BasicTurbulenceModel>
void ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::referencesComputation()
{
    bool debug=1;
    if (debug)
        Info << "Entered to referencesComputation" << endl;
    
    // Set reference values for internal fields
    volScalarField coord = (centres_ + z0_)/z0_;
    volScalarField coord2 = (centres_ + z0_);
    dimensionSet coord2d(coord2.dimensions());
    coord2.dimensions().reset(dimless);
    
    uref_ = ustar_*log(coord)/kappa_;
    kref_ = A_*log(coord) + B_*sqr(coord) + C_*coord + D_ + E_*log(coord2);
    
    coord2.dimensions().reset(coord2d);
    oref_ = (1/(kappa_*ustar_*coord2))*kref_;
    
    hombetaStar_ = homogeneousbetaStar();

// Boundary fields
    const fvPatchList& patches = this->mesh_.boundary();
    volScalarField::Boundary& centresBf_ = this->centres_.boundaryFieldRef();
    volScalarField::Boundary& krefBf_ = this->kref_.boundaryFieldRef();
    volScalarField::Boundary& orefBf_ = this->oref_.boundaryFieldRef();
    volScalarField::Boundary& urefBf_ = this->uref_.boundaryFieldRef();
    volScalarField::Boundary& hombetaStarBf_ = this->hombetaStar_.boundaryFieldRef();


	scalar coordt=0;
    	scalar coord2t=0;
    
    forAll(patches, patchi)
    {
        fvPatchScalarField& pcentres_  = centresBf_[patchi];
        fvPatchScalarField& pkref_  = krefBf_[patchi];
        fvPatchScalarField& puref_  = urefBf_[patchi];
        fvPatchScalarField& phombetaStar_  = hombetaStarBf_[patchi];
        fvPatchScalarField& poref_  = orefBf_[patchi];


        const fvPatch& curPatch = patches[patchi];
        
        /*if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];
                wSwitch_[faceCelli] = 0.0;
            }
        }*/
        
        forAll(curPatch, facei)
        {

            coordt = (pcentres_[facei] + z0_.value())/z0_.value();
            coord2t = (pcentres_[facei] + z0_.value());
            
            puref_[facei] = ustar_.value()*log(coordt)/kappa_.value();
            pkref_[facei] = A_.value()*log(coordt) + B_.value()*sqr(coordt) + C_.value()*coordt + D_.value() +E_.value()*log(coord2t);
	    poref_[facei] = (1/(kappa_.value()*ustar_.value()*coord2t))*pkref_[facei];
	}
    }
}


template<class TurbulenceModel, class BasicTurbulenceModel>
void ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::errorsAndbetaStarComputation()
{
    bool debug=1;
    if (debug)
        Info << "Entered to computeErrorsAndbetaStar" << endl;
    Uerr_ = uError();
    Kerr_ = kError();
    Oerr_ = oError();
    Hybriderr_ = hybridError();
    BlendingFun_ = blendingFunction();
    hombetaStar_ = homogeneousbetaStar();

}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField> ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::uError()
{
    bool debug=1;
    if (debug)
        Info << "Entered to uError" << endl;
        
    volScalarField Udiff = mag((uref_ - this->U_.component(0))/max(uref_,usmall))/constantU_;
        
    forAll(Udiff,celli)
    {

        if(Udiff[celli] < thresholdU_.value())
        {
            Udiff[celli] = scalar(0);
        }
    }

    return min(Udiff,scalar(1));
}
 
template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField> ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::kError()
{
    bool debug=1;
    if (debug)
        Info << "Entered to kError" << endl;
    volScalarField Kdiff = mag((kref_ - k_)/max(kref_,ksmall))/constantK_;
        
    forAll(Kdiff,celli)
    {

        if(Kdiff[celli] < thresholdK_.value())
        {
            Kdiff[celli] = scalar(0);
        }
    }
        
    return min(Kdiff,scalar(1));
}
    
template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField> ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::oError()
{
    bool debug=1;
    if (debug)
        Info << "Entered to oError" << endl;

    volScalarField Odiff = mag((oref_ - omega_)/max(oref_,osmall))/constantO_;
        
    forAll(Odiff,celli)
    {

        if(Odiff[celli] < thresholdO_.value())
        {
            Odiff[celli] = scalar(0);
        }
    }

    return min(Odiff,scalar(1));
}
    
template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField> ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::hybridError()
{
    bool debug=1;
    if (debug)
        Info << "Entered to hybridError" << endl;

    return max(max(Uerr_,Kerr_),Oerr_);
}
    
template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField> ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::blendingFunction()
const
{
    bool debug=1;
    if (debug)
        Info << "Entered to blendingFunction" << endl;

    volScalarField marker = msel_.component(0)*Uerr_ + msel_.component(1)*Kerr_ + msel_.component(2)*Hybriderr_;

    volScalarField transErr = Foam::constant::mathematical::pi*max(marker- 0.5,-0.5);

    // Calculate and return blender value
    return pow(1.0 - 0.5*(1.0 + sin(transErr)),blendCoeff_);
}
    
// Smooth betaStar
template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField> ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::homogeneousbetaStar()
{
    bool debug=1;
    if (debug)
        Info << "Entered to homogeneousbetaStar" << endl;
    
    return min((pow(ustar_,4)/sqr(kref_)),betaStarMax_);
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField>
ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::ABLkOmegaSST::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1))*blendingFunction();
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField>
ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::ABLkOmegaSST::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2))*blendingFunction();
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField>
ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::ABLkOmegaSST::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return (1 - tanh(pow4(arg3)))*blendingFunction();
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField>
ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::ABLkOmegaSST::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class TurbulenceModel, class BasicTurbulenceModel>
void ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& F2
)
{
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F2*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class TurbulenceModel, class BasicTurbulenceModel>
void ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))), F23());
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField::Internal>
ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField::Internal>
ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::epsilonByk
(
    const volScalarField::Internal& F1,
    const volScalarField::Internal& F2
) const
{
    return betaStar_*omega_();
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<fvScalarMatrix>
ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<fvScalarMatrix>
ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<fvScalarMatrix> ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceModel, class BasicTurbulenceModel>
ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::ABLkOmegaSST
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    TurbulenceModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

    betaStarMax_
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
    

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    referencesComputation();
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);
    errorsAndbetaStarComputation();    
    correctNut();

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TurbulenceModel, class BasicTurbulenceModel>
bool ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::read()
{
    if (TurbulenceModel::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        betaStarMax_.readIfPresent(this->coeffDict());
        
        z0_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());
        ustar_.readIfPresent(this->coeffDict());
        zDir_.readIfPresent(this->coeffDict());
        flowDir_.readIfPresent(this->coeffDict());
        
        A_.readIfPresent(this->coeffDict());
        B_.readIfPresent(this->coeffDict());
        C_.readIfPresent(this->coeffDict());
        D_.readIfPresent(this->coeffDict());
        E_.readIfPresent(this->coeffDict());
        blendCoeff_.readIfPresent(this->coeffDict());
        
        thresholdU_.readIfPresent(this->coeffDict());
        thresholdK_.readIfPresent(this->coeffDict());
        thresholdO_.readIfPresent(this->coeffDict());
        constantU_.readIfPresent(this->coeffDict());
        constantK_.readIfPresent(this->coeffDict());
        constantO_.readIfPresent(this->coeffDict());
        
        msel_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class TurbulenceModel, class BasicTurbulenceModel>
void ABLkOmegaSST<TurbulenceModel, BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const surfaceScalarField& phi = this->phi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    BasicTurbulenceModel::correct();

    if (this->mesh_.changing())
    {
        referencesComputation();
    }

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu(dev(twoSymm(tgradU()())) && tgradU()());
    volScalarField::Internal G(this->GName(), nut()*GbyNu);
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    volScalarField blendNL = blendingFunction();

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));


        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(omega_)
          + fvm::div(phi, omega_)
          - fvm::laplacian(DomegaEff(F1), omega_)
         ==
           gamma
           *min
            (
                GbyNu,
                (c1_/a1_)*betaStar_*omega_()
               *max(a1_*omega_(), b1_*F23()*sqrt(S2()))
            )
          - fvm::Sp(beta*omega_(), omega_)
          - fvm::SuSp
            (
                (F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          - fvc::laplacian(DomegaEff(F1), omega_)*blendNL

          - (gamma
            *min
             (
                GbyNu,
                (c1_/a1_)*betaStar_*omega_()
               *max(a1_*omega_(), b1_*F23()*sqrt(S2()))
             ))*blendNL

          + (beta*omega_()*omega_())*blendNL

	  + ((F1()-scalar(1))*CDkOmega())*blendNL

	  //+ Qsas(S2(), gamma, beta)
          //+ omegaSource()
          //+ fvOptions(alpha, rho, omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi, k_)
      - fvm::laplacian(DkEff(F1), k_)
     ==
       Pk(G)
      - fvm::Sp(epsilonByk(F1, F23), k_)
      - fvc::laplacian(DkEff(F1), k_)*blendNL
      - Pk(G)*blendNL
      + (epsilonByk(F1, F23)*k_())*blendNL
      //+ fvc::Sp(epsilonByk(F1, F23), k_)*blendNL
      //+ kSource()
      //+ fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(S2, F23);

    errorsAndbetaStarComputation();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
