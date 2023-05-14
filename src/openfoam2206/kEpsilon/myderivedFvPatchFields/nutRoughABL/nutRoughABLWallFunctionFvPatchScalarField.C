/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "nutRoughABLWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
scalar NutRoughABLWallFunctionFvPatchScalarField::blendingFunction
(
    const scalar BIA,
    const scalar blendCoeff
) const
{
    scalar transErr = Foam::constant::mathematical::pi*max(BIA-0.5,-0.5);
        /* Calculate and return blender value */
    return pow(1.0 - 0.5*(1.0 + sin(transErr)),blendCoeff);
}
    
tmp<scalarField> NutRoughABLWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchi = patch().index();

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const tmp<volScalarField> te = turbModel.epsilon();
    const volScalarField& epsilon = te();

    tmp<scalarField> tnutw(new scalarField(*this));
    scalarField& nutw = tnutw.ref();

    // Get first cell centroid
    scalarField ccentre = zDir_ & patch().Cn();
    forAll(nutw, facei)
    {
        if(roughWall_[facei] > 0)
        {
         ccentre = y[facei];
        }
    }
    
    const scalarField coord = (ccentre + z0_)/z0_;
    const scalarField coord2 = (ccentre + z0_);
    
    // Calculate references
    const scalarField uref = ustar_*log(coord)/kappaWF_;
    const scalarField kref = A_*log(coord) + B_*sqr(coord) + C_*(coord)+ D_ +Etke_*log(coord2);
    const scalarField homCmu = min(pow(ustar_,4)/sqr(kref),CmuMax_);
    const scalarField eref = (sqrt(homCmu)*ustar_*kref)/(kappaWF_*coord2);

    // Calculate errors
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const volVectorField& Ut = turbModel.U();
    vectorField Uwt = Uw;
    forAll(nutw, facei)
    {
      label faceCelli = patch().faceCells()[facei];
      Uwt[facei] = Ut[faceCelli];
    }
    
    scalarField Uerr  = mag((uref - Uwt.component(0))/max(uref,uref*1e-14))/constantU_;
    scalarField kerr = mag((kref - k)/max(kref,kref*1e-14))/constantK_;
    scalarField Eerr = mag((eref - epsilon)/max(eref,eref*1e-14))/constantE_;
    
    forAll(Uerr, facei)
  
    {
        if(Uerr[facei] < thresholdU_[facei])
        {
            Uerr[facei] = scalar(0.0);
        }
        else if(Uerr[facei] > scalar(1))
        {
            Uerr[facei] = scalar(1.0);;
        }
    }

    forAll(kerr, facei)
  
    {
        if(kerr[facei] < thresholdK_[facei])
        {
            kerr[facei] = scalar(0);
        }
        if(kerr[facei] > scalar(1))
        {
            kerr[facei] = scalar(1);
        }
    }
    
    forAll(Eerr, facei)
  
    {
        if(Eerr[facei] < thresholdE_[facei])
        {
            Eerr[facei] = scalar(0);
        }
        if(Eerr[facei] > scalar(1))
        {
            Eerr[facei] = scalar(1);
        }
    }

    scalarField herr = max(max(Uerr,kerr),Eerr);
    
    /* Non Linear CMU computation */
    volSymmTensorField S = symm(fvc::grad(Ut));
    volTensorField W = skew(fvc::grad(Ut));

    const scalarField Xx = k/epsilon;
    const scalarField St = (Xx)*sqrt(2.0)*mag(S);
    const scalarField Ot = (Xx)*sqrt(2.0)*mag(W);

    const scalarField NLCmu = min(1/(0.9*pow(St,1.4)+0.4*pow(Ot,1.4)+3.5),0.15);
    
    forAll(nutw, facei)
    {
        label celli = patch().faceCells()[facei];
        
        scalar bCmu = NLCmu[facei] + (homCmu[facei] - NLCmu[facei])*blendingFunction(herr[facei],blendCoeff_[facei]);
        
        const scalar Cmu25 = pow025(bCmu);

        scalar uStar = Cmu25*sqrt(k[celli]);
        
        if(roughWall_[facei] > 0) // Rough case
        {
            // Original implementaiton from Leo/Parante code (not working)
            // scalar yPlus = uStar*(y[facei]+z0_[facei])/nuw[facei];
            // scalar Edash = nuw[facei]/(z0_[facei]*uStar);

            // nutw[facei] =
            // nuw[facei]*(yPlus*kappaWF_[facei]/log(Edash*yPlus) - 1);


            // YA - based on OF2206 impelmentation atmNutkWallFunction
            scalar yPlus = uStar*y[facei]/nuw[facei];
            scalar Edash = (y[facei] + z0_[facei])/z0_[facei];

            nutw[facei] = nuw[facei]*(yPlus*kappaWF_[facei]/log(max(Edash, 1 + 1e-4)) - 1);

            if (debug)
            {
                Info<< "yPlus = " << yPlus
                    << ", Edash = " << Edash
                    << ", nutw = " << nutw[facei]
                    << endl;
            }
        }
        else
        {
            scalar yPlus = uStar*y[facei]/nuw[facei];
            scalar Edash = wallCoeffs_.E();//E_; for openfoam2112
            
            // // Original implementaiton from Leo/Parante code
            // nutw[facei] =
            // nuw[facei]*(yPlus*kappaWF_[facei]/log(Edash*yPlus) - 1);
            
            // YA - based on OF2206 impelmentation nutkRoughWallFunction (improved stability)
            const scalar limitingNutw = max(nutw[facei], nuw[facei]);
            // To avoid oscillations limit the change in the wall viscosity
            // which is particularly important if it temporarily becomes zero
            nutw[facei] =
                max
                (
                    min
                    (
                        nuw[facei]
                        *(yPlus*kappaWF_[facei]/log(max(Edash*yPlus, 1 + 1e-4)) - 1),
                        2*limitingNutw
                    ), 0.5*limitingNutw
                );

            if (debug)
            {
                Info<< "yPlus = " << yPlus
                    << ", Edash = " << Edash
                    << ", nutw = " << nutw[facei]
                    << endl;
            }
        }
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

NutRoughABLWallFunctionFvPatchScalarField::
NutRoughABLWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    z0_(p.size(), 0.0),
    ustar_(p.size(), 0.0),
    kappaWF_(p.size(), 0.41),
    A_(p.size(), 0.0),
    B_(p.size(), 0.0),
    C_(p.size(), 0.0),
    D_(p.size(), 1.37),
    Etke_(p.size(), 0.0),
    CmuMax_(p.size(), 0.0),
    blendCoeff_(p.size(), 0.0),
    zDir_(p.size(), vector(0, 0, 1)),
    roughWall_(p.size(), 1.0),
    thresholdU_(p.size(), 0.3),
    thresholdE_(p.size(), 0.5),
    thresholdK_(p.size(), 0.5),
    constantU_(p.size(), 1),
    constantE_(p.size(), 3),
    constantK_(p.size(), 10)
    
{}


NutRoughABLWallFunctionFvPatchScalarField::
NutRoughABLWallFunctionFvPatchScalarField
(
    const NutRoughABLWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    z0_(mapper(ptf.z0_)),
    ustar_(mapper(ptf.ustar_)),
    kappaWF_(mapper(ptf.kappaWF_)),
    A_(mapper(ptf.A_)),
    B_(mapper(ptf.B_)),
    C_(mapper(ptf.C_)),
    D_(mapper(ptf.D_)),
    Etke_(mapper(ptf.Etke_)),
    CmuMax_(mapper(ptf.CmuMax_)),
    blendCoeff_(mapper(ptf.blendCoeff_)),
    zDir_(mapper(ptf.zDir_)),
    roughWall_(mapper(ptf.roughWall_)),
    thresholdU_(mapper(ptf.thresholdU_)),
    thresholdE_(mapper(ptf.thresholdE_)),
    thresholdK_(mapper(ptf.thresholdK_)),
    constantU_(mapper(ptf.constantU_)),
    constantE_(mapper(ptf.constantE_)),
    constantK_(mapper(ptf.constantK_))
{}


NutRoughABLWallFunctionFvPatchScalarField::
NutRoughABLWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    z0_("z0", dict, p.size()),
    ustar_("ustar", dict, p.size()),
    kappaWF_("kappaWF", dict, p.size()),
    A_("A", dict, p.size()),
    B_("B", dict, p.size()),
    C_("C", dict, p.size()),
    D_("D", dict, p.size()),
    Etke_("Etke", dict, p.size()),
    CmuMax_("CmuMax", dict, p.size()),
    blendCoeff_("blendCoeff", dict, p.size()),
    zDir_("zDir", dict, p.size()),
    roughWall_("roughWall", dict, p.size()),
    thresholdU_("thresholdU", dict, p.size()),
    thresholdE_("thresholdE", dict, p.size()),
    thresholdK_("thresholdK", dict, p.size()),
    constantU_("constantU", dict, p.size()),
    constantE_("constantE", dict, p.size()),
    constantK_("constantK", dict, p.size())
{}


NutRoughABLWallFunctionFvPatchScalarField::
NutRoughABLWallFunctionFvPatchScalarField
(
    const NutRoughABLWallFunctionFvPatchScalarField& rwfpsf
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf),
    z0_(rwfpsf.z0_),
    ustar_(rwfpsf.ustar_),
    kappaWF_(rwfpsf.kappaWF_),
    A_(rwfpsf.A_),
    B_(rwfpsf.B_),
    C_(rwfpsf.C_),
    D_(rwfpsf.D_),
    Etke_(rwfpsf.Etke_),
    CmuMax_(rwfpsf.CmuMax_),
    blendCoeff_(rwfpsf.blendCoeff_),
    zDir_(rwfpsf.zDir_),
    roughWall_(rwfpsf.roughWall_),
    thresholdU_(rwfpsf.thresholdU_),
    thresholdE_(rwfpsf.thresholdE_),
    thresholdK_(rwfpsf.thresholdK_),
    constantU_(rwfpsf.constantU_),
    constantE_(rwfpsf.constantE_),
    constantK_(rwfpsf.constantK_)
{}


NutRoughABLWallFunctionFvPatchScalarField::
NutRoughABLWallFunctionFvPatchScalarField
(
    const NutRoughABLWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf, iF),
    z0_(rwfpsf.z0_),
    ustar_(rwfpsf.ustar_),
    kappaWF_(rwfpsf.kappaWF_),
    A_(rwfpsf.A_),
    B_(rwfpsf.B_),
    C_(rwfpsf.C_),
    D_(rwfpsf.D_),
    Etke_(rwfpsf.Etke_),
    CmuMax_(rwfpsf.CmuMax_),
    blendCoeff_(rwfpsf.blendCoeff_),
    zDir_(rwfpsf.zDir_),
    roughWall_(rwfpsf.roughWall_),
    thresholdU_(rwfpsf.thresholdU_),
    thresholdE_(rwfpsf.thresholdE_),
    thresholdK_(rwfpsf.thresholdK_),
    constantU_(rwfpsf.constantU_),
    constantE_(rwfpsf.constantE_),
    constantK_(rwfpsf.constantK_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void NutRoughABLWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    nutkWallFunctionFvPatchScalarField::autoMap(m);
    z0_.autoMap(m);
    ustar_.autoMap(m);
    kappaWF_.autoMap(m);
    A_.autoMap(m);
    B_.autoMap(m);
    C_.autoMap(m);
    D_.autoMap(m);
    Etke_.autoMap(m);
    CmuMax_.autoMap(m);
    blendCoeff_.autoMap(m);
    zDir_.autoMap(m);
    roughWall_.autoMap(m);
    thresholdU_.autoMap(m);
    thresholdE_.autoMap(m);
    thresholdK_.autoMap(m);
    constantU_.autoMap(m);
    constantE_.autoMap(m);
    constantK_.autoMap(m);
}


void NutRoughABLWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    nutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const NutRoughABLWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const NutRoughABLWallFunctionFvPatchScalarField>(ptf);

    z0_.rmap(nrwfpsf.z0_, addr);
    ustar_.rmap(nrwfpsf.ustar_, addr);
    kappaWF_.rmap(nrwfpsf.kappaWF_, addr);
    A_.rmap(nrwfpsf.A_, addr);
    B_.rmap(nrwfpsf.B_, addr);
    C_.rmap(nrwfpsf.C_, addr);
    D_.rmap(nrwfpsf.D_, addr);
    Etke_.rmap(nrwfpsf.Etke_, addr);
    CmuMax_.rmap(nrwfpsf.CmuMax_, addr);
    blendCoeff_.rmap(nrwfpsf.blendCoeff_, addr);
    zDir_.rmap(nrwfpsf.zDir_, addr);
    roughWall_.rmap(nrwfpsf.roughWall_, addr);
    thresholdU_.rmap(nrwfpsf.thresholdU_, addr);
    thresholdE_.rmap(nrwfpsf.thresholdE_, addr);
    thresholdK_.rmap(nrwfpsf.thresholdK_, addr);
    constantU_.rmap(nrwfpsf.constantU_, addr);
    constantE_.rmap(nrwfpsf.constantE_, addr);
    constantK_.rmap(nrwfpsf.constantK_, addr);
}


void NutRoughABLWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
    os.writeEntry("z0", z0_);
    os.writeEntry("ustar", ustar_);
    os.writeEntry("kappaWF", kappaWF_);
    os.writeEntry("A", A_);
    os.writeEntry("B", B_);
    os.writeEntry("C", C_);
    os.writeEntry("D", D_);
    os.writeEntry("Etke", Etke_);
    os.writeEntry("CmuMax", CmuMax_);
    os.writeEntry("blendCoeff", blendCoeff_);
    os.writeEntry("zDir", zDir_);
    os.writeEntry("roughWall", roughWall_);
    os.writeEntry("thresholdU", thresholdU_);
    os.writeEntry("thresholdE", thresholdE_);
    os.writeEntry("thresholdK", thresholdK_);
    os.writeEntry("constantU", constantU_);
    os.writeEntry("constantE", constantE_);
    os.writeEntry("constantK", constantK_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    NutRoughABLWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
