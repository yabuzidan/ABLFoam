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
scalar nutRoughABLWallFunctionFvPatchScalarField::blendingFunction
(
    const scalar BIA,
    const scalar blendCoeff
) const
{
    scalar transErr = Foam::constant::mathematical::pi*max(BIA-0.5,-0.5);
        /* Calculate and return blender value */
    return pow(1.0 - 0.5*(1.0 + sin(transErr)),blendCoeff);
}


tmp<scalarField> nutRoughABLWallFunctionFvPatchScalarField::nut() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
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
    const scalarField hombetaStar = min(pow(ustar_,4)/sqr(kref),betaStarMax_);
    const scalarField oref = (ustar_/sqrt(hombetaStar))/(kappaWF_*coord2);

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
    scalarField Oerr = mag((oref - (epsilon/(hombetaStar*k)))/max(oref,oref*1e-14))/constantO_;
    
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
    
    forAll(Oerr, facei)
  
    {
        if(Oerr[facei] < thresholdO_[facei])
        {
            Oerr[facei] = scalar(0);
        }
        if(Oerr[facei] > scalar(1))
        {
            Oerr[facei] = scalar(1);
        }
    }

    scalarField herr = max(max(Uerr,kerr),Oerr);
    
   
    forAll(nutw, facei)
    {
        label celli = patch().faceCells()[facei];
        
        
        
        const scalar betaStar25 = pow025(hombetaStar[facei]);

        scalar uStar = betaStar25*sqrt(k[celli]);
        
        if(roughWall_[facei] > 0) // Rough case
        {
            scalar yPlus = uStar*(y[facei]+z0_[facei])/nuw[facei];
            scalar Edash = nuw[facei]/(z0_[facei]*uStar);
            nutw[facei] =
            nuw[facei]*(yPlus*kappaWF_[facei]/log(Edash*yPlus) - 1);
            
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
            scalar Edash = E_;
            nutw[facei] =
            nuw[facei]*(yPlus*kappaWF_[facei]/log(Edash*yPlus) - 1);
            
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

nutRoughABLWallFunctionFvPatchScalarField::
nutRoughABLWallFunctionFvPatchScalarField
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
    betaStarMax_(p.size(), 0.0),
    blendCoeff_(p.size(), 0.0),
    zDir_(p.size(), vector(0, 0, 1)),
    roughWall_(p.size(), 1.0),
    thresholdU_(p.size(), 0.3),
    thresholdO_(p.size(), 0.5),
    thresholdK_(p.size(), 0.5),
    constantU_(p.size(), 1),
    constantO_(p.size(), 3),
    constantK_(p.size(), 10)
    
{}


nutRoughABLWallFunctionFvPatchScalarField::
nutRoughABLWallFunctionFvPatchScalarField
(
    const nutRoughABLWallFunctionFvPatchScalarField& ptf,
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
    betaStarMax_(mapper(ptf.betaStarMax_)),
    blendCoeff_(mapper(ptf.blendCoeff_)),
    zDir_(mapper(ptf.zDir_)),
    roughWall_(mapper(ptf.roughWall_)),
    thresholdU_(mapper(ptf.thresholdU_)),
    thresholdO_(mapper(ptf.thresholdO_)),
    thresholdK_(mapper(ptf.thresholdK_)),
    constantU_(mapper(ptf.constantU_)),
    constantO_(mapper(ptf.constantO_)),
    constantK_(mapper(ptf.constantK_))
{}


nutRoughABLWallFunctionFvPatchScalarField::
nutRoughABLWallFunctionFvPatchScalarField
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
    betaStarMax_("betaStarMax", dict, p.size()),
    blendCoeff_("blendCoeff", dict, p.size()),
    zDir_("zDir", dict, p.size()),
    roughWall_("roughWall", dict, p.size()),
    thresholdU_("thresholdU", dict, p.size()),
    thresholdO_("thresholdO", dict, p.size()),
    thresholdK_("thresholdK", dict, p.size()),
    constantU_("constantU", dict, p.size()),
    constantO_("constantO", dict, p.size()),
    constantK_("constantK", dict, p.size())
{}


nutRoughABLWallFunctionFvPatchScalarField::
nutRoughABLWallFunctionFvPatchScalarField
(
    const nutRoughABLWallFunctionFvPatchScalarField& rwfpsf
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
    betaStarMax_(rwfpsf.betaStarMax_),
    blendCoeff_(rwfpsf.blendCoeff_),
    zDir_(rwfpsf.zDir_),
    roughWall_(rwfpsf.roughWall_),
    thresholdU_(rwfpsf.thresholdU_),
    thresholdO_(rwfpsf.thresholdO_),
    thresholdK_(rwfpsf.thresholdK_),
    constantU_(rwfpsf.constantU_),
    constantO_(rwfpsf.constantO_),
    constantK_(rwfpsf.constantK_)
{}


nutRoughABLWallFunctionFvPatchScalarField::
nutRoughABLWallFunctionFvPatchScalarField
(
    const nutRoughABLWallFunctionFvPatchScalarField& rwfpsf,
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
    betaStarMax_(rwfpsf.betaStarMax_),
    blendCoeff_(rwfpsf.blendCoeff_),
    zDir_(rwfpsf.zDir_),
    roughWall_(rwfpsf.roughWall_),
    thresholdU_(rwfpsf.thresholdU_),
    thresholdO_(rwfpsf.thresholdO_),
    thresholdK_(rwfpsf.thresholdK_),
    constantU_(rwfpsf.constantU_),
    constantO_(rwfpsf.constantO_),
    constantK_(rwfpsf.constantK_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nutRoughABLWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    nutkWallFunctionFvPatchScalarField::autoMap(m);
    m(z0_, z0_);
    m(ustar_, ustar_);
    m(kappaWF_, kappaWF_);
    m(A_,A_);
    m(B_,B_);
    m(C_,C_);
    m(D_,D_);
    m(Etke_,Etke_);
    m(betaStarMax_,betaStarMax_);
    m(blendCoeff_,blendCoeff_);
    m(zDir_,zDir_);
    m(roughWall_,roughWall_);
    m(thresholdU_,thresholdU_);
    m(thresholdO_,thresholdO_);
    m(thresholdK_,thresholdK_);
    m(constantU_,constantU_);
    m(constantO_,constantO_);
    m(constantK_,constantK_);
}


void nutRoughABLWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    nutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const nutRoughABLWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const nutRoughABLWallFunctionFvPatchScalarField>(ptf);

    z0_.rmap(nrwfpsf.z0_, addr);
    ustar_.rmap(nrwfpsf.ustar_, addr);
    kappaWF_.rmap(nrwfpsf.kappaWF_, addr);
    A_.rmap(nrwfpsf.A_, addr);
    B_.rmap(nrwfpsf.B_, addr);
    C_.rmap(nrwfpsf.C_, addr);
    D_.rmap(nrwfpsf.D_, addr);
    Etke_.rmap(nrwfpsf.Etke_, addr);
    betaStarMax_.rmap(nrwfpsf.betaStarMax_, addr);
    blendCoeff_.rmap(nrwfpsf.blendCoeff_, addr);
    zDir_.rmap(nrwfpsf.zDir_, addr);
    roughWall_.rmap(nrwfpsf.roughWall_, addr);
    thresholdU_.rmap(nrwfpsf.thresholdU_, addr);
    thresholdO_.rmap(nrwfpsf.thresholdO_, addr);
    thresholdK_.rmap(nrwfpsf.thresholdK_, addr);
    constantU_.rmap(nrwfpsf.constantU_, addr);
    constantO_.rmap(nrwfpsf.constantO_, addr);
    constantK_.rmap(nrwfpsf.constantK_, addr);
}


void nutRoughABLWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry(os, "z0", z0_);
    writeEntry(os, "ustar", ustar_);
    writeEntry(os, "kappaWF", kappaWF_);
    writeEntry(os, "value", *this);
    writeEntry(os,"A", A_);
    writeEntry(os,"B", B_);
    writeEntry(os,"C", C_);
    writeEntry(os,"D", D_);
    writeEntry(os,"Etke", Etke_);
    writeEntry(os,"betaStarMax", betaStarMax_);
    writeEntry(os,"blendCoeff", blendCoeff_);
    writeEntry(os,"zDir", zDir_);
    writeEntry(os,"roughWall", roughWall_);
    writeEntry(os,"thresholdU", thresholdU_);
    writeEntry(os,"thresholdO", thresholdO_);
    writeEntry(os,"thresholdK", thresholdK_);
    writeEntry(os,"constantU", constantU_);
    writeEntry(os,"constantO", constantO_);
    writeEntry(os,"constantK", constantK_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutRoughABLWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
