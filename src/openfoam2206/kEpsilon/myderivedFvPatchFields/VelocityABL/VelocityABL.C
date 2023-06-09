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

#include "VelocityABL.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

VelocityABL::
VelocityABL
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchVectorField(p, iF),
    mainABL()
{}


VelocityABL::
VelocityABL
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchVectorField(p, iF),
    mainABL(patch().Cf(), dict)
{
    phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    refValue() = U(patch().Cf());
    refGrad() = Zero;
    valueFraction() = 1;

    if (dict.found("value"))
    {
        vectorField::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        vectorField::operator=(refValue());
    }
}


VelocityABL::
VelocityABL
(
    const VelocityABL& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchVectorField(pvf, p, iF, mapper),
    mainABL(pvf, mapper)
{}


VelocityABL::
VelocityABL
(
    const VelocityABL& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchVectorField(pvf, iF),
    mainABL(pvf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void VelocityABL::autoMap
(
    const fvPatchFieldMapper& m
)
{
    inletOutletFvPatchVectorField::autoMap(m);
    mainABL::autoMap(m);
}


void VelocityABL::rmap
(
    const fvPatchVectorField& pvf,
    const labelList& addr
)
{
    inletOutletFvPatchVectorField::rmap(pvf, addr);

    const VelocityABL& blpvf =
        refCast<const VelocityABL>(pvf);

    mainABL::rmap(blpvf, addr);
}


void VelocityABL::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    mainABL::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    VelocityABL
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
