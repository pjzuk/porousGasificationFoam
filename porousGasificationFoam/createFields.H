Info<< nl << "Reading thermophysicalProperties" << endl;

autoPtr<psiChemistryModel> pChemistry
(
    psiChemistryModel::New(mesh)
);

psiChemistryModel& chemistry = pChemistry();

hsCombustionThermo& thermo = chemistry.thermo();

basicMultiComponentMixture& composition = thermo.composition();	

PtrList<volScalarField>& Y = composition.Y();

word inertSpecie(thermo.lookup("inertSpecie"));

Switch solveEnergy(thermo.lookupOrDefault("solveEnergy",true));

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();
volScalarField& hs = thermo.hs();
const volScalarField& T = thermo.T();

#include "compressibleCreatePhi.H"

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(hs);

volScalarField chemistrySh
(
    IOobject
    (
        "gasChemistrySh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("chemistrySh", dimEnergy/dimTime/dimVolume, 0.0)
);

volScalarField ddtKinEn
(
    IOobject
    (
        "ddtKinEn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", chemistrySh.dimensions(), 0.0)
);

volScalarField divPhiKinEn
(
    IOobject
    (
        "divPhiKinEn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", chemistrySh.dimensions(), 0.0)
);


Info << "Reading field Df\n" << endl;

volTensorField Df
(
        IOobject
        (
            "Df",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
);

Info << "Reading field porosityF\n" << endl;

volScalarField porosityF
(
        IOobject
        (
            "porosityF",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
);

volScalarField radiationF
(
        IOobject
        (
            "radiationF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	    dimensionedScalar("zero",dimMass/dimLength/pow3(dimTime),0.0)
);


volScalarField cellVolume
(
	IOobject
	(
		"cellVolume",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("zero", dimVolume, 0.0)
);

forAll(cellVolume,cellI)
{
	cellVolume[cellI] = mesh.V()[cellI];
}


tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.schemesDict().divScheme("div(phi,Yi)")
    )
);

mesh.schemesDict().setFluxRequired(p.name());
