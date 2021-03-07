Info<< "\nCreating pyrolysis model" << endl;

typedef heterogeneousPyrolysisModels::heterogeneousPyrolysisModel pyrolysisModelType;

autoPtr<heterogeneousPyrolysisModels::heterogeneousPyrolysisModel> pyrolysis
(
    pyrolysisModelType::New(mesh,chemistry,porosityF,radiationF)
);

pyrolysisModelType& pyrolysisZone = pyrolysis();

volScalarField Srho
(
    IOobject
    (
        "Srho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pyrolysisZone.Srho()()
);

volScalarField heatTransfer
(
    IOobject
    (
        "heatTransfer",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pyrolysisZone.heatTransfer()()
);

