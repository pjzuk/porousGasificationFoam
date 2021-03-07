    autoPtr<radiation::heterogeneousRadiationModel> radiation
    (
        radiation::heterogeneousRadiationModel::New(thermo.T(),porosityF,pyrolysisZone.surf(),pyrolysisZone.T())
    );
