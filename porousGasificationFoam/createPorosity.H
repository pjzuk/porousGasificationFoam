    Info<< "Creating porous zones" << nl << endl;
    porousReactingZone pZone(mesh,porosityF);
    Switch pressureImplicitPorosity(false);
