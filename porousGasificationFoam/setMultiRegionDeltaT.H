if (adjustTimeStep)
{

    if (CoNum == -GREAT)
    {
        CoNum = SMALL;
    }

    if (DiNum == -GREAT)
    {
        DiNum = SMALL;
    }

    const scalar TFactorFluid = maxCo/(CoNum + SMALL);
    const scalar TFactorSolid = maxDi/(DiNum + SMALL);
    const scalar dt0 = runTime.deltaT().value()*
       min(
           min(
	       min(TFactorFluid,1.0+0.1*TFactorFluid),
               min(TFactorSolid,1.0+0.1*TFactorSolid)
	       ),
           1.2
	  );
    runTime.setDeltaT
    (
        min(dt0,maxDeltaT)
    );
}
