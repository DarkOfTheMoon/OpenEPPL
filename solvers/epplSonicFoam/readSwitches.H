const dictionary& epplDict = eppl.dict();

const scalar pressureRDC
(
    epplDict.lookupOrDefault<scalar>("pressureRDC", 1.)
);

const scalar momentumResCritia
(
    epplDict.lookupOrDefault<scalar>("momentumResCritia", 1E-02)
);
