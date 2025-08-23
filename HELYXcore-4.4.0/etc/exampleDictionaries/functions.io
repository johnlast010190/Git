partialWrite
{
    // Write some registered objects more often than others.
    // Above writeControl determines most frequent dump.

    type            partialWrite;

    // Where to load it from
    functionObjectLibs ("libIOFunctionObjects.so");

    // Optional mesh region to operate on. Only one partialWrite per
    // region allowed.
    region wallFilmRegion;

    // Execute upon options:
    //  timeStep
    //  outputTime
    //  adjustableTime
    //  runTime
    //  clockTime
    //  cpuTime
    outputControl   outputTime;

    // Objects (fields or lagrangian fields in any of the clouds)
    // to write every outputTime
    objectNames    (p positions nParticle);

    // Write as normal every writeInterval'th outputTime.
    outputInterval  1; // (timeStep, outputTime)

    // Interval of time (sec) to write down(
    writeInterval   10.5 //(adjustableTime, runTime, clockTime, cpuTime)
}

dumpObjects
{
    // Forcibly write registered objects

    type            writeRegisteredObject;

    // Where to load it from
    functionObjectLibs ("libIOFunctionObjects.so");

    // When to write:
    //  timeStep            (with optional outputInterval)
    //  outputTime          (with optional outputInterval)
    //  adjustableTime
    //  runTime
    //  clockTime
    //  cpuTime
    outputControl   outputTime;

    // Write every writeInterval (only valid for timeStemp, outputTime)
    outputInterval  1;

    // Interval of time (valid for adjustableTime, runTime, clockTime,
    //  cpuTime)
    writeInterval   10.5;

    // Objects to write
    objectNames    ();


    // Is the object written by this function Object alone
    // (default is false)
    //exclusiveWriting       true;
}
