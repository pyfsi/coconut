labelList globalFaceZones(0);

if (Pstream::parRun())
{
    SLList<label> globalFaceZonesSet;

    // philipc: directly lookup globalFaceZones from decomposeParDict
    IOdictionary decompDict
    (
        IOobject
        (
            "decomposeParDict",
            runTime.time().system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    if (decompDict.found("globalFaceZones"))
    {
        wordList globalFaceZoneNames(decompDict.lookup("globalFaceZones"));

        const faceZoneMesh& faceZones = mesh.faceZones();

        forAll(globalFaceZoneNames, nameI)
        {
            const label zoneID =
                faceZones.findZoneID(globalFaceZoneNames[nameI]);

            if (zoneID == -1)
            {
                FatalError
                    << "findGlobalFaceZones.H: cannot find globalFaceZone:"
                    << " " << globalFaceZoneNames[nameI] << abort(FatalError);
            }

            globalFaceZonesSet.insert(zoneID);
        }

        globalFaceZones = labelList(globalFaceZonesSet);
    }
}

// fix: read/write map for restarts, as comparing points might not work on
// restart due to write tolerance
// labelListList globalToLocalFaceZonePointMap
IOList<labelList> globalToLocalFaceZonePointMap
(
    IOobject
    (
        "globalToLocalFaceZonePointMap",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    globalFaceZones.size()
);

if (Pstream::parRun() && globalFaceZones.size() > 0)
{
    // If any of the map sizes are zero then we assume that the maps were
    // not read
    if (globalToLocalFaceZonePointMap[0].size() == 0)
    {
        forAll(globalFaceZones, zoneI)
        {
            label curZoneID = globalFaceZones[zoneI];

            Info<< "Creating faceMap for globalFaceZones "
                << mesh.faceZones()[curZoneID].name()<< endl;

            labelList curMap(mesh.faceZones()[curZoneID]().nPoints(), -1);

            vectorField fzGlobalPoints =
                mesh.faceZones()[curZoneID]().localPoints();

            // Set all slave points to zero because only the master order is
            // used
            if (!Pstream::master())
            {
                fzGlobalPoints *= 0.0;
            }

            //- pass points to all procs
            reduce(fzGlobalPoints, sumOp<vectorField>());

            // Now every proc has the master's list of FZ points
            // every proc must now find the mapping from their local FZ points
            // to the global FZ points

            const vectorField& fzLocalPoints =
                mesh.faceZones()[curZoneID]().localPoints();

            const edgeList& fzLocalEdges =
                mesh.faceZones()[curZoneID]().edges();

            const labelListList& fzPointEdges =
                mesh.faceZones()[curZoneID]().pointEdges();

            scalarField minEdgeLength(fzLocalPoints.size(), GREAT);

            forAll(minEdgeLength, pI)
            {
                const labelList& curPointEdges = fzPointEdges[pI];

                forAll(curPointEdges, eI)
                {
                    scalar Le =
                        fzLocalEdges[curPointEdges[eI]].mag(fzLocalPoints);
                    if (Le < minEdgeLength[pI])
                    {
                        minEdgeLength[pI] = Le;
                    }
                }
            }

            forAll(fzGlobalPoints, globalPointI)
            {
                //scalar minDist = GREAT;
                bool pointFound = false;

                forAll(fzLocalPoints, procPointI)
                {
                    scalar curDist =
                        mag
                        (
                            fzLocalPoints[procPointI]
                            - fzGlobalPoints[globalPointI]
                        );

                    // if (curDist < minDist)
                    // {
                    //     minDist = curDist;
                    // }

                    if (curDist < 1e-4*minEdgeLength[procPointI])
                    {
                        curMap[globalPointI] = procPointI;
                        //Info<< "curDist " << curDist
                        //<< " minDist " << minDist << endl;
                        //break;
                        if (pointFound)
                        {
                            FatalError
                                << "findGlobalFaceZones: point found twice!"
                                << abort(FatalError);
                        }
                        pointFound = true;
                    }
                }

                // if (curMap[globalPointI] == -1)
                // {
                //     Pout << "minDist: " << minDist << endl;
                // }
            }

            forAll(curMap, globalPointI)
            {
                if (curMap[globalPointI] == -1)
                {
                    FatalErrorIn(args.executable())
                        << "local to global face zone point map is not correct"
                        << " for zone " << zoneI
                        << abort(FatalError);
                }
            }

            globalToLocalFaceZonePointMap[zoneI] = curMap;
        }
    }
    else
    {
        Info<< "globalToLocalFaceZonePointMap read from file" << endl;
    }
}
