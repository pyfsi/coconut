void waitForSync(word message)
{
    // using separate blocks ensures that the labels are restricted to their own scope
    // the sequential use of nProcStart and nProcEnd avoids that for one of the processors the label from a previous
    // call still lives and is mistaken for the current one
    {
        label nProcsStart = 1; // use different label names for start and end, to avoid mix up!
        label totNprocs = returnReduce(nProcsStart, sumOp<label>()); // sync (gather): wait until all processors reach this point

        if (Pstream::master() && (totNprocs == Pstream::nProcs()))
        {
            word extension = ".coco";
            word ready = "_ready";
            fileName file;

            file = message + ready + extension;
            OFstream outfile (file);

            file = message + extension;
            rm(file);
        }
    }

    {
        label nProcsEnd = 1;
        returnReduce(nProcsEnd, sumOp<label>()); // sync (gather): wait until all processors reach this point
    }
}
