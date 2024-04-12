#include "OFstream.H"

void sendMessage(word message)
{
    if (Pstream::master())
    {
        OFstream outfile (message);
    }
}

void removeMessage(word message)
{
    if (Pstream::master())
    {
       rm(message);
    }
}

void waitForSync(word message)
{
    // using separate blocks ensures that the labels are restricted to their own scope
    // the sequential use of nProcStart and nProcEnd avoids that for one of the processors the label from a previous
    // call still lives and is mistaken for the current one
    {
        label nProcsStart = 1; // use different label names for start and end, to avoid mix up
        label totNprocs = returnReduce(nProcsStart, sumOp<label>()); // sync (gather): wait until all processors reach this point

        if (Pstream::master() && (totNprocs == Pstream::nProcs()))
        {
            sendMessage(message + "_ready.coco");
            removeMessage(message + ".coco");
        }
    }

    {
        label nProcsEnd = 1;
        returnReduce(nProcsEnd, sumOp<label>()); // sync (gather): wait until all processors reach this point
    }
}
