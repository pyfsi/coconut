void waitForSync(word message)
{
   label nProcs = 1;
   label totNprocs = returnReduce(nProcs, sumOp<label>()); // sync (gather): wait until all processors reach this point

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
