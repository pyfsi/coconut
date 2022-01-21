void waitForSync(word message)
{
   int dummy;
   Pstream stream(Pstream::commsTypes::blocking);
   stream.scatter(dummy);  // sync: wait until all processors reach this point

   if (Pstream::master())
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
