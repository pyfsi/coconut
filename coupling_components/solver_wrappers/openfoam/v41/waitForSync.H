void waitForSync(word message)
{
   int dummy;
   Pstream::scatter(dummy, Pstream::blocking);  // sync: wait until all processors reach this point

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
