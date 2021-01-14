cd CSM
ls -1 | egrep -v '(makeHostFile.sh|AbaqusHosts.txt|Base.inp)' | xargs -I files rm -r "files"

