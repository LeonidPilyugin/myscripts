curdir=$PWD

for subdir in $(ls -d -- */); do
    export PATH=$PATH:$PWD/$subdir

    if [ -f $subdir/add-to-path.sh ]; then
        cd $subdir
        source ./add-to-path.sh
        cd ..
    fi
done

cd $curdir
