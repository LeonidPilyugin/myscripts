curdir=$PWD

for subdir in $PWD; do
    export PATH=$PATH:$PWD/$subdir

    if [ -f $PWD/$subdir/add-to-path.sh ]; then
        cd $PWD/$subdir
        source ./add-to-path.sh
        cd ..
    fi
done

cd $curdir
