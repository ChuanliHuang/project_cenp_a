cd master || exit
python initialize.py
cd ..
cd results || exit
# find . -name "simulate.py" | parallel python {}
for d in */ ; do
    python $d/"simulate.py" &
done
wait
# to use:
# chmod +x run.sh
# time ./run.sh
