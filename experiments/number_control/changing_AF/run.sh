cd master || exit
python initialize.py
cd ..
cd results || exit
find . -name "simulate.py" | parallel python {}

# to use:
# chmod +x run.py
# time ./run.py
