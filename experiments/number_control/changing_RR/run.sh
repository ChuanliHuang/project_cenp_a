cd master || exit
python initialize.py
cd ..
cd results || exit
find . -name "simulate.py" | parallel python {}

