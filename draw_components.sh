for i in `ls -Sr *.dot`; do
    echo $i
    neato -Gorientation=landscape -s0.1 -Tpng -O -Kcirco $i
done
