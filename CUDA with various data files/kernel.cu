#include <stdio.h>
#include <conio.h>

#include "TFcore.cuh"

int main()
{
    TFcore TF;
    TF.loadData("data.csv");

    while (!TF.isEmpty()) {
        TF.getSample();
        TF.proceed();
    }

    return 0;
}