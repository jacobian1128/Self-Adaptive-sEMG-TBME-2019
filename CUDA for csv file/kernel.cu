#include <stdio.h>
#include <conio.h>

#include <Windows.h>

#include "TFcore.cuh"

int main()
{
    TFcore TF;
    TF.loadData("data.csv");

    while (TF.getRemaining() > 0) {
        printf("process %5.2f [%%]: M %d ",
            100 * (1 - TF.getRemaining() / (double)TF.getSamples()), TF.getPatterns());

        TF.getSample();

        TF.proceed();
        printf("\n");
    }

    Sleep(10 * 1000);
    return 0;
}