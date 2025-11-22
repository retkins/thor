#include <stdio.h>
#include <stdlib.h>

struct Test {
    double *x;
};

int main() {

    struct Test *test = calloc(1, sizeof(struct Test)); 

    if (test->x == NULL) {
        printf("Calloc made *x NULL\n");
    }
    else {
        printf("*x is not NULL\n");
    }

    free(test);

    return 0;

}