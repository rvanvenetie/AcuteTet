#include <stdio.h>
#include <Judy.h>

#define MAXLINE 1000000                 // max string (line) length

uint8_t   Index[MAXLINE];               // string to insert

int     // Usage:  JudySort < file_to_sort
main()
{
Word_t   Index;                     // array index
Word_t   Value;                     // array element value
Word_t * PValue;                    // pointer to array element value
int      Rc_int;                    // return code

Pvoid_t  PJLArray = (Pvoid_t) NULL; // initialize JudyL array

while (scanf("%lu %lu", &Index, &Value))
{
    JLI(PValue, PJLArray, Index);

    *PValue = Value;                 // store new value
}
// Next, visit all the stored indexes in sorted order, first ascending,
// then descending, and delete each index during the descending pass.

Index = 0;
JLF(PValue, PJLArray, Index);
while (PValue != NULL)
{
    printf("%lu %lu\n", Index, *PValue);
    JLN(PValue, PJLArray, Index);
}

Index = -1;
JLL(PValue, PJLArray, Index);
while (PValue != NULL)
{
    printf("%lu %lu\n", Index, *PValue);

    JLD(Rc_int, PJLArray, Index);


    JLP(PValue, PJLArray, Index);
}
}
