#include "AddMe.h"

int main(void)
{
    int num1, num2, sum;
    cout << "\n\n Enters two numbers to add: " << endl;
    cin >> num1 >> num2;

    // Function call
    sum = add(num1, num2);
    cout << "\n Sum: " << num1 << " and "
         << num2 << " is " <<  sum << "\n\n";

    return 0;
}
