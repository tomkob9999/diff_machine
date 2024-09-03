# diff_machine
Calculates difference equations.  It can handle both additive and multiplicative equations.

This can handle only differential equations that converges after rounds of differentiations or divisions with the previous value.  This cannot solve difference equations like Fibonacci like y(x) = y(x-1)+y(x-2) which does not seem to converge after rounds of subtractions or divisions with the previous value.

For example, "y=4x^3+3x^2" converges as below after 3 orders of differencing.  It can be expressed with either in the form of using only y(x), y(x-1), y(x-2)... or with combinations of using x.

| 0 |  0  |      |      |     |
|---|-----|------|------|-----|
| 1 |  7  |  7   |      |     |
| 2 | 44  |  37  |  30  |     |
| 3 | 135 |  91  |  54  | 24  |
| 4 | 304 | 169  |  78  | 24  |
| 5 | 575 | 271  | 102  | 24  |
| 6 | 972 | 397  | 126  | 24  |

Here is the comparison of drawing y=x**2 by using differential equation y'=2x and the differential equation.  The difference equation does not lose precision even if the discretization step is wide whereas the derivative loses precision significantly, even though the difference equation does not use the x value to find y value.

### Differential equation approach
![derivative](https://github.com/user-attachments/assets/46755d52-0006-4c5a-b716-ea574dcb2a77)

### Difference equation approach
![difference](https://github.com/user-attachments/assets/f192faaf-a4a7-40e4-8f0e-92e820449ee6)
