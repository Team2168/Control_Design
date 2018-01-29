funtion [Ders] = DersBasisFuns(u, j, p, n, U)

% Inputs: u - value of the independent variable
% j - index of the knot span, which includes u
% p - degree of the spline
% n - max degree of differentiation of B-spline
% basis functions
% U[] - Knot vector
% Output: Ders[][] - values of B-spline basis functions and theirs
% derivatives at u

{




double DR[MAX_P], DL[MAX_P];
Matrix Du, a;
double acc, temp, d;
int i, r, k, s1, s2, rk, pk, i1, i2;
Du[0][0] = 1.0;
for (j=1; j<=p; j++)
{
DL[j] = u - U[i+1-j];
DR[j] = U[i+j]-u;
acc = 0.0;
for (r=0;r<j;r++)
{
    Du[j][r] = DR[r+1] + DL[j-r];
temp = Du[r][j-1] / Du[j][r];
Du[r][j] = acc + DR[r+1] * temp;
acc = DL[j-r] * temp;
}
Du[j][j] = acc;
}
for (j=0; j<=p; j++)
Ders[0][j] = Du[j][p];

for (r=0; r<=p; r++)
{
s1=0;
s2=1;
a[0][0] = 1.0;
for (k=1; k<=n; k++)
{
d = 0.0;
rk = r - k;
pk = p - k;
if (r >= k)
{
a[s2][0] = a[s1][0] / Du[pk+1][rk];
d = a[s2][0] * Du[rk][pk];
}
if (rk >= -1)
j1 = 1;
else
j1 = -rk;
if (r-1 <= pk)
j2 = k - 1;
else
j2 = p - r;
for (j=j1; j<=j2; j++)
{
a[s2][j] = (a[s1][j] - a[s1][j-1]) / Du[pk+1][rk+j];
d += a[s2][j] * Du[rk+j][pk];
}
if (r <= pk)
{
a[s2][k] = -a[s1][k-1] / Du[pk+1][r];
d += a[s2][k] * Du[r][pk];
}
Ders[k][r] = d;
j = s1; s1 = s2; s2 = j;
}
}
r = p;
for (k=1; k<=n; k++)
{
for (j=0; j<=p; j++) Ders[k][j] *= r;
r *= (p-k);
}
}