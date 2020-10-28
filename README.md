# ZetaPlaneCurves

Computes L-polynomial of the projective normalisation of the curve C defined by the zero locus of f in the projective plane, where C has at most nodal singularities.
The L-polynomial is obtained by point counts of normalisation, deduced from the point counts on the plane model obtained via Harvey's trace formula, Theorem 3.1 in "Computing zeta functions of arithmetic schemes".


This package defines the following intrinsics:
  - `LPolynomial(f::RngMPolElt : KnownFactor:=false, corrections:=false) -> RngUPolElt`
  - `LPolynomial(f::RngMPolElt, p::RngIntElt : KnownFactor:=false, corrections:=false) -> RngUPolElt`


For genus at least 3 and moderate prime should outperform Magma under generic conditions.
The current implementation computes f^(p-1) and its powers in a naive fashion.
This takes about p^(2 + o(1)) deg(f)^(g/2 + o(1))) time and space, which makes this algorithm unpractical for large characteristic.


```
> AttachSpec(spec");
> SetIgnorePrompttrue); // so you can copy paste without worrying about the ">"
> R<a,b,c> := PolynomialRing(Integers(), 3);
//genus 5 curve that is a quartic cover of the elliptic curve 330e2
> f := 900*a^6 - 8236*a^4*c^2 + 784*a^3*b^2*c + 11*a^2*b^4 + 1228*a^2*c^4 - 16*a*b^2*c^3 - 11*b^4*c^2 - 36*c^6;
> time LPolynomial(f, 19);
2476099*T^10 + 1042568*T^9 + 541861*T^8 + 173280*T^7 + 52478*T^6 + 12976*T^5 + 2762*T^4 + 480*T^3 + 79*T^2 + 8*T + 1
Time: 0.390
> time LPolynomial(f, 53);
418195493*T^10 - 110466734*T^9 + 4912941*T^8 + 786520*T^7 + 246874*T^6 - 80052*T^5 + 4658*T^4 + 280*T^3 + 33*T^2 - 14*T + 1
Time: 0.890
> time LPolynomial(f, 109);
15386239549*T^10 - 1976214254*T^9 + 446785005*T^8 - 50660584*T^7 + 7023306*T^6 - 630644*T^5 + 64434*T^4 - 4264*T^3 + 345*T^2 - 14*T + 1
Time: 3.720
> time LPolynomial(f, 211);
418227202051*T^10 + 47570866584*T^9 + 4048784261*T^8 + 235070880*T^7 + 25301854*T^6 + 1763856*T^5 + 119914*T^4 + 5280*T^3 + 431*T^2 + 24*T + 1
Time: 17.880
```
One may speed up the computation by providing a known factor of the L-polynomial
```
> time LPolynomial(f, 19: KnownFactor:=EulerFactor(EllipticCurve("330e2"), 19));
2476099*T^10 + 1042568*T^9 + 541861*T^8 + 173280*T^7 + 52478*T^6 + 12976*T^5 + 2762*T^4 + 480*T^3 + 79*T^2 + 8*T + 1
Time: 0.150
> time LPolynomial(f, 53: KnownFactor:=EulerFactor(EllipticCurve("330e2"), 53));
418195493*T^10 - 110466734*T^9 + 4912941*T^8 + 786520*T^7 + 246874*T^6 - 80052*T^5 + 4658*T^4 + 280*T^3 + 33*T^2 - 14*T + 1
Time: 0.280
> time LPolynomial(f, 109: KnownFactor:=EulerFactor(EllipticCurve("330e2"), 109));
15386239549*T^10 - 1976214254*T^9 + 446785005*T^8 - 50660584*T^7 + 7023306*T^6 - 630644*T^5 + 64434*T^4 - 4264*T^3 + 345*T^2 - 14*T + 1
Time: 1.340
> time LPolynomial(f, 211: KnownFactor:=EulerFactor(EllipticCurve("330e2"), 211));
418227202051*T^10 + 47570866584*T^9 + 4048784261*T^8 + 235070880*T^7 + 25301854*T^6 + 1763856*T^5 + 119914*T^4 + 5280*T^3 + 431*T^2 + 24*T + 1
Time: 6.600
```
In this case, since the defining polynomial is quite sparse, one may even consider significantly larger.
Here are some timings:
```
> time LPolynomial(f, 307);
2727042318307*T^10 + 25433375397*T^8 + 136890686*T^6 + 445898*T^4 + 879*T^2 + 1
Time: 38.930
> time LPolynomial(f, 401);
10368641602001*T^10 + 775708848030*T^9 + 69962103085*T^8 + 4380219240*T^7 + 279565170*T^6 + 15992948*T^5 + 697170*T^4 + 27240*T^3 + 1085*T^2 + 30*T + 1
Time: 31.520
> time LPolynomial(f, 503);
32198817702743*T^10 - 1536325297944*T^9 + 234546680261*T^8 - 8136769440*T^7 + 763728038*T^6 - 20263824*T^5 + 1518346*T^4 - 32160*T^3 + 1843*T^2 - 24*T + 1
Time: 50.120
> time LPolynomial(f, 601);
78410163603001*T^10 + 782796974406*T^9 + 202537320333*T^8 + 5990157384*T^7 + 494754018*T^6 + 7825316*T^5 + 823218*T^4 + 16584*T^3 + 933*T^2 + 6*T + 1
Time: 80.990
> time LPolynomial(f, 701);
169273934903501*T^10 + 14971446453662*T^9 + 1171549615501*T^8 + 67070339688*T^7 + 3383386314*T^6 + 132459668*T^5 + 4826514*T^4 + 136488*T^3 + 3401*T^2 + 62*T + 1
Time: 147.610
> time LPolynomial(f, 809);
346531411903049*T^10 - 12850361380830*T^9 + 782034765533*T^8 - 13115799240*T^7 + 472638834*T^6 - 1314420*T^5 + 584226*T^4 - 20040*T^3 + 1477*T^2 - 30*T + 1
Time: 170.200
> time LPolynomial(f, 907);
613813499121307*T^10 - 18949038561628*T^9 + 2022792705173*T^8 - 41553646288*T^7 + 2960413534*T^6 - 45576552*T^5 + 3263962*T^4 - 50512*T^3 + 2711*T^2 - 28*T + 1
Time: 175.660
> time LPolynomial(f, 1009);
1045817322864049*T^10 - 26948711986586*T^9 + 1229610743613*T^8 - 4976379928*T^7 + 337153314*T^6 - 3755804*T^5 + 334146*T^4 - 4888*T^3 + 1197*T^2 - 26*T + 1
Time: 269.430
> time LPolynomial(f, 2003);
2240721080810243*T^10 - 579463783778916*T^9 + 39079330733301*T^8 - 413076446640*T^7 + 19873978318*T^6 - 150966936*T^5 + 9922106*T^4 - 102960*T^3 + 4863*T^2 - 36*T + 1
Time: 1685.860
```


