# ZetaPlaneCurves

Computes L-polynomial of the projective normalisation of the curve C defined by the zero locus of f in the projective plane, where C has at most nodal singularities.
The L-polynomial is obtained by point counts of normalisation, deduced from the point counts on the plane model obtained via Harvey's trace formula, Theorem 3.1 in "Computing zeta functions of arithmetic schemes".


This package defines the following intrinsics:
  - `LPolynomial(f::RngMPolElt : KnownFactor:=false, corrections:=false) -> RngUPolElt`
  - `LPolynomial(f::RngMPolElt, p::RngIntElt : KnownFactor:=false, corrections:=false) -> RngUPolElt`


For genus at least 3 and moderate prime should outperform Magma under generic conditions.

The current implementation computes `f^(p-1)` and its powers in a naive fashion.
This takes about `p^(2 + o(1)) deg(f)^(g/2 + o(1)))` time and space, which makes this algorithm unpractical for large characteristic.

# Installation
You can enable the functionality of this package in Magma by attaching the `spec` file with AttachSpec.
```
> AttachSpec("<replace with the installation directory>/ZetaPlaneCurves/spec");
```
this package is part of [CHIMP](https://github.com/edgarcosta/CHIMP/).

# Usage


## Genus 5, degree 8

X is a genus 5 curve that is a quartic cover of the elliptic curve with [LMFDB](https://www.lmfdb.org) label [65.a1](https://www.lmfdb.org/EllipticCurve/Q/65/a/1) and Cremona label 65a1.
```
> R<a,b,c> := PolynomialRing(Integers(), 3);
> f := 12352*a^8 - 15792*a^7*b + 93440*a^7*c - 479*a^6*b^2 - 33520*a^6*b*c + 85696*a^6*c^2 + 1691*a^5*b^3 - 49410*a^5*b^2*c + 169000*a^5*b*c^2 - 303680*a^5*c^3 + 1507*a^4*b^4 - 4740*a^4*b^3*c + 4983*a^4*b^2*c^2 - 176880*a^4*b*c^3 + 352352*a^4*c^4 + 595*a^3*b^5 - 2928*a^3*b^4*c + 12401*a^3*b^3*c^2 - 10758*a^3*b^2*c^3 + 9410*a^3*b*c^4 - 238720*a^3*c^5 + 124*a^2*b^6 + 674*a^2*b^5*c + 812*a^2*b^4*c^2 - 17094*a^2*b^3*c^3 + 43019*a^2*b^2*c^4 + 90590*a^2*b*c^5 + 122600*a^2*c^6 - a*b^7 + 262*a*b^6*c - 1173*a*b^5*c^2 + 2342*a*b^4*c^3 + 1815*a*b^3*c^4 - 45420*a*b^2*c^5 - 63450*a*b*c^6 - 39000*a*c^7 + 3*b^8 + 48*b^7*c - 116*b^6*c^2 + 616*b^5*c^3 + 172*b^4*c^4 + 3240*b^3*c^5 + 11900*b^2*c^6 + 12000*b*c^7 + 5625*c^8;
> X := Curve(ProjectiveSpace(R), f);
> time LPolynomial(f, 11);
161051*T^10 - 29282*T^9 + 6655*T^8 + 1452*T^7 - 1672*T^6 + 172*T^5 - 152*T^4 + 12*T^3 + 5*T^2 - 2*T + 1
Time: 1.140
> time LPolynomial(f, 17);
1419857*T^10 - 167042*T^9 + 142477*T^8 - 6936*T^7 + 4114*T^6 - 76*T^5 + 242*T^4 - 24*T^3 + 29*T^2 - 2*T + 1
Time: 1.300
> time LPolynomial(f, 19);
2476099*T^10 + 781926*T^9 + 253783*T^8 + 38988*T^7 + 6384*T^6 - 36*T^5 + 336*T^4 + 108*T^3 + 37*T^2 + 6*T + 1
Time: 1.320
> time LPolynomial(f, 53);
418195493*T^10 + 15780962*T^9 + 12654545*T^8 + 719104*T^7 + 22790*T^6 + 16028*T^5 + 430*T^4 + 256*T^3 + 85*T^2 + 2*T + 1
Time: 2.740
> time LPolynomial(f, 109);
15386239549*T^10 - 846948966*T^9 + 99717233*T^8 + 3801920*T^7 - 1365770*T^6 + 129772*T^5 - 12530*T^4 + 320*T^3 + 77*T^2 - 6*T + 1
Time: 12.880
```
One may speed up the computation by providing a known factor of the L-polynomial.
```
> time LPolynomial(f, 11: KnownFactor:=EulerFactor(EllipticCurve("65a1"), 11));
161051*T^10 - 29282*T^9 + 6655*T^8 + 1452*T^7 - 1672*T^6 + 172*T^5 - 152*T^4 + 12*T^3 + 5*T^2 - 2*T + 1
Time: 1.070
> time LPolynomial(f, 17: KnownFactor:=EulerFactor(EllipticCurve("65a1"), 17));
1419857*T^10 - 167042*T^9 + 142477*T^8 - 6936*T^7 + 4114*T^6 - 76*T^5 + 242*T^4 - 24*T^3 + 29*T^2 - 2*T + 1
Time: 0.390
> time LPolynomial(f, 19: KnownFactor:=EulerFactor(EllipticCurve("65a1"), 19));
2476099*T^10 + 781926*T^9 + 253783*T^8 + 38988*T^7 + 6384*T^6 - 36*T^5 + 336*T^4 + 108*T^3 + 37*T^2 + 6*T + 1
Time: 0.390
> time LPolynomial(f, 53: KnownFactor:=EulerFactor(EllipticCurve("65a1"), 53));
418195493*T^10 + 15780962*T^9 + 12654545*T^8 + 719104*T^7 + 22790*T^6 + 16028*T^5 + 430*T^4 + 256*T^3 + 85*T^2 + 2*T + 1
Time: 1.350
```
Magma does not handle p = 11 or 17 and struggles with 19.
```
> time LPolynomial(ChangeRing(X, GF(11)));
...
Runtime error in 'SingularPoints': Singular locus of argument is not zero dimensional
> time LPolynomial(ChangeRing(X, GF(17)));
...
Runtime error in '/': Division by zero
> // and struggles with p=19
> time LPolynomial(ChangeRing(X, GF(19)));
2476099*x^10 + 781926*x^9 + 253783*x^8 + 38988*x^7 + 6384*x^6 - 36*x^5 + 336*x^4 + 108*x^3 + 37*x^2 + 6*x + 1
Time: 55.580
> time LPolynomial(ChangeRing(X, GF(53)));
418195493*x^10 + 15780962*x^9 + 12654545*x^8 + 719104*x^7 + 22790*x^6 + 16028*x^5 + 430*x^4 + 256*x^3 + 85*x^2 + 2*x + 1
Time: 58.050
> time LPolynomial(ChangeRing(X, GF(109)));
15386239549*x^10 - 846948966*x^9 + 99717233*x^8 + 3801920*x^7 - 1365770*x^6 + 129772*x^5 - 12530*x^4 + 320*x^3 + 77*x^2 - 6*x + 1
Time: 87.190
```


## Genus 5, degree 6
Y is a genus 5 curve that is a quartic cover of the elliptic curve with [LMFDB](https://www.lmfdb.org) label [330.b2](https://beta.lmfdb.org/EllipticCurve/Q/330/b/2) and Cremona label 330e2.
```
//genus 5 curve that is a quartic cover of the elliptic curve 330e2
> g := 900*a^6 - 8236*a^4*c^2 + 784*a^3*b^2*c + 11*a^2*b^4 + 1228*a^2*c^4 - 16*a*b^2*c^3 - 11*b^4*c^2 - 36*c^6;
> time LPolynomial(g, 19);
2476099*T^10 + 1042568*T^9 + 541861*T^8 + 173280*T^7 + 52478*T^6 + 12976*T^5 + 2762*T^4 + 480*T^3 + 79*T^2 + 8*T + 1
Time: 0.390
> time LPolynomial(g, 53);
418195493*T^10 - 110466734*T^9 + 4912941*T^8 + 786520*T^7 + 246874*T^6 - 80052*T^5 + 4658*T^4 + 280*T^3 + 33*T^2 - 14*T + 1
Time: 0.890
> time LPolynomial(g, 109);
15386239549*T^10 - 1976214254*T^9 + 446785005*T^8 - 50660584*T^7 + 7023306*T^6 - 630644*T^5 + 64434*T^4 - 4264*T^3 + 345*T^2 - 14*T + 1
Time: 3.720
> time LPolynomial(g, 211);
418227202051*T^10 + 47570866584*T^9 + 4048784261*T^8 + 235070880*T^7 + 25301854*T^6 + 1763856*T^5 + 119914*T^4 + 5280*T^3 + 431*T^2 + 24*T + 1
Time: 17.880
```
As before, one may speed up the computation by providing a known factor of the L-polynomial
```
> time LPolynomial(g, 19: KnownFactor:=EulerFactor(EllipticCurve("330e2"), 19));
2476099*T^10 + 1042568*T^9 + 541861*T^8 + 173280*T^7 + 52478*T^6 + 12976*T^5 + 2762*T^4 + 480*T^3 + 79*T^2 + 8*T + 1
Time: 0.150
> time LPolynomial(g, 53: KnownFactor:=EulerFactor(EllipticCurve("330e2"), 53));
418195493*T^10 - 110466734*T^9 + 4912941*T^8 + 786520*T^7 + 246874*T^6 - 80052*T^5 + 4658*T^4 + 280*T^3 + 33*T^2 - 14*T + 1
Time: 0.280
> time LPolynomial(g, 109: KnownFactor:=EulerFactor(EllipticCurve("330e2"), 109));
15386239549*T^10 - 1976214254*T^9 + 446785005*T^8 - 50660584*T^7 + 7023306*T^6 - 630644*T^5 + 64434*T^4 - 4264*T^3 + 345*T^2 - 14*T + 1
Time: 1.340
> time LPolynomial(g, 211: KnownFactor:=EulerFactor(EllipticCurve("330e2"), 211));
418227202051*T^10 + 47570866584*T^9 + 4048784261*T^8 + 235070880*T^7 + 25301854*T^6 + 1763856*T^5 + 119914*T^4 + 5280*T^3 + 431*T^2 + 24*T + 1
Time: 6.600
```

In this case, since the polynomial `g` is quite sparse, one may even consider significantly larger primes.
Here are some timings:
```
> Y := Curve(ProjectiveSpace(R), g);
> time LPolynomial(g, 307);
2727042318307*T^10 + 25433375397*T^8 + 136890686*T^6 + 445898*T^4 + 879*T^2 + 1
Time: 38.930
> time LPolynomial(g, 401);
10368641602001*T^10 + 775708848030*T^9 + 69962103085*T^8 + 4380219240*T^7 + 279565170*T^6 + 15992948*T^5 + 697170*T^4 + 27240*T^3 + 1085*T^2 + 30*T + 1
Time: 31.520
> time LPolynomial(g, 503);
32198817702743*T^10 - 1536325297944*T^9 + 234546680261*T^8 - 8136769440*T^7 + 763728038*T^6 - 20263824*T^5 + 1518346*T^4 - 32160*T^3 + 1843*T^2 - 24*T + 1
Time: 50.120
> time LPolynomial(g, 601);
78410163603001*T^10 + 782796974406*T^9 + 202537320333*T^8 + 5990157384*T^7 + 494754018*T^6 + 7825316*T^5 + 823218*T^4 + 16584*T^3 + 933*T^2 + 6*T + 1
Time: 80.990
> time LPolynomial(g, 701);
169273934903501*T^10 + 14971446453662*T^9 + 1171549615501*T^8 + 67070339688*T^7 + 3383386314*T^6 + 132459668*T^5 + 4826514*T^4 + 136488*T^3 + 3401*T^2 + 62*T + 1
Time: 147.610
> time LPolynomial(g, 809);
346531411903049*T^10 - 12850361380830*T^9 + 782034765533*T^8 - 13115799240*T^7 + 472638834*T^6 - 1314420*T^5 + 584226*T^4 - 20040*T^3 + 1477*T^2 - 30*T + 1
Time: 170.200
> time LPolynomial(g, 907);
613813499121307*T^10 - 18949038561628*T^9 + 2022792705173*T^8 - 41553646288*T^7 + 2960413534*T^6 - 45576552*T^5 + 3263962*T^4 - 50512*T^3 + 2711*T^2 - 28*T + 1
Time: 175.660
> time LPolynomial(g, 1009);
1045817322864049*T^10 - 26948711986586*T^9 + 1229610743613*T^8 - 4976379928*T^7 + 337153314*T^6 - 3755804*T^5 + 334146*T^4 - 4888*T^3 + 1197*T^2 - 26*T + 1
Time: 269.430
> time LPolynomial(g, 2003);
2240721080810243*T^10 - 579463783778916*T^9 + 39079330733301*T^8 - 413076446640*T^7 + 19873978318*T^6 - 150966936*T^5 + 9922106*T^4 - 102960*T^3 + 4863*T^2 - 36*T + 1
Time: 1685.860
```

In this example, Magma starts to outperform for p > 500 , but still fails at some primes.
```
> time LPolynomial(ChangeRing(Y, GF(19)));
2476099*x^10 + 1042568*x^9 + 541861*x^8 + 173280*x^7 + 52478*x^6 + 12976*x^5 + 2762*x^4 + 480*x^3 + 79*x^2 + 8*x + 1
Time: 2.210
> time LPolynomial(ChangeRing(Y, GF(53)));
...
Runtime error in '/': Division by zero
> time LPolynomial(ChangeRing(Y, GF(109)));
15386239549*x^10 - 1976214254*x^9 + 446785005*x^8 - 50660584*x^7 + 7023306*x^6 - 630644*x^5 + 64434*x^4 - 4264*x^3 + 345*x^2 - 14*x + 1
Time: 8.530
> time LPolynomial(ChangeRing(Y, GF(503)));
32198817702743*x^10 - 1536325297944*x^9 + 234546680261*x^8 - 8136769440*x^7 + 763728038*x^6 - 20263824*x^5 + 1518346*x^4 - 32160*x^3 + 1843*x^2 - 24*x + 1
Time: 41.110
> time LPolynomial(ChangeRing(Y, GF(601)));
78410163603001*x^10 + 782796974406*x^9 + 202537320333*x^8 + 5990157384*x^7 + 494754018*x^6 + 7825316*x^5 + 823218*x^4 + 16584*x^3 + 933*x^2 + 6*x + 1
Time: 49.020
> time LPolynomial(ChangeRing(Y, GF(1009)));
1045817322864049*x^10 - 26948711986586*x^9 + 1229610743613*x^8 - 4976379928*x^7 + 337153314*x^6 - 3755804*x^5 + 334146*x^4 - 4888*x^3 + 1197*x^2 - 26*x + 1
Time: 90.150
```

