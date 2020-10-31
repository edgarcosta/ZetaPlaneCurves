AttachSpec("spec");


// returns the LPolynomial of f(x,y,z)=0 wth f homogeneous in GF(p)[x,y,z]
function ZetaPlaneCurve(f)
  return LPolynomial(Curve(ProjectiveSpace(Parent(f)), f));
end function;

function ZetaPlaneCurve2(f, p)
  return ZetaPlaneCurve(ChangeRing(f, GF(p)));
end function;


R<a,b,c> := PolynomialRing(Integers(), 3);
print "Testing Genus 5 curve defined by a degree 6 polynomial over GF(107)";
f := a^5*b + 10*a^4*b^2 + 105*a^3*b^3 + 62*a^2*b^4 + 66*a*b^5 + 65*b^6 + 45*a^5*c + 47*a^4*b*c + 41*a^3*b^2*c + 72*a^2*b^3*c + 99*a*b^4*c + 88*b^5*c + 78*a^4*c^2 + 29*a^3*b*c^2 + 53*a^2*b^2*c^2 + 9*a*b^3*c^2 + 8*b^4*c^2 + a^3*c^3 + 84*a^2*b*c^3 + 7*a*b^2*c^3 + 60*b^3*c^3 + 20*a^2*c^4 + 99*a*b*c^4 + 104*b^2*c^4 + 65*a*c^5 + 14*b*c^5 + 62*c^6;
time assert Coefficients(LPolynomial(f, 107)) eq [1, 4, 247, 560, 23642, 34648, 2529694, 6411440, 302585621, 524318404, 14025517307];
print "Passed!\n";

print "Testing Genus 5 curve defined by a degree 5 polynomial over GF(107)";
g := 13*a^5 + 3*a^4*b + 11*a^4*c + 11*a^3*b^2 + 14*a^3*b*c + 12*a^3*c^2 + 14*a^2*b^3 + 12*a^2*b^2*c + 4*a^2*b*c^2 + 8*a^2*c^3 + 14*a*b^4 + 16*a*b^3*c + a*b^2*c^2 + 2*a*b*c^3 + 10*a*c^4 + 9*b^5 + 3*b^4*c + 6*b^3*c^2 + 6*b^2*c^3 + 4*b*c^4 + 12*c^5;
time assert Coefficients(LPolynomial(g, 17)) eq [1, 4, 17, 72, 393, 2774, 6681, 20808, 83521, 334084, 1419857 ];
print "Passed!\n";


print "Testing different models for genus 6 curves over GF(7)";
p:=7;
print "Example 1";
f := 4*a^6 + 6*a^4*b^2 + 6*a^3*b^3 + 5*b^6 + 3*a^5*c + 3*a^4*b*c + 6*a^3*b^2*c + 6*a*b^4*c + 2*b^5*c + 3*a^4*c^2 + 3*a^3*b*c^2 + a^2*b^2*c^2 + 5*a*b^3*c^2 + 4*b^4*c^2 + a^3*c^3 + 2*a^2*b*c^3 + 3*a*b^2*c^3 + 2*b^3*c^3 + 5*a^2*c^4 + 2*a*b*c^4 + 4*b^2*c^4 + 3*a*c^5 + 2*b*c^5 + c^6;
g := a^10 + 4*a^9*b + a^8*b^2 + 5*a^7*b^3 + 6*a^6*b^4 + 4*a^5*b^5 + 6*a^4*b^6 + 6*a^2*b^8 + 4*a*b^9 + 6*b^10 + 3*a^8*b*c + 4*a^7*b^2*c + 3*a^6*b^3*c + 3*a^5*b^4*c + 4*a^4*b^5*c + 4*a^3*b^6*c + a^2*b^7*c + 3*a*b^8*c + 4*b^9*c + 6*a^8*c^2 + 2*a^7*b*c^2 + 6*a^6*b^2*c^2 + 3*a^5*b^3*c^2 + 5*a^4*b^4*c^2 + 5*a^3*b^5*c^2 + 2*a^2*b^6*c^2 + 4*a*b^7*c^2 + a^7*c^3 + 6*a^6*b*c^3 + 3*a^4*b^3*c^3 + 2*a^3*b^4*c^3 + 6*a^2*b^5*c^3 + 6*a*b^6*c^3 + 4*a^5*b*c^4 + a^2*b^4*c^4 + 6*a*b^5*c^4 + 5*b^6*c^4 + 4*a^5*c^5 + 4*a^4*b*c^5 + 4*a^3*b^2*c^5 + 6*a^2*b^3*c^5 + 4*b^5*c^5 + 5*a^4*c^6 + a^3*b*c^6 + 4*a^2*b^2*c^6 + a*b^3*c^6 + 3*b^4*c^6 + a^3*c^7 + 6*a^2*b*c^7 + 6*a*b^2*c^7 + 4*b^3*c^7 + 6*a^2*c^8 + 2*a*b*c^8 + b^2*c^8 + a*c^9 + b*c^9 + 4*c^10;
h := a^10 + 6*a^9*b + 5*a^8*b^2 + 5*a^7*b^3 + a^6*b^4 + 4*a^5*b^5 + a^4*b^6 + a^3*b^7 + 4*a^2*b^8 + 2*a*b^9 + 4*b^10 + 4*a^9*c + 3*a^8*b*c + 4*a^7*b^2*c + 3*a^6*b^3*c + 2*a^4*b^5*c + 3*a^3*b^6*c + a*b^8*c + 5*b^9*c + 5*a^8*c^2 + 6*a^7*b*c^2 + 3*a^6*b^2*c^2 + 2*a^5*b^3*c^2 + 4*a^3*b^5*c^2 + a^2*b^6*c^2 + 2*a*b^7*c^2 + 4*b^8*c^2 + 6*a^7*c^3 + 2*a^6*b*c^3 + 5*a^5*b^2*c^3 + 4*a^4*b^3*c^3 + a^2*b^5*c^3 + 5*a*b^6*c^3 + b^7*c^3 + 4*a^6*c^4 + a^5*b*c^4 + a^4*b^2*c^4 + 2*a^3*b^3*c^4 + 2*a^2*b^4*c^4 + 5*a^4*b*c^5 + 4*a^3*b^2*c^5 + 4*a^2*b^3*c^5 + 5*b^5*c^5 + 6*a^3*b*c^6 + 3*a^2*b^2*c^6 + 4*a*b^3*c^6 + 4*b^4*c^6 + 6*a^3*c^7 + 6*a^2*b*c^7 + 2*b^3*c^7 + 5*a^2*c^8 + 3*a*b*c^8 + 2*b^2*c^8 + 4*a*c^9 + 3*b*c^9 + 6*c^10;
printf "nodal model degree 6: ";
time L1 := LPolynomial(f, p);
printf "non nodal degree 10: ";
time L2 := LPolynomial(g, p : NewModelTries:=0);
assert L1 eq L2;
printf "non nodal degree 10, find new model: ";
time L2 := LPolynomial(g, p);
assert L1 eq L2;
printf "nodal model degree 10: ";
time L2 := LPolynomial(h, p);
assert L1 eq L2;
printf "Magma with nodal model degree 6: ";
time L2 := ZetaPlaneCurve2(f, p);
assert L1 eq L2;
printf "Magma with non nodal model degree 10: ";
time L2 := ZetaPlaneCurve2(g, p);
assert L1 eq L2;
printf "Magma with nodal model degree 10: ";
time L2 := ZetaPlaneCurve2(h, p);
assert L1 eq L2;
print "Passed!\n";


print "Example 2";
f := 3*a^5*b + 3*a^4*b^2 + 2*a^3*b^3 + 6*a^2*b^4 + a*b^5 + 3*b^6 + 2*a^4*b*c + 6*a^3*b^2*c + 5*a^2*b^3*c + 4*a*b^4*c + 4*b^5*c + a^4*c^2 + 6*a^3*b*c^2 + 4*a^2*b^2*c^2 + a*b^3*c^2 + b^4*c^2 + 5*a^2*b*c^3 + b^3*c^3 + 3*a^2*c^4 + 6*a*b*c^4 + b^2*c^4 + 2*b*c^5 + 4*c^6;
g := a^9*b + 3*a^8*b^2 + 4*a^7*b^3 + 4*a^6*b^4 + 4*a^4*b^6 + 6*a^3*b^7 + 2*a^2*b^8 + 5*b^10 + 3*a^8*b*c + 5*a^7*b^2*c + 2*a^6*b^3*c + 2*a^5*b^4*c + 5*a^4*b^5*c + 3*a^3*b^6*c + 5*a^2*b^7*c + a*b^8*c + 5*b^9*c + 2*a^8*c^2 + 4*a^7*b*c^2 + a^6*b^2*c^2 + 3*a^5*b^3*c^2 + 3*a^4*b^4*c^2 + 4*a^3*b^5*c^2 + 6*a^2*b^6*c^2 + 4*a*b^7*c^2 + 2*a^7*c^3 + a^6*b*c^3 + 3*a^5*b^2*c^3 + 5*a^4*b^3*c^3 + 3*a^2*b^5*c^3 + 6*a*b^6*c^3 + 5*b^7*c^3 + 3*a^6*c^4 + 5*a^5*b*c^4 + 5*a^4*b^2*c^4 + 3*a^2*b^4*c^4 + 3*a*b^5*c^4 + 5*b^6*c^4 + 6*a^5*c^5 + 6*a^3*b^2*c^5 + 2*a^2*b^3*c^5 + 3*a*b^4*c^5 + 6*b^5*c^5 + a^4*c^6 + 2*a^3*b*c^6 + 6*a^2*b^2*c^6 + 4*a*b^3*c^6 + 6*b^4*c^6 + 5*a^3*c^7 + 5*a^2*b*c^7 + 4*a*b^2*c^7 + 5*b^3*c^7 + 3*a^2*c^8 + a*b*c^8 + b^2*c^8 + 6*b*c^9 + 3*c^10;
h := a^9*b + 5*a^8*b^2 + 5*a^7*b^3 + 6*a^6*b^4 + 4*a^5*b^5 + 5*a^4*b^6 + 4*a^3*b^7 + 4*a^2*b^8 + 3*a*b^9 + 2*b^10 + 5*a^8*b*c + a^7*b^2*c + 2*a^6*b^3*c + 4*a^5*b^4*c + 4*a^4*b^5*c + 4*a^3*b^6*c + a^2*b^7*c + 3*a*b^8*c + 5*b^9*c + a^8*c^2 + 2*a^7*b*c^2 + 5*a^6*b^2*c^2 + 4*a^5*b^3*c^2 + 3*a^3*b^5*c^2 + 6*a^2*b^6*c^2 + 6*a*b^7*c^2 + b^8*c^2 + 5*a^7*c^3 + 3*a^6*b*c^3 + 5*a^5*b^2*c^3 + 6*a^4*b^3*c^3 + 3*a^3*b^4*c^3 + 3*a*b^6*c^3 + b^7*c^3 + 2*a^6*c^4 + 6*a^5*b*c^4 + 3*a^4*b^2*c^4 + 6*a^2*b^4*c^4 + 6*a*b^5*c^4 + b^6*c^4 + 5*a^5*c^5 + 4*a^4*b*c^5 + 2*a^3*b^2*c^5 + 4*a^2*b^3*c^5 + 2*a*b^4*c^5 + 4*a^3*b*c^6 + a^3*c^7 + 5*a^2*b*c^7 + 4*b^3*c^7 + 6*a^2*c^8 + 6*a*b*c^8 + 6*b^2*c^8 + 4*a*c^9 + 3*b*c^9 + 5*c^10;
printf "nodal model degree 6: ";
time L1 := LPolynomial(f, p);
printf "non nodal degree 10: ";
time L2 := LPolynomial(g, p : NewModelTries:=0);
assert L1 eq L2;
printf "non nodal degree 10, find new model: ";
time L2 := LPolynomial(g, p);
assert L1 eq L2;
printf "nodal model degree 10: ";
time L2 := LPolynomial(h, p);
assert L1 eq L2;
printf "Magma with nodal model degree 6: ";
time L2 := ZetaPlaneCurve2(f, p);
assert L1 eq L2;
printf "Magma with non nodal model degree 10: ";
time L2 := ZetaPlaneCurve2(g, p);
assert L1 eq L2;
printf "Magma with nodal model degree 10: ";
time L2 := ZetaPlaneCurve2(h, p);
assert L1 eq L2;
print "Passed!\n";
print "Testing different models for genus 6 curves over GF(13)";
p := 13;
f := 6*a^6 + 5*a^5*b + 11*a^4*b^2 + 6*a^3*b^3 + 11*a^2*b^4 + 2*b^6 + 6*a^5*c + 12*a^4*b*c + a^3*b^2*c + 7*a^2*b^3*c + 11*a*b^4*c + 11*a^4*c^2 + 12*a^3*b*c^2 + 8*a^2*b^2*c^2 + 12*a*b^3*c^2 + 4*b^4*c^2 + 2*a^3*c^3 + 10*a^2*b*c^3 + 8*a*b^2*c^3 + 8*b^3*c^3 + 11*a^2*c^4 + 7*a*b*c^4 + 4*b^2*c^4 + 7*a*c^5 + c^6;
g := a^10 + 4*a^9*b + 2*a^8*b^2 + 11*a^7*b^3 + 3*a^6*b^4 + 10*a^5*b^5 + 3*a^4*b^6 + 7*a^3*b^7 + 4*a^2*b^8 + 11*a*b^9 + 6*a^9*c + 7*a^8*b*c + 9*a^7*b^2*c + 8*a^6*b^3*c + 11*a^5*b^4*c + 3*a^4*b^5*c + 4*a^3*b^6*c + 11*a^2*b^7*c + 10*a*b^8*c + 9*b^9*c + 12*a^8*c^2 + 11*a^7*b*c^2 + 3*a^6*b^2*c^2 + 8*a^5*b^3*c^2 + 4*a^4*b^4*c^2 + 3*a^3*b^5*c^2 + 12*a^2*b^6*c^2 + 4*a*b^7*c^2 + 11*b^8*c^2 + 10*a^7*c^3 + 8*a^6*b*c^3 + 5*a^5*b^2*c^3 + 6*a^4*b^3*c^3 + 6*a^3*b^4*c^3 + 7*a^2*b^5*c^3 + a*b^6*c^3 + b^7*c^3 + 6*a^6*c^4 + 9*a^4*b^2*c^4 + 4*a^3*b^3*c^4 + 9*a^2*b^4*c^4 + 10*a*b^5*c^4 + 6*b^6*c^4 + 8*a^4*b*c^5 + 11*a^3*b^2*c^5 + 10*a^2*b^3*c^5 + 8*a*b^4*c^5 + 6*a^4*c^6 + 9*a^3*b*c^6 + 4*a^2*b^2*c^6 + 9*a*b^3*c^6 + b^4*c^6 + 12*a^3*c^7 + 10*a^2*b*c^7 + 9*a*b^2*c^7 + 7*b^3*c^7 + 8*a^2*c^8 + 9*a*b*c^8 + 11*b^2*c^8 + 9*a*c^9 + 9*b*c^9 + 5*c^10;
h := a^10 + 4*a^9*b + 11*a^8*b^2 + 3*a^7*b^3 + 10*a^6*b^4 + a^5*b^5 + 4*a^3*b^7 + 5*a^2*b^8 + 9*a*b^9 + 7*b^10 + 3*a^9*c + 10*a^8*b*c + 3*a^7*b^2*c + 4*a^6*b^3*c + 5*a^5*b^4*c + 9*a^4*b^5*c + 11*a^3*b^6*c + 2*a^2*b^7*c + 7*a*b^8*c + 8*b^9*c + 12*a^8*c^2 + 11*a^7*b*c^2 + 8*a^6*b^2*c^2 + 5*a^5*b^3*c^2 + 10*a^3*b^5*c^2 + 2*a^2*b^6*c^2 + 7*a*b^7*c^2 + 2*b^8*c^2 + 11*a^7*c^3 + a^6*b*c^3 + 3*a^5*b^2*c^3 + 8*a^4*b^3*c^3 + 10*a^3*b^4*c^3 + 9*a^2*b^5*c^3 + 11*a*b^6*c^3 + 3*b^7*c^3 + 6*a^5*b*c^4 + 8*a^4*b^2*c^4 + 3*a^3*b^3*c^4 + 2*a^2*b^4*c^4 + 3*a*b^5*c^4 + 2*a^4*b*c^5 + 8*a^3*b^2*c^5 + 9*a^2*b^3*c^5 + 10*a*b^4*c^5 + 6*b^5*c^5 + a^4*c^6 + 8*a^3*b*c^6 + 3*a^2*b^2*c^6 + 11*a*b^3*c^6 + 7*b^4*c^6 + 5*a^3*c^7 + 6*a^2*b*c^7 + 2*b^3*c^7 + 4*a^2*c^8 + 12*a*b*c^8 + 12*b^2*c^8 + 11*a*c^9 + b*c^9 + 9*c^10;
printf "nodal model degree 6: ";
time L1 := LPolynomial(f, p);
printf "non nodal degree 10: ";
time L2 := LPolynomial(g, p : NewModelTries:=0);
assert L1 eq L2;
printf "non nodal degree 10, find new model: ";
time L2 := LPolynomial(g, p);
assert L1 eq L2;
printf "nodal model degree 10: ";
time L2 := LPolynomial(h, p);
assert L1 eq L2;
// 60.260s
printf "Magma with nodal model degree 6: ";
time L2 := ZetaPlaneCurve2(f, p);
assert L1 eq L2;
p := 13;
// too slow to test, doesn't add much
//printf "Magma with non nodal model degree 10: ";
//time L2 := ZetaPlaneCurve2(g, p); // 146.950s
//assert L1 eq L2;
//printf "Magma with nodal model degree 10: ";
//time L2 := ZetaPlaneCurve2(h, p); // 149.580
//assert L1 eq L2;

f := 8*a^6 + 3*a^5*b + 9*a^4*b^2 + 6*a^3*b^3 + 7*a^2*b^4 + a*b^5 + 5*b^6 + 7*a^5*c + 5*a^4*b*c + 8*a^3*b^2*c + 2*a^2*b^3*c + 4*a*b^4*c + 6*a^4*c^2 + 11*a^3*b*c^2 + 11*a^2*b^2*c^2 + 11*b^4*c^2 + 10*a^3*c^3 + 10*a^2*b*c^3 + 9*a*b^2*c^3 + 2*b^3*c^3 + 8*a^2*c^4 + 2*a*b*c^4 + 10*b^2*c^4 + 4*a*c^5 + b*c^5 + 11*c^6;
g := a^12 + 4*a^11*b + 8*a^10*b^2 + 9*a^9*b^3 + 12*a^8*b^4 + 7*a^7*b^5 + 12*a^5*b^7 + 9*a^4*b^8 + a^3*b^9 + 11*a^2*b^10 + 2*a*b^11 + 3*b^12 + 8*a^11*c + 7*a^10*b*c + 3*a^9*b^2*c + 10*a^8*b^3*c + 9*a^7*b^4*c + 4*a^6*b^5*c + 11*a^4*b^7*c + a^3*b^8*c + 6*a^2*b^9*c + 5*a*b^10*c + 5*b^11*c + 5*a^10*c^2 + 5*a^9*b*c^2 + 4*a^8*b^2*c^2 + 4*a^7*b^3*c^2 + 3*a^6*b^4*c^2 + 12*a^5*b^5*c^2 + 11*a^4*b^6*c^2 + a^3*b^7*c^2 + 12*a^2*b^8*c^2 + 10*a*b^9*c^2 + 7*b^10*c^2 + 6*a^9*c^3 + 10*a^8*b*c^3 + 3*a^7*b^2*c^3 + a^6*b^3*c^3 + 6*a^5*b^4*c^3 + 11*a^4*b^5*c^3 + 11*a^3*b^6*c^3 + 12*a^2*b^7*c^3 + 9*a*b^8*c^3 + 2*b^9*c^3 + 6*a^8*c^4 + a^7*b*c^4 + 4*a^6*b^2*c^4 + 11*a^5*b^3*c^4 + a^4*b^4*c^4 + 12*a^3*b^5*c^4 + 12*a^2*b^6*c^4 + 6*a*b^7*c^4 + 3*b^8*c^4 + 11*a^7*c^5 + 7*a^6*b*c^5 + 9*a^4*b^3*c^5 + 12*a^3*b^4*c^5 + 4*a^2*b^5*c^5 + 7*a*b^6*c^5 + 2*b^7*c^5 + a^6*c^6 + 8*a^5*b*c^6 + 5*a^4*b^2*c^6 + 9*a^3*b^3*c^6 + 3*a^2*b^4*c^6 + 11*a*b^5*c^6 + 12*b^6*c^6;
h := a^12 + 5*a^11*b + 12*a^10*b^2 + 2*a^9*b^3 + 12*a^8*b^4 + 8*a^7*b^5 + 6*a^6*b^6 + 10*a^5*b^7 + 8*a^4*b^8 + 7*a^3*b^9 + 4*a^2*b^10 + 2*a*b^11 + 4*b^12 + 8*a^11*c + 12*a^10*b*c + a^9*b^2*c + 10*a^7*b^4*c + 3*a^6*b^5*c + 9*a^5*b^6*c + 2*a^4*b^7*c + a^3*b^8*c + 10*a^2*b^9*c + 11*a*b^10*c + 11*b^11*c + 3*a^10*c^2 + 10*a^9*b*c^2 + 3*a^8*b^2*c^2 + a^7*b^3*c^2 + 7*a^6*b^4*c^2 + 5*a^5*b^5*c^2 + a^4*b^6*c^2 + 8*a^3*b^7*c^2 + 5*a^2*b^8*c^2 + 9*a*b^9*c^2 + 8*b^10*c^2 + 3*a^9*c^3 + 12*a^8*b*c^3 + 4*a^7*b^2*c^3 + 3*a^6*b^3*c^3 + 6*a^5*b^4*c^3 + 7*a^4*b^5*c^3 + 11*a^3*b^6*c^3 + 9*a^2*b^7*c^3 + 12*a*b^8*c^3 + 6*a^8*c^4 + 4*a^7*b*c^4 + 12*a^6*b^2*c^4 + 8*a^5*b^3*c^4 + 10*a^4*b^4*c^4 + 11*a^3*b^5*c^4 + 4*a^2*b^6*c^4 + 7*a*b^7*c^4 + 3*a^7*c^5 + 7*a^6*b*c^5 + 10*a^5*b^2*c^5 + 4*a^4*b^3*c^5 + 6*a^3*b^4*c^5 + 5*a^2*b^5*c^5 + a*b^6*c^5 + 9*b^7*c^5 + 3*a^6*c^6 + 12*a^5*b*c^6 + 12*a^4*b^2*c^6 + 12*a^3*b^3*c^6 + 7*a^2*b^4*c^6 + 2*a*b^5*c^6 + 11*b^6*c^6 + 8*a^4*b*c^7 + 4*a^3*b^2*c^7 + 5*a^2*b^3*c^7 + 3*a*b^4*c^7 + 12*b^5*c^7 + 3*a^4*c^8 + 2*a^3*b*c^8 + 11*a^2*b^2*c^8 + a*b^3*c^8 + 9*b^4*c^8 + 11*a^3*c^9 + a^2*b*c^9 + 11*a*b^2*c^9 + 7*b^3*c^9 + 8*a^2*c^10 + 11*a*b*c^10 + 3*b^2*c^10 + 11*a*c^11 + 5*b*c^11 + 5*c^12;
printf "nodal model degree 6: ";
time L1 := LPolynomial(f, p);
printf "non nodal degree 12: ";
time L2 := LPolynomial(g, p : NewModelTries:=0);
assert L1 eq L2;
printf "non nodal degree 12, find new model: ";
time L2 := LPolynomial(g, p);
assert L1 eq L2;
printf "nodal model degree 12: ";
time L2 := LPolynomial(h, p);
assert L1 eq L2;
// too slow to compare against magma


print "Passed!\n";

print "Testing different models for genus 6 curves over GF(17)";
p := 17;
f := 13*a^6 + a^5*b + 4*a^4*b^2 + 3*a^2*b^4 + 8*a*b^5 + b^6 + a^5*c + 5*a^4*b*c + 9*a^3*b^2*c + 7*a^2*b^3*c + 6*a*b^4*c + 8*b^5*c + 4*a^4*c^2 + 10*a^3*b*c^2 + 9*a^2*b^2*c^2 + 2*a*b^3*c^2 + 11*b^4*c^2 + 4*a^3*c^3 + 3*a^2*b*c^3 + 4*a*b^2*c^3 + b^3*c^3 + 11*a^2*c^4 + 14*a*b*c^4 + 16*b^2*c^4 + 11*a*c^5 + 5*b*c^5 + 16*c^6;
g := a^10 + 13*a^9*b + 2*a^8*b^2 + 14*a^7*b^3 + 7*a^6*b^4 + 2*a^5*b^5 + 10*a^4*b^6 + 15*a^3*b^7 + 9*a^2*b^8 + 11*a*b^9 + 16*b^10 + 10*a^9*c + 2*a^8*b*c + 6*a^7*b^2*c + 16*a^6*b^3*c + 16*a^5*b^4*c + 15*a^4*b^5*c + 11*a^3*b^6*c + 6*a^2*b^7*c + 5*a*b^8*c + 3*b^9*c + 10*a^8*c^2 + 13*a^7*b*c^2 + 5*a^6*b^2*c^2 + 16*a^5*b^3*c^2 + 15*a^4*b^4*c^2 + 4*a^3*b^5*c^2 + 12*a^2*b^6*c^2 + 14*a*b^7*c^2 + 5*b^8*c^2 + 13*a^6*b*c^3 + 5*a^5*b^2*c^3 + 3*a^4*b^3*c^3 + 11*a^3*b^4*c^3 + 9*a^2*b^5*c^3 + 4*a*b^6*c^3 + 8*b^7*c^3 + 7*a^6*c^4 + 14*a^5*b*c^4 + 5*a^4*b^2*c^4 + 9*a^3*b^3*c^4 + 7*a^2*b^4*c^4 + 16*a*b^5*c^4 + 11*b^6*c^4 + 8*a^5*c^5 + 8*a^4*b*c^5 + 15*a^3*b^2*c^5 + 5*a^2*b^3*c^5 + 16*a*b^4*c^5 + 11*b^5*c^5 + 15*a^4*c^6 + 14*a^3*b*c^6 + 2*a^2*b^2*c^6 + 8*a*b^3*c^6 + 10*b^4*c^6 + 11*a^3*c^7 + 16*a^2*b*c^7 + 11*a*b^2*c^7 + b^3*c^7 + 5*a^2*c^8 + a*b*c^8 + 4*b^2*c^8 + 6*a*c^9 + 14*b*c^9 + c^10;
h := a^10 + 11*a^9*b + 8*a^8*b^2 + 13*a^7*b^3 + 3*a^6*b^4 + 13*a^5*b^5 + 14*a^3*b^7 + 6*a^2*b^8 + 14*a*b^9 + 9*b^10 + 14*a^9*c + 4*a^8*b*c + 11*a^7*b^2*c + 12*a^6*b^3*c + 4*a^5*b^4*c + 12*a^4*b^5*c + 10*a^3*b^6*c + 9*a^2*b^7*c + a*b^8*c + 13*b^9*c + 10*a^8*c^2 + 11*a^7*b*c^2 + 10*a^6*b^2*c^2 + 3*a^5*b^3*c^2 + 2*a^3*b^5*c^2 + 2*a^2*b^6*c^2 + 9*a*b^7*c^2 + 5*b^8*c^2 + 12*a^6*b*c^3 + 6*a^5*b^2*c^3 + 2*a^4*b^3*c^3 + 9*a^3*b^4*c^3 + 12*a^2*b^5*c^3 + 2*a*b^6*c^3 + 14*b^7*c^3 + 13*a^6*c^4 + 3*a^5*b*c^4 + 15*a^4*b^2*c^4 + 10*a^3*b^3*c^4 + 9*a^2*b^4*c^4 + 14*a*b^5*c^4 + 12*b^6*c^4 + 11*a^5*c^5 + 10*a^4*b*c^5 + 6*a^3*b^2*c^5 + 4*a^2*b^3*c^5 + 16*a*b^4*c^5 + b^5*c^5 + 16*a^4*c^6 + 2*a^2*b^2*c^6 + 14*a*b^3*c^6 + 5*b^4*c^6 + 5*a^3*c^7 + 13*a^2*b*c^7 + 10*a*b^2*c^7 + 15*b^3*c^7 + 4*a^2*c^8 + 2*a*b*c^8 + 13*b^2*c^8 + 10*a*c^9 + 15*b*c^9 + 5*c^10;
printf "nodal model degree 6: ";
time L1 := LPolynomial(f, p);
printf "non nodal degree 10: ";
time L2 := LPolynomial(g, p : NewModelTries:=0);
assert L1 eq L2;
printf "non nodal degree 10 find new model: ";
time L2 := LPolynomial(g, p);
assert L1 eq L2;
printf "nodal model degree 10: ";
time L2 := LPolynomial(h, p);
assert L1 eq L2;
print "Passed!\n";

print "Testing different models for genus 6 curves over GF(23)";
p := 23;
f := 14*a^6 + 14*a^5*b + 4*a^4*b^2 + 20*a^3*b^3 + 5*a^2*b^4 + 7*a*b^5 + 15*b^6 + 15*a^5*c + 16*a^4*b*c + 21*a^3*b^2*c + 20*a^2*b^3*c + 13*a*b^4*c + 8*b^5*c + 13*a^4*c^2 + 14*a^3*b*c^2 + 13*a^2*b^2*c^2 + 9*a*b^3*c^2 + 12*b^4*c^2 + 5*a^3*c^3 + 20*a^2*b*c^3 + 11*a*b^2*c^3 + 17*b^3*c^3 + 15*a^2*c^4 + 13*a*b*c^4 + 11*b^2*c^4 + 15*a*c^5 + 5*b*c^5 + 8*c^6;
g := a^10 + 12*a^9*b + 6*a^8*b^2 + 8*a^7*b^3 + 8*a^6*b^4 + 14*a^5*b^5 + 7*a^4*b^6 + 5*a^3*b^7 + 3*a^2*b^8 + 5*a*b^9 + 4*b^10 + 6*a^9*c + 2*a^8*b*c + a^7*b^2*c + 8*a^6*b^3*c + 16*a^5*b^4*c + a^4*b^5*c + 6*a^3*b^6*c + 9*a^2*b^7*c + 7*b^9*c + 19*a^8*c^2 + 16*a^7*b*c^2 + 21*a^6*b^2*c^2 + 6*a^5*b^3*c^2 + 20*a^4*b^4*c^2 + 4*a^3*b^5*c^2 + 19*a^2*b^6*c^2 + 16*a*b^7*c^2 + 22*b^8*c^2 + 11*a^7*c^3 + 10*a^6*b*c^3 + 5*a^5*b^2*c^3 + a^4*b^3*c^3 + 13*a^3*b^4*c^3 + 5*a^2*b^5*c^3 + 3*a*b^6*c^3 + 6*b^7*c^3 + 20*a^6*c^4 + 14*a^5*b*c^4 + 7*a^4*b^2*c^4 + 9*a^3*b^3*c^4 + 9*a^2*b^4*c^4 + 10*a*b^5*c^4 + 19*b^6*c^4 + 14*a^5*c^5 + a^4*b*c^5 + 9*a^3*b^2*c^5 + 2*a^2*b^3*c^5 + 13*a*b^4*c^5 + 6*b^5*c^5 + 6*a^4*c^6 + 8*a^3*b*c^6 + 10*a^2*b^2*c^6 + 4*a*b^3*c^6 + 13*a^3*c^7 + 14*a^2*b*c^7 + 16*a*b^2*c^7 + 21*b^3*c^7 + 16*a^2*c^8 + 21*a*b*c^8 + 13*b^2*c^8 + 13*a*c^9 + 3*b*c^9 + 22*c^10;
h := a^10 + 22*a^8*b^2 + 12*a^7*b^3 + 14*a^6*b^4 + 21*a^5*b^5 + 20*a^4*b^6 + 17*a^3*b^7 + 8*a^2*b^8 + 19*a*b^9 + 6*b^10 + 2*a^9*c + 16*a^8*b*c + 22*a^7*b^2*c + 16*a^6*b^3*c + 22*a^5*b^4*c + 12*a^4*b^5*c + 11*a^3*b^6*c + 10*a^2*b^7*c + 2*a*b^8*c + 18*b^9*c + 9*a^8*c^2 + 7*a^7*b*c^2 + 9*a^6*b^2*c^2 + 7*a^5*b^3*c^2 + 21*a^4*b^4*c^2 + 15*a^3*b^5*c^2 + 4*a^2*b^6*c^2 + 15*a*b^7*c^2 + 21*b^8*c^2 + a^6*b*c^3 + 2*a^5*b^2*c^3 + a^4*b^3*c^3 + 3*a^3*b^4*c^3 + 16*a^2*b^5*c^3 + 22*a*b^6*c^3 + 7*b^7*c^3 + 8*a^6*c^4 + 5*a^5*b*c^4 + a^4*b^2*c^4 + 3*a^3*b^3*c^4 + 19*a^2*b^4*c^4 + 14*a*b^5*c^4 + 21*b^6*c^4 + 18*a^5*c^5 + 2*a^4*b*c^5 + 22*a^3*b^2*c^5 + 9*a^2*b^3*c^5 + 7*a*b^4*c^5 + 21*b^5*c^5 + 22*a^4*c^6 + 22*a^3*b*c^6 + 19*a^2*b^2*c^6 + 6*a*b^3*c^6 + 9*a^3*c^7 + 20*a^2*b*c^7 + 19*a*b^2*c^7 + 22*b^3*c^7 + 22*a^2*c^8 + 11*a*b*c^8 + 20*b^2*c^8 + 4*a*c^9 + 4*c^10;
printf "nodal model degree 6: ";
time L1 := LPolynomial(f, p);
printf "non nodal degree 10: ";
time L2 := LPolynomial(g, p : NewModelTries:=0);
assert L1 eq L2;
printf "non nodal degree 10, find new model: ";
time L2 := LPolynomial(g, p);
assert L1 eq L2;
printf "nodal model degree 10: ";
time L2 := LPolynomial(h, p);
assert L1 eq L2;
print "Passed!\n";

print "Testing different models for genus 6 curves over GF(107)";
p := 107;
f := 76*a^6 + 62*a^5*b + 17*a^4*b^2 + 66*a^3*b^3 + 7*a^2*b^4 + 100*a*b^5 + 55*b^6 + 83*a^5*c + 93*a^4*b*c + 80*a^3*b^2*c + 62*a^2*b^3*c + 87*a*b^4*c + 38*b^5*c + 10*a^4*c^2 + 64*a^3*b*c^2 + 23*a^2*b^2*c^2 + 64*a*b^3*c^2 + 25*b^4*c^2 + 24*a^3*c^3 + a^2*b*c^3 + 27*a*b^2*c^3 + 91*b^3*c^3 + 24*a^2*c^4 + 71*a*b*c^4 + 6*b^2*c^4 + 69*a*c^5 + 83*b*c^5 + 88*c^6;
g := a^10 + 63*a^9*b + 36*a^8*b^2 + 50*a^7*b^3 + 47*a^6*b^4 + 17*a^4*b^6 + 88*a^3*b^7 + 63*a^2*b^8 + 82*a*b^9 + 94*b^10 + 58*a^9*c + 98*a^8*b*c + 84*a^7*b^2*c + 19*a^6*b^3*c + 49*a^5*b^4*c + 101*a^4*b^5*c + 71*a^3*b^6*c + 65*a^2*b^7*c + 7*a*b^8*c + 75*b^9*c + 101*a^8*c^2 + 65*a^7*b*c^2 + 39*a^6*b^2*c^2 + 103*a^5*b^3*c^2 + 15*a^4*b^4*c^2 + 98*a^3*b^5*c^2 + 39*a^2*b^6*c^2 + 22*a*b^7*c^2 + 56*b^8*c^2 + 32*a^7*c^3 + 45*a^6*b*c^3 + 76*a^5*b^2*c^3 + 22*a^4*b^3*c^3 + 90*a^3*b^4*c^3 + 46*a^2*b^5*c^3 + 15*a*b^6*c^3 + 7*b^7*c^3 + 35*a^6*c^4 + 97*a^5*b*c^4 + 85*a^4*b^2*c^4 + 99*a^3*b^3*c^4 + 45*a^2*b^4*c^4 + 20*a*b^5*c^4 + 49*b^6*c^4 + 35*a^5*c^5 + 83*a^4*b*c^5 + 53*a^3*b^2*c^5 + 64*a^2*b^3*c^5 + 42*a*b^4*c^5 + 59*b^5*c^5 + 104*a^4*c^6 + 52*a^3*b*c^6 + 46*a^2*b^2*c^6 + 57*a*b^3*c^6 + 100*b^4*c^6 + 99*a^3*c^7 + 92*a^2*b*c^7 + 104*a*b^2*c^7 + 60*b^3*c^7 + 106*a^2*c^8 + 74*a*b*c^8 + 4*b^2*c^8 + 66*a*c^9 + 5*b*c^9 + 22*c^10;
h := a^10 + 40*a^9*b + 48*a^8*b^2 + 16*a^7*b^3 + 37*a^6*b^4 + 44*a^5*b^5 + 11*a^4*b^6 + 34*a^3*b^7 + 82*a^2*b^8 + 54*a*b^9 + 80*b^10 + 67*a^9*c + 77*a^8*b*c + 75*a^7*b^2*c + a^6*b^3*c + 93*a^5*b^4*c + 13*a^4*b^5*c + 52*a^3*b^6*c + 44*a^2*b^7*c + 51*a*b^8*c + 69*b^9*c + 42*a^8*c^2 + 86*a^7*b*c^2 + 4*a^6*b^2*c^2 + 13*a^5*b^3*c^2 + 64*a^4*b^4*c^2 + 74*a^3*b^5*c^2 + 76*a^2*b^6*c^2 + 97*a*b^7*c^2 + 77*b^8*c^2 + 102*a^7*c^3 + 44*a^6*b*c^3 + 17*a^5*b^2*c^3 + 7*a^4*b^3*c^3 + 9*a^3*b^4*c^3 + 42*a^2*b^5*c^3 + 31*a*b^6*c^3 + 61*b^7*c^3 + 89*a^6*c^4 + a^5*b*c^4 + 42*a^4*b^2*c^4 + 53*a^3*b^3*c^4 + 56*a^2*b^4*c^4 + 24*a*b^5*c^4 + 17*b^6*c^4 + 103*a^5*c^5 + 95*a^4*b*c^5 + 49*a^3*b^2*c^5 + 51*a^2*b^3*c^5 + 10*a*b^4*c^5 + 95*b^5*c^5 + 101*a^4*c^6 + 47*a^3*b*c^6 + 106*a^2*b^2*c^6 + 102*a*b^3*c^6 + 105*b^4*c^6 + 65*a^3*c^7 + 77*a^2*b*c^7 + 42*a*b^2*c^7 + 81*b^3*c^7 + 16*a^2*c^8 + 31*a*b*c^8 + 12*b^2*c^8 + 23*a*c^9 + 61*b*c^9 + 45*c^10;
printf "nodal model degree 6: ";
time L1 := LPolynomial(f, p);
printf "non nodal degree 12: ";
time L2 := LPolynomial(g, p : NewModelTries:=0);
assert L1 eq L2;
printf "non nodal degree 10, find new model: ";
time L2 := LPolynomial(g, p);
assert L1 eq L2;
printf "nodal model degree 10: ";
time L2 := LPolynomial(h, p);
assert L1 eq L2;
print "Passed!\n";





print "Testing genus 5 curve that is a quartic cover of the elliptic curve 330e2";
h := 900*a^6 - 8236*a^4*c^2 + 784*a^3*b^2*c + 11*a^2*b^4 + 1228*a^2*c^4 - 16*a*b^2*c^3 - 11*b^4*c^2 - 36*c^6;
E := EllipticCurve("330e2");
for p in [19, 107, 211] do
  printf "p=%o\n", p;
  time L1 := LPolynomial(h, p : KnownFactor := EulerFactor(E, p));
  L2 := LPolynomial(h, p);
  try
    t0 := Cputime();
    L3 := ZetaPlaneCurve2(h, p);
    print "Magma time:", Cputime() - t0;
  catch E
    print "Magma failed to compute L-polynomial";
    L3 := 1;
  end try;

  assert L1 eq L2;
  assert L1 eq L3 or L3 eq 1;
end for;
print "Passed!\n";


print "Testing against random curves and small p.\nVery long test!";
global_errors := [* *];
for p in PrimesInInterval(7, 20) do
    P<x,y,z> := ProjectivePlane(GF(p));
    for d in [3..6] do
        maxg := Min((d-1)*(d-2) div 2, 6);
        if p ge 17 then
          maxg := Min(maxg, 5);
        end if;
        for g in [1..maxg] do
            e := Ceiling(Log(p,4*g*p^(g/2)));
            if p gt e then
                ct := 0;
                t0 := Cputime();
                errors := [* *];
                tt1 := 0;
                tt2 := 0;
                while ct lt 10 do
                    C := RandomNodalCurve(d, g, P);
                    f := DefiningEquation(C);
                    t1 := Cputime();
                    L1 := LPolynomial(f);
                    t1 := Cputime() - t1;
                    try
                        t2 := Cputime();
                        L2 := LPolynomial(C);
                        tt2 +:= Cputime() - t2;
                        tt1 +:= t1;
                        ct +:= 1;
                    catch E
                        Append(~errors, <p, f>);
                        L2 := 1;
                    end try;
                    assert L2 eq 1 or L2 eq L1;
                    if (ct gt 0 and (ct + 1)*(Cputime() - t0)/ct gt 10) or Cputime() - t0 gt 10 then
                        break;
                    end if;
                end while;
                printf "p=%o d=%o g=%o : tested %o curves in %os (%os vs Magma %os)\n", p, d, g, ct, Cputime() - t0, tt1, tt2;
                if #errors gt 0 then
                    printf "in this set Magma failed %o times to compute the L-polynomial\n", #errors;
                    global_errors := global_errors cat errors;
                end if;
            end if;
        end for;
    end for;
end for;
print "Passed!\n";
