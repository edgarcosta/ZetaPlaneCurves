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
