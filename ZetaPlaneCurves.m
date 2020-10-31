freeze;
/***********************************************************************************************

  ZetaPlaneCurves.

  Computes LPolynomial of plane curves, potentially singular.

  Code is based on Theorem 3.1 of D. Harvey's "Computing zeta functions of arithmetic schemes"

  intrinsics in this file:

    - LPolynomial(f::RngMPolElt : KnownFactor:=false, corrections:=false) -> RngUPolElt
    - LPolynomial(f::RngMPolElt, p::RngIntElt : KnownFactor:=false, corrections:=false) -> RngUPolElt


  Edgar Costa, 2020
***********************************************************************************************/

declare verbose ZetaPlaneCurve, 2;



// Given f defining a plane curve in P2
// computes the diferences of point counts between the plane model and its normalisation for GF(p^r)
// at the moment we only support nodal singularities
function NodalCorrections(f, up_to_r)
  t0 := Cputime();
  corrections := [0 : _ in [1..up_to_r]];
  R3 := Parent(f);
  assert Rank(R3) eq 3;
  assert IsPrimeField(BaseRing(R3));
  p := Characteristic(BaseRing(R3));
  C := Curve(ProjectiveSpace(R3), f);
  assert IsNodalCurve(C); // Is possible to generalize to other types of singularities
  if not IsSingular(C) then
    return corrections;
  end if;

  // one could work over the algebraic closure, but that sometime fails
  for m in [1..up_to_r] do
    Fm := GF(p^m);
    Cm := ChangeRing(C, Fm);
    for pt in SingularPoints(Cm) do
      deg := LCM([Degree(MinimalPolynomial(elt)) : elt in Coordinates(pt)]);
      if deg eq m then
        eqTC := DefiningPolynomials(TangentCone(pt));
        assert #eqTC eq 1;
        // the degree of the extension where the 2 points will arise
        // We care about how it factors over the field of definition of the singular pt
        ndeg := (3 - #Factorisation(ChangeRing(eqTC[1], GF(p^deg)))) * deg;
        // replace the counts of singular points by the counts of their normalisations
        extension := deg;
        while extension le up_to_r do
          corrections[extension] -:= 1;
          extension +:= deg;
        end while;
        extension := ndeg;
        while extension le up_to_r do
          corrections[extension] +:= 2;
          extension +:= ndeg;
        end while;
      end if;
    end for;
  end for;
  vprint ZetaPlaneCurve, 1: "NodalCorrections time:", Cputime() - t0;
  return corrections;
end function;


// Given a divisor counts the number of points on its cluster over various base extensions
function points_divisor(D, up_to_r)
  p := Characteristic(BaseRing(Curve(D)));
  CD := [Cluster(elt) : elt in Support(D)];
  AS := AmbientSpace(Curve(D));
  return Vector([#Set(Flat(
    [ [ASi!Coordinates(pt): pt in Points(ChangeRing(elt, F))] : elt in CD ]
    )) where ASi := ChangeRing(AS, F) where F := GF(p^i)
    : i in [1..up_to_r]]);
end function;
// Given a plane curve in P2 and its canonical image
// computes the diferences of point counts between the plane model and its normalisation for GF(p^r)
// this can be quite expensive!
function GenericCorrections(C, CM, up_to_r)
  t0 := Cputime();
  SC := SingularSubscheme(C);
  PC := Divisor(C, SC);
  // the direction is important
  t1 := Cputime();
  b, iso := IsIsomorphic(CM, C);
  vprint ZetaPlaneCurve, 1: "isomorphism time: ", Cputime() - t1;
  t1 := Cputime();
  PM := Pullback(iso, PC);
  vprint ZetaPlaneCurve, 1: "pullback time: ", Cputime() - t1;
  // We would like to do this, but it is to slow
  //corrections := [#Points(ChangeRing(Y, GF(p^i))) - #Points(ChangeRing(X, GF(p^i)))
  //  : i in [1..up_to_r]] where X := SC where Y := Cluster(PM);
  t1 := Cputime();
  PM_counts := points_divisor(PM, up_to_r);
  vprint ZetaPlaneCurve, 1: "point counts on the pullback time:", Cputime() - t1;
  t1 := Cputime();
  p := Characteristic(BaseRing(Curve(C)));
  corrections := [PM_counts[i] - #Points(ChangeRing(SC, GF(p^i))) : i in  [1..up_to_r]];
  vprint ZetaPlaneCurve, 1: "point counts on singular subscheme time::", Cputime() - t1;
  vprint ZetaPlaneCurve, 1: "GenericCorrections time:", Cputime() - t0;
  return corrections;
end function;


// auxiliar functions to find  a new plane model with only nodal singularities
function RandomPoint(P)
  affine_chart := Random(Integers(Dimension(P)));
  return P![i eq affine_chart select 1 else Random(BaseRing(P)) : i in [0..Dimension(P)]];
end function;
function RandomProjection(X)
  P := AmbientSpace(X);
  g := Genus(X);
  newg := -1;
  for k in [1..100] do
    Y, proj := Projection(X, RandomPoint(AmbientSpace(X)));
    Y := Curve(Y);
    newg := Genus(Y);
    if g eq newg then
      return Y, proj;
    end if;
  end for;
  return false, false;
end function;
function RandomPlaneModel(X)
  Y, proj := RandomProjection(X);
  while Dimension(AmbientSpace(Y)) ne 2 do
    Y, proj0 := RandomProjection(Y);
    proj := proj*proj0;
  end while;
  return Y, proj;
end function;
// Find a new plane model with only nodal singularities
// X should be the image of the canonical map
function FindNewPlaneNodalModel(X, tries)
  for _ in [1..tries] do
    Y := RandomPlaneModel(X);
    if IsNodalCurve(Y) then
      return Y;
    end if;
  end for;
  return false;
end function;



// Counts points on f(x,y,z)=0 with one or more of the coordinates zero
// by Andrew V. Sutherland
function MissingPoints(f)
  assert IsHomogeneous(f) and Rank(Parent(f)) eq 3;
  R<t> := PolynomialRing(BaseRing(Parent(f)));
  rot := func<v, i | [v[1 + ((i+j) mod 3)] : j in [0..2]]>;
  N := &+[#[a : a in Roots(Evaluate(f,rot([0, t, 1], i))) | a[1] ne 0] : i in [0..2]];
  N +:= &+[Evaluate(f,rot([0,0,1],i)) eq 0 select 1 else 0 :i in [0..2]];
  return N;
end function;

// Return the L-polynomial of the nice curve X/Fp given [p^r+1-#C(Fp^r):r in [1..g]]
function TracesToLPolynomial(ts, p, g: KnownFactor:=1)
  R<T> := PolynomialRing(Integers());
  possible := FrobeniusTracesToWeilPolynomials(ts, p, 1, 2*g :KnownFactor:=ReciprocalPolynomial(R!KnownFactor));
  assert #possible eq 1;
  return ReciprocalPolynomial(R!possible[1]);
end function;


function TracesToLPolynomial(ts, p, prec)
  lift := func<m, n| 2*a gt n select a-n else a where a:= (Integers()!m mod n)>;
  R := Integers(p^prec);
  // first reduce them mod p^e
  // then lift them to ZZ
  ts := ChangeUniverse(ts, R);
  g := #ts;
  // as if we had computed the characteristic polynomial modulo p^e
  // must apply PowerSumsToPolynomial over ZZ
  e := Coefficients(Reverse(PowerSumsToPolynomial(ChangeUniverse(ts, Integers()))));
  // convert to elementary symmetric and drop the first element
  e := [* k lt #e select (-1)^k * Integers()!R!e[k+1] else 0: k in [1..g] *];
  s := [*0 : _ in e*];
  // lift the trace
  e[1] := lift(e[1], Modulus(R));
  s[1] := e[1];
  res := [0 : i in [0..2*g]];
  res[1] := 1;
  res[2*g+1] := p^g;
  res[2] := -e[1];
  res[2*g] := res[2]*p^(g-1);
  for k in [2..g] do
    // assume that s[i] and e[i] are correct for i < k
    // thus S = sum (-1)^i e[k-i] * s[i] is correct
    S := &+[(-1)^i * e[k - i] * s[i] : i in [1..k-1]];
    // and e[k] is correct mod p^e
    // s[k] = (-1)^(k-1) (k*e[k] + S) ==> (-1)^(k-1) s[k] - S = k*e[k]
    // hence s[k] is correct modulo k*p^e
    s[k] := lift((-1)^(k - 1) * (S + k * e[k]), k*Modulus(R));
    // now correct e[k] with:
    // (-1)^(k-1) s[k] - S = k*e[k];
    e[k] := (-S + (-1)^(k - 1) * s[k]) div k;
    res[k+1] := (-1)^k * e[k];
    res[2*g - k + 1] := res[k+1]*p^(g-k);
  end for;
  ZZT<T> := PolynomialRing(Integers());
  return ZZT!res;
end function;



// Compute the matrices A_{F^s} = M_s mod p^e (a=1) for s in [0..e], via Lemma 3.2 using brute force
// The dimension of the returned matrices is Binomial(d*s+2,2), where d=Degree(f)
function mats(f, p, e)
    t0 := Cputime();
    R1 := PolynomialRing(Integers(p^e));
    R2 := PolynomialRing(R1);
    F := Evaluate(ChangeRing(ChangeRing(f,Integers()),Integers(p^e)), [R1.1, R2.1, 1]);
    // compute all the powers of f necessary
    F := F^(p-1);
    fp1s := [1, F];
    for s in [#fp1s..e] do
      // trying to keep the mults balanced
      k1 := s div 2;
      k2 := s - k1;
      Append(~fp1s, fp1s[k1 + 1]*fp1s[k2 + 1]);
    end for;
    assert #fp1s eq e + 1;
    vprint ZetaPlaneCurve, 1: "Powering time:", Cputime() - t0;


    // list of exponent vectors of monomial basis for R[x,y,z]_(d*s)
    R := PolynomialRing(Integers(),3);
    d := Degree(f);
    Bds := [[[Degree(m,R.i) : i in [1..3]]: m in MonomialsOfDegree(R,d*s)] : s in [0..e]];


    // for g := f^((p-1)s) return M_s
    function mat(g, s, B, p)
        assert s le e;
        cc := [Coefficients(a) : a in Coefficients(g)];
        M :=  Matrix([[
                    (Min(ee) lt 0) or (ee[1]+1 gt #cc) or (ee[2] + 1 gt #cc[ee[1]+1])
                          select 0 else cc[ee[1]+1][ee[2]+1]
                          where ee is [p * u[k] - v[k] : k in [1..3]]
                      : u in B] : v in B]) where cc := [Coefficients(a) : a in Coefficients(g)];
        return M;
    end function;
    t1 := Cputime();
    res := [*mat(fp1s[s+1], s, Bds[s+1], p) : s in [0..e] *];
    vprint ZetaPlaneCurve, 1: "Converting powers to matrices time:", Cputime() - t1;
    vprint ZetaPlaneCurve, 1: "mats time:", Cputime() - t0;
    return res;
end function;




// Compute #X_f(F_p^r) mod p^e, where X_f is the plane curve f(x,y,z) = 0 using the trace formula (Theorem 3.1)
// see ApproximateNumberOfPoints for possible improvments
function points_trace_formula(f, p, r, e : Ms:=[])
  assert p ge 1 + e/r;
  if Ms eq [* *] then
      Ms := mats(f, p, e);
  end if;

  missing := MissingPoints(ChangeRing(f,GF(p^r)));
  vprint ZetaPlaneCurve, 2: "r = ", r;
  vprint ZetaPlaneCurve, 2: "MissingPoints = ", missing;
  formula := (p^r-1)^2*&+[(-1)^s*Binomial(e,s)*Trace(Ms[s+1]^r):s in [0..e]];
  vprint ZetaPlaneCurve, 2: "Trace sum = ", formula;
  vprint ZetaPlaneCurve, 2: "total = ", missing + formula;
  return missing + formula;
end function;



// Computes the zeta function of f(x,y,z) = 0 over Fp
intrinsic LPolynomial(f::RngMPolElt : KnownFactor:=false, Corrections:=false, NewModelTries:=10) -> RngUPolElt
{
  The L-polynomial of the projective normalisation of the curve C defined by the zero locus of f in P^2_Q.
  If a factor of the polynomial is known it can be added as KnownFactor.
  If the point difference between the singular model and its normalisation is know, it can passed as Corrections;
  If the curve is not nodal, we will attempt to find a new nodal model. The maximum number of attempts to find such a model can be cassed via NewModelTries.

  The L-polynomial is obtained by point counts of the normalisation, by combining:
    - point counts on the plane model obtained via Harvey's trace formula, Theorem 3.1 in "Computing zeta functions of arithmetic schemes". where the matrices are computed in a naive fashion
    - corrections by solving singularities, currently only optimised for nodal curves
  }
  R3 := Parent(f);
  require Rank(R3) eq 3 : "f is expected to be a multivariate polynomial in 3 variables";
  require IsPrimeField(BaseRing(R3)) : "curve must be over prime field";
  p := Characteristic(BaseRing(R3));
  require not IsDivisibleBy(Degree(f), p) : "p divides the degree of f"; // possible to work around
  C := Curve(ProjectiveSpace(R3), f);
  g := Genus(Curve(ProjectiveSpace(Parent(f)), f));
  require not IsDivisibleBy(g, p) : "p divides the genus of C"; // possible to work around
  if KnownFactor cmpeq false then KnownFactor := PolynomialRing(Integers())!1; end if;
  up_to_r := (2*g - Degree(KnownFactor)) div 2;
  vprint ZetaPlaneCurve, 1: "Using up_to_r = ", up_to_r;
  if Corrections cmpeq false then
    if IsNodalCurve(C) then
      Corrections := NodalCorrections(f, up_to_r);
    else
      t0 := Cputime();
      canmap := CanonicalMap(C);
      CM := CanonicalImage(C, canmap);
      vprint ZetaPlaneCurve, 1: "canonical map time: ", Cputime() - t0;
      newC := FindNewPlaneNodalModel(CM, NewModelTries);
      if newC cmpeq false then
        vprint ZetaPlaneCurve, 1: "Could not find nodal model!";
        Corrections := GenericCorrections(C, CM, up_to_r);
      else
        f := R3!DefiningEquation(newC);
        Corrections := NodalCorrections(f, up_to_r);
      end if;
    end if;
    assert #Corrections ge up_to_r;
  end if;
  vprint ZetaPlaneCurve, 1: "Corrections = ", Corrections;
  // one could decrease 2*g to 2*up_to_r by correcting the point counts with known factor and just computing the unkown factor
  e := Ceiling(Log(p,4*up_to_r*p^(up_to_r/2)));
  e := Max([Ceiling(Log(p,4*up_to_r*p^(r/2)/r)) : r in [1..up_to_r]]);
  vprint ZetaPlaneCurve, 1: "Using prec = ", e;
  require p gt 1 : "p is too small"; // only implemented the simpler version of the trace formula
  // Compute M_s for s in [0..e]
  Ms := mats(f, p, e);
  KnownFrob := KnownFactor eq 1 select Matrix([[0]]) else CompanionMatrix(Reverse(KnownFactor));
  tmodpe := [
    p^r + 1
    -Integers()!points_trace_formula(f, p, r, e : Ms:=Ms)
    - Corrections[r]
    - Trace(KnownFrob^r)
    : r in [1..up_to_r]];
  return TracesToLPolynomial(tmodpe, p, e) * KnownFactor;
end intrinsic;


intrinsic LPolynomial(f::RngMPolElt, p::RngIntElt : KnownFactor:=false, Corrections:=false, NewModelTries:=10) -> RngUPolElt
{The L-polynomial of the of the projective normalisation of the curve C defined by the zero locus of f in P^3_GF(p), where C has at most nodal singularities and f in Q[x,y,z].}
 return  LPolynomial(PolynomialRing(GF(p),3)!f: KnownFactor:=KnownFactor, Corrections:=Corrections, NewModelTries:=NewModelTries);
end intrinsic;





