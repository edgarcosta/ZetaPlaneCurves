/***********************************************************************************************

  ZetaPlaneCurves.

  Computes LPolynomial of plane curves, potentially singular.

  Code is base on Theorem 3.1 of D. Harvey's "Computing zeta functions of arithmetic schemes"

  intrinsics in this file:

    - LPolynomial(f::RngMPolElt : KnownFactor:=false, corrections:=false) -> RngUPolElt
    - LPolynomial(f::RngMPolElt, p::RngIntElt : KnownFactor:=false, corrections:=false) -> RngUPolElt


  Edgar Costa, 2020
***********************************************************************************************/

declare verbose ZetaPlaneCurve, 2;
// Given f defining a plane curve in P2
// computes the diferences of point counts between the plane model and its normalisation for GF(p^r)
// at the moment we only support nodal singularities
function NormalizationCorrections(f, up_to_r)
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
  return corrections;
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
    return [*mat(fp1s[s+1], s, Bds[s+1], p) : s in [0..e] *];
end function;




// Compute #X_f(F_p^r) mod p^e, where X_f is the plane curve f(x,y,z) = 0 using the trace formula (Theorem 3.1)
// see ApproximateNumberOfPoints for possible improvments
function points_trace_formula(f, p, r, e : Ms:=[])
  assert p ge 1 + e/r;
  if Ms eq [* *] then
      Ms := mats(f, p, e);
  end if;

  missing := MissingPoints(ChangeRing(f,GF(p^r)));
  vprint ZetaPlaneCurve, 1: "r = ", r;
  vprint ZetaPlaneCurve, 1: "MissingPoints = ", missing;
  formula := (p^r-1)^2*&+[(-1)^s*Binomial(e,s)*Trace(Ms[s+1]^r):s in [0..e]];
  vprint ZetaPlaneCurve, 1: "Trace sum = ", formula;
  vprint ZetaPlaneCurve, 1: "total = ", missing + formula;
  return missing + formula;
end function;

// Computes the zeta function of f(x,y,z) = 0 over Fp
intrinsic LPolynomial(f::RngMPolElt : KnownFactor:=false, corrections:=false) -> RngUPolElt
{
  The L-polynomial of the projective normalisation of the curve C defined by the zero locus of f in P^2_Q, where C has at most nodal singularities.
  The L-polynomial is obtained by point counts of normalisation, deduced from the point counts on the plane model obtained via Harvey's trace formula, Theorem 3.1 in "Computing zeta functions of arithmetic schemes". where the matrices are computed in a naive fashion.
  }
  R3 := Parent(f);
  require Rank(R3) eq 3 : "f is expected to be a multivariate polynomial in 3 variables";
  require IsPrimeField(BaseRing(R3)) : "curve must be over prime field";
  p := Characteristic(BaseRing(R3));
  require not IsDivisibleBy(Degree(f), p) : "p divides the degree of f"; // possible to work around
  C := Curve(ProjectiveSpace(R3), f);
  g := Genus(Curve(ProjectiveSpace(Parent(f)), f));
  if KnownFactor cmpeq false then KnownFactor := PolynomialRing(Integers())!1; end if;
  up_to_r := (2*g - Degree(KnownFactor)) div 2;
  vprint ZetaPlaneCurve, 1: "Using up_to_r = ", up_to_r;
  if corrections cmpeq false then
    require IsNodalCurve(C) : "the plane curve C does not have  only nodes as singularities, please provide point counts corrections for the normalisation";
    corrections := NormalizationCorrections(ChangeRing(f, GF(p)), g);
  else
    require #corrections ge up_to_r : "need corrections up to GF(p^" cat IntegerToString(up_to_r) cat ")";
  end if;
  vprint ZetaPlaneCurve, 1: "corrections = ", corrections;
  // one could decrease 2*g to 2*up_to_r by correcting the point counts with known factor and just computing the unkown factor
  e := Ceiling(Log(p,4*up_to_r*p^(up_to_r/2)));
  vprint ZetaPlaneCurve, 1: "Using prec = ", e;
  require p gt 1 : "p is too small"; // only implemented the simpler version of the trace formula
  // Compute M_s for s in [0..e]
  Ms := mats(f, p, e);
  KnownFrob := KnownFactor eq 1 select Matrix([[0]]) else CompanionMatrix(Reverse(KnownFactor));
  tmodpe := [
    p^r + 1
    -Integers()!points_trace_formula(f, p, r, e : Ms:=Ms) 
    - corrections[r]
    - Trace(KnownFrob^r)
    : r in [1..up_to_r]];
  lift := func<m, n| 2*a gt n select a-n else a where a:= (m mod n)>;
  t := [lift(elt, p^e) : elt in tmodpe];
  return TracesToLPolynomial(t, p, up_to_r) * KnownFactor;
end intrinsic;


intrinsic LPolynomial(f::RngMPolElt, p::RngIntElt : KnownFactor:=false, corrections:=false) -> RngUPolElt
{The L-polynomial of the of the projective normalisation of the curve C defined by the zero locus of f in P^3_GF(p), where C has at most nodal singularities and f in Q[x,y,z].}
 return  LPolynomial(PolynomialRing(GF(p),3)!f: KnownFactor := KnownFactor, corrections := false);
end intrinsic;





