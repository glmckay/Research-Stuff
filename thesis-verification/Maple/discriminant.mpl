interface(quiet=true):
read("partial_colouring.mpl"):

# Verification used for section 6.5 (titled "Roots of $Q(G;k,y)$") of my master's thesis

# Amount of error we consider acceptable
ACCEPTABLE_ERROR = 0.0000001
# Returns true if z is approximately a real number, false otherwise
IsReal := proc(z::complex)
    return evalb(abs(Im(z)) < ACCEPTABLE_ERROR);
end proc:


# Precompute the stirling polynomials.
MAX_STIRLING := 12: # Max degree we will compute
stirlingY[0] := 1:
for k from 1 to MAX_STIRLING do:
  stirlingY[k] := sort(expand(y*stirlingY[k-1] + y*diff(stirlingY[k-1],y))):
od:


# Maps sum_{k=0}^n a_k*(t)_k to sum_{k=0}^n a_k*y^k
FallingFactsToYPowers := proc(p::polynom(integer), t::name, $)::polynom(integer);
    local P,Q,k,d;
    Q := 0;
    P := expand(p):
    d := degree(p,x);
    for k from 0 to d do;
        Q := Q + coeff(p,t,k)*stirlingY[k];
    end do;
    return simplify(Q);
end proc:


# Compute discriminant of quadratic (in y) and check if it is positive as x goes to infinity
QuadDiscr := proc(p::polynom(integer), y::name, x::name, $)::integer;
    local A,B,C,disc,m;
    A := coeff(p,y,2);
    B := coeff(p,y,1);
    C := coeff(p,y,0);
    disc := B^2 - 4*A*C;

    m := Re(max(select(IsReal, [fsolve(disc=0, complex)])));
    # print(m, fsolve(disc=0, complex));
    # print(fsolve(disc=0));
    if (eval(disc,x=m+1) > 0) then
        return m;
    else
        return -1;
    end if;
end proc:


# Compute discriminant of cubic (in y) and check if it is positive as x goes to infinity
CubicDiscr := proc(p::polynom(integer), y::name, x::name, $)::integer;
    local A,B,C,D,disc,m;
    A := coeff(p,y,3);
    B := coeff(p,y,2);
    C := coeff(p,y,1);
    D := coeff(p,y,0);
    disc := B^2*C^2 - 4*A*C^3 - 4*B^3*D - 27*A^2*D^2 + 18*A*B*C*D;

    m := Re(max(select(IsReal, [fsolve(disc=0, complex)])));
    if (eval(disc,x=m+1) > 0) then
        return m;
    else
        print("discriminant is not positive");
        print(simplify(disc));
        print(select(IsReal, [fsolve(disc=0, complex)]));
        print(m, eval(disc,x=m+1));
        return -1;
    end if;
end proc:


# Check if the discriminant (with respect to y) is positive as x tends to infinity
CheckDiscr := proc(p::polynom(integer), y::name, x::name, $)
    local d,discr,m;
    d := degree(p,y);

    if (d <= 1) then
        return -infinity;
    elif (d = 2) then
        discr := QuadDiscr(p,y,x);
    elif (d = 3) then
        discr := CubicDiscr(p,y,x);
    else
        return FAIL; # Should not occur for the graphs we're checking
    end if;

    m := Re(max(select(IsReal, [fsolve(disc=0, complex)])));

    if (eval(disc,x=m+1) > 0) then
        return m;
    else
        return FAIL;
    end if;
end proc:


n := 6;
NextGraph := ConnectedGraphIter(n):

count := 0:
displaycount := 0:
max_km := 0:
flag_bad := false:

do:
    count := count + 1;
    displaycount := displaycount + 1;
    if (displaycount >= 250) then
        print(count);
        displaycount := 0;
    end if;

    G := NextGraph();
    if G = FAIL then
        break;
    end if;

    P := BivariateChromaticPoly(G,x,y);
    sig := FallingFactsToYPowers(P,y);

    Gkm := CheckDiscr(sig,y,x);
    if (Gkm = FAIL) then
        flag_bad = true;
        print(GetGraphAttribute(G, "g6string"), "BAD");
        # quit;
    else
        if (Gkm > max_km) then
            max_km := Gkm;
            print(max_km);
        end if;
    end if;

end do:

print("MAX:", max_km):
if (flag_bad = true) then
    print("There is a bad graph!");
end if;

quit:
