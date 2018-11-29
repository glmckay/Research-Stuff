interface(quiet=true):
read("partial_colouring.mpl"):

# Verification used for section 6.3 (titled "Stability of $B_k(G;y)$") of my master's thesis
# Notes:
# - This is an updated version of what was used in my thesis which uses Maple's sturm function
#   to count roots instead of the approximations returned by fsolve. I've left the previous version
#   commented out.
# - Sturm returns the number of distinct real roots, but we only remove the duplicate roots at -1
#   since those are the only ones we know about. So this test also verifies that all other roots
#   arise with multiplicity 1 for all graphs tested.


# # Amount of error we consider acceptable
# ACCEPTABLE_ERROR = 0.0000001
# # Returns true if z is approximately a real number, false otherwise
# IsReal := proc(z::complex)
#     return evalb(abs(Im(z)) < ACCEPTABLE_ERROR);
# end proc:


# Remove duplicate roots at -1
RemoveDuplicateRoots := proc(bk::polynom(integer, x))::polynom(integer,x);
    local newBk;
    if (bk = 0) then
        return bk;
    end;
    newBk := bk;
    while (divide(newBk, (1+x)^2)) do
        newBk := quo(newBk, 1+x, x);
    end do;
    return newBk;
end proc:



n := 9;
kmax := n + 1;
NextGraph := ConnectedGraphIter(n):

count := 0:
displaycount := 0:

flag := false:

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

    for k from 1 to kmax do
        # if (not andmap(IsReal, [fsolve(Bkpoly(G,k)=1,complex,x)])) then
        Bk := RemoveDuplicateRoots(Bkpoly(G,k));
        if (sturm(Bk, x, -infinity, -1) <> degree(Bk,x)) then
            flag := true;
            print(GetGraphAttribute(G, "g6string"));
            print(Bk);
            print(sturm(Bk, x, -infinity, -1), degree(Bk));
            print(fsolve(Bk=0, complex));
        end if;
    end do;
end do:

print(flag);

quit:
