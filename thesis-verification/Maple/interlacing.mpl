interface(quiet=true):
read("partial_colouring.mpl"):

# Verification used for section 6.4 (titled "Interlacing of Roots) of my master's thesis

# Amount of error we consider acceptable
ACCEPTABLE_ERROR = 0.0000001

ApproxEqual := proc(x,y)::boolean;
    return (abs(x-y) < ACCEPTABLE_ERROR);
end proc:


# Returns true if the elements of L1 and L2 interlace each other
Interlace := proc(L1::list, L2::list, $)::boolean;
    local i,L1max,L2max;
    i := 1;
    L1max := nops(L1);
    L2max := nops(L2);
    while(i <= L1max and i <= L2max and ApproxEqual(L1[i],L2[i])) do
        i := i + 1;
    end do;

    if (i > L1max) then return i >= L2max; end if;
    if (i > L2max) then return i >= L1max; end if;
    if (L1[i] > L2[i]) then return Interlace(L2, L1); end if;

    while (i < L1max) do
        if (i > L2max) then return false; end if;
        if (not (ApproxEqual(L1[i], L2[i]) or ApproxEqual(L2[i],L1[i+1])))
           and (L1[i] > L2[i] or L2[i] > L1[i+1]) then
            return false;
        end if;
        i := i + 1;
    end do;

    if (i = L2max) then
        return ApproxEqual(L1[i],L2[i]) or L1[i] <= L2[i];
    else
        return (i > L2max);
    end if;
end proc:


n := 8;

kmax := n + 2:
NextGraph := ConnectedGraphIter(n):

count := 0:
displaycount := 0:

do:
    count := count + 1;
    displaycount := displaycount + 1;
    if (displaycount >= 100) then
        print(count);
        displaycount := 0;
    end if;

    G := NextGraph();
    if G = FAIL then
        break;
    end if;

    E := Edges(G);

    for e in E do
        Gcon := SafeContract(G, e);
        # Gdel := SafeDeleteEdge(G, e);
        Gext := SafeDeleteVertex(G, convert(e, list));
        InitGraph(Gcon);
        # InitGraph(Gdel);
        InitGraph(Gext);

        for k from 1 to kmax do
            # Bk := Bkpoly(G, k);
            # Bkdel := Bkpoly(Gdel, k);
            Bkcon := Bkpoly(Gcon, k);
            Bkext := Bkpoly(Gext, k);
            A := k*Bkext - x*diff(Bkext, x);

            Broots := [fsolve(Bkcon=0, complex)];
            Aroots := [fsolve(A=0, complex)];

            if not Interlace(Broots,Aroots) then
                print(GetGraphAttribute(G, "g6string"), e, k);
                print(Broots);
                print(Aroots);
            end if;
        end do;
    end do;
end do:

quit: