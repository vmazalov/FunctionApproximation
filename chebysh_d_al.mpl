readlib(readdata); points := readdata(`/scl/people/vadim/invariants/pointssimple`, integer, 2);with(orthopoly, T);
al[1] := 0; for i from 2 to min(100, nops(points)) do al[i] := al[i-1]+sqrt((points[i, 1]-points[i-1, 1])^2+(points[i, 2]-points[i-1, 2])^2) end do;

Digits := 50; a := 1; inter := 2; b := a+inter; curerror := 20; totalcoefs := 0; totalpoints := 1;
NULL;
for d from 4 to 20 do print("d", d); while nops(points) > inter do finalx := 0; finaly := 0; for m from 0 to d do integx := 0; integy := 0; for k from a to b-1 do fxa := points[k, 1]; fxb := points[k+1, 1]; fya := points[k, 2]; fyb := points[k+1, 2]; integx := integx+round(evalf(int((fxa+(fxb-fxa)*(x-al[k])/(al[k+1]-al[k]))*subs(y = (2*x-al[b]-al[a])/(al[b]-al[a]), T(m, y))/sqrt(1-((2*x-al[b]-al[a])/(al[b]-al[a]))^2), x = al[k] .. al[k+1]))); integy := integy+round(evalf(int((fya+(fyb-fya)*(x-al[k])/(al[k+1]-al[k]))*subs(y = (2*x-al[b]-al[a])/(al[b]-al[a]), T(m, y))/sqrt(1-((2*x-al[b]-al[a])/(al[b]-al[a]))^2), x = al[k] .. al[k+1]))) end do; if m = 0 then coefsx[m] := integx*evalf(2/(Pi*(al[b]-al[a]))); coefsy[m] := integy*evalf(2/(Pi*(al[b]-al[a]))) else coefsx[m] := integx*evalf(4/(Pi*(al[b]-al[a]))); coefsy[m] := integy*evalf(4/(Pi*(al[b]-al[a]))) end if; finalx := finalx+evalf(coefsx[m])*subs(x = (2*y-al[b]-al[a])/(al[b]-al[a]), T(m, x)); finaly := finaly+evalf(coefsy[m])*subs(x = (2*y-al[b]-al[a])/(al[b]-al[a]), T(m, x)) end do; RMS := 0; for i from a to b do RMS := RMS+(points[i, 1]-evalf(subs(y = i, finalx)))^2+(points[i, 2]-evalf(subs(y = i, finaly)))^2 end do; RMS := round(sqrt(RMS/(b-a+1))); print("RMS", RMS); if RMS >= curerror then totalcoefs := totalcoefs+m; totalpoints := totalpoints+b; print("totalpoints", totalpoints); print("totalcoefs", totalcoefs); for we from a to b do points := subsop(1 = NULL, points) end do; b := a+inter; al[1] := 0; for i from 2 to min(100, nops(points)) do al[i] := al[i-1]+sqrt((points[i, 1]-points[i-1, 1])^2+(points[i, 2]-points[i-1, 2])^2) end do else b := b+inter end if end do end do;










