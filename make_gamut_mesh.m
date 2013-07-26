% Makes things okay to make a meshgrid
function gmesh = make_gamut_mesh(g)

% First, eliminate any lightness values near the edges for which we have
% not found a value for every hue
HHH = unique(g.lch(:,3));
LLL = unique(g.lch(:,1));
Liscomplete = false(size(LLL));
LLL = unique(g.lch(:,1));
for i=1:length(LLL)
    LL=LLL(i);
    if sum(g.lch(:,1)==LL) == length(HHH)
        Liscomplete(i) = true;
    end
end
% LLL = setdiff(unique(g.lch(:,1)),Lincomplete);
LLL = LLL(Liscomplete);
li  = ismember(g.lch(:,1),LLL);

% By excluding these, we can make a complete mesh of chroma values
% -------------------------------------------------------------------------
% It would be better to duplicate some values to use at the edges, since h
% is an angle
% -------------------------------------------------------------------------
CCC = reshape(g.lch(li,2),[length(HHH) length(LLL)]);
[Lgrid,Hgrid] = meshgrid(LLL,HHH);

gmesh.Lvec  = LLL;
gmesh.hvec  = HHH;
gmesh.Lgrid = Lgrid;
gmesh.hgrid = Hgrid;
gmesh.cgrid = CCC;

end