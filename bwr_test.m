g = fetch_cielchab_gamut('srgb');

%%

% CIELCH      [  L*    c    h]
  lchblue   = [  39,  83, 292];
% lchbluew  = [ 100,   0, 292];
  lchred    = [  39,  83,  40];
% lchredw   = [ 100,   0,  36];
wp = [100 0 0];


figure;
hold on;
ghb = g.lch(g.lch(:,3)==lchblue(3),:);
plot(ghb(:,2),ghb(:,1),'b-');
plot([wp(2) lchblue(2)],[wp(1) lchblue(1)],'k-');

ghr = g.lch(g.lch(:,3)==lchred(3),:);
plot(ghr(:,2),ghr(:,1),'r-');
plot([wp(2) lchred(2)],[wp(1) lchred(1)],'ko');
set(gca,'Color',[.48 .48 .48]);
box on;

title('Initial')

%%

Cmins = min(ghb(:,2),ghr(:,2));
[Cmax,I] = max(Cmins);
L = ghr(:,1);
L_Cmax = L(I);


figure;
plot(Cmins,ghb(:,1))
title('Cmins');

figure;
D = diff(Cmins);
plot(ghb(1:end-1,1),diff(Cmins))
title('Difs');


% CIELCH      [  L*    c    h]
  lchblue   = [  L_Cmax,  Cmax, 292];
% lchbluew  = [ 100,   0, 292];
  lchred    = [  L_Cmax,  Cmax,  40];
% lchredw   = [ 100,   0,  36];

figure;
hold on;
ghb = g.lch(g.lch(:,3)==lchblue(3),:);
plot(ghb(:,2),ghb(:,1),'b-');
plot([wp(2) lchblue(2)],[wp(1) lchblue(1)],'k-');

ghr = g.lch(g.lch(:,3)==lchred(3),:);
plot(ghr(:,2),ghr(:,1),'r-');
plot([wp(2) lchred(2)],[wp(1) lchred(1)],'ko');
set(gca,'Color',[.48 .48 .48]);
box on;

title('Cmax')


%%

k = 0.4;
k = min(D(D>0));
m = k/(ghb(2,1)-ghb(1,1));

Cm = m*L;
I = find(Cm<Cmins,1,'last');
L_lineend = L(I);
C_lineend = Cm(I);

% CIELCH      [  L*    c    h]
  lchblue   = [  L_lineend,  C_lineend, 292];
% lchbluew  = [ 100,   0, 292];
  lchred    = [  L_lineend,  C_lineend,  40];
% lchredw   = [ 100,   0,  36];

figure;
hold on;
ghb = g.lch(g.lch(:,3)==lchblue(3),:);
plot(ghb(:,2),ghb(:,1),'b-');
plot([wp(2) lchblue(2)],[wp(1) lchblue(1)],'k-');

ghr = g.lch(g.lch(:,3)==lchred(3),:);
plot(ghr(:,2),ghr(:,1),'r-');
plot([wp(2) lchred(2)],[wp(1) lchred(1)],'ko');
set(gca,'Color',[.48 .48 .48]);
box on;

title('C_lineend, m = mindif')

%%

m = Cmax/L_Cmax;

% Find a line with this gradient which is below all the Cmins
% C = (100-L-Lstart)*m;
% C = (100-L)*m-Lstart*m;
C = (100-L)*m-Coff <= Cmins;



% ms = (Cmins(L>=L_Cmax))./(100-L(L>=L_Cmax));

Cm = m*L;
I = find(Cm<C,1,'last');
L_lineend = L(I);
C_lineend = Cm(I);

% CIELCH      [  L*    c    h]
  lchblue   = [  L_lineend,  C_lineend, 292];
% lchbluew  = [ 100,   0, 292];
  lchred    = [  L_lineend,  C_lineend,  40];
% lchredw   = [ 100,   0,  36];

figure;
hold on;
ghb = g.lch(g.lch(:,3)==lchblue(3),:);
plot(ghb(:,2),ghb(:,1),'b-');
plot([wp(2) lchblue(2)],[wp(1) lchblue(1)],'k-');

ghr = g.lch(g.lch(:,3)==lchred(3),:);
plot(ghr(:,2),ghr(:,1),'r-');
plot([wp(2) lchred(2)],[wp(1) lchred(1)],'ko');
set(gca,'Color',[.48 .48 .48]);
box on;

title('C_lineend, m = total gradient')


%%
% REVAMP
%
% bwr_test2.m