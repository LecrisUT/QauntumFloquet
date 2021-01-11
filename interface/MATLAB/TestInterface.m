t=1; v0=2; v1=5; U=10; w=1.5; nF_max=20; nFH=10;
n2F_max=2*nF_max+1;
Psi1A=[1/2; sqrt(3)/2; 0];
Psi2A=[-1/sqrt(2); 1/sqrt(2); 0];
Psi1B=zeros(3*n2F_max,1);
Psi1B(nF_max*3+(1:3))=[1/2; sqrt(3)/2; 0];
Psi2B=zeros(3*n2F_max,1);
Psi2B(nF_max*3+(1:3))=[-1/sqrt(2); 1/sqrt(2); 0];
Psi1C=[1/2; sqrt(3)/2];
Psi2C=[-1/sqrt(2); 1/sqrt(2)];
Psi1D=zeros(2*n2F_max,1);
Psi1D(nF_max*2+(1:2))=[1/2; sqrt(3)/2];
Psi2D=zeros(2*n2F_max,1);
Psi2D(nF_max*2+(1:2))=[-1/sqrt(2); 1/sqrt(2)];
A=clib.Dimer.QuanFloq.dDimer(v0,U,t);
B=clib.Dimer.QuanFloq.dFloqDimer(nF_max,v0,v1,U,w,t);
C=clib.Dimer.QuanFloq.dHFDimer(v0,U,t);
D=clib.Dimer.QuanFloq.dFloqHFDimer(nFH,nF_max,v0,v1,U,w,t);
%% GetPsi
GetPsiA=A.getPsi(3,3);
GetPsiB=B.getPsi(3*n2F_max,3);
GetPsiC=C.getPsi(2,2);
GetPsiD=D.getPsi(2*n2F_max,2);
%% GetE
GetEA=A.getE(3);
GetEB=B.getE(3);
GetEC=C.getE(2);
GetED=D.getE(2);
%% GetH
GetHA=A.getH(1,6);
GetHB=B.getH(1,9*2);
GetHtB=permute(B.getH(3,3,2),[2 3 1]);
GetHC=C.getH(1,3);
GetHD=D.getH(1,4*(nFH+1));
GetHtD=permute(D.getH(2,2,nFH+1),[2 3 1]);
%% GetUEx
GetUExC=C.getUEx(1,3);
GetUExD=D.getUEx(2*n2F_max,2);
%% GetUEx
C.UpdateH(Psi1C);
D.UpdateH(Psi1D);
%% GetH after Update
GetHC2=C.getH(1,3);
GetHD2=D.getH(1,4*(nFH+1));
GetHtD2=permute(D.getH(2,2,nFH+1),[2 3 1]);
%% GetUEx
GetUExC2=C.getUEx(1,3);
GetUExD2=D.getUEx(2*n2F_max,2);
%% Geth
GethC=C.geth(1,3);
GethD=D.geth(1,4*2);
GethtD=permute(D.geth(2,2,2),[2 3 1]);
%% Overlap
SA=A.Overlap(Psi1A, Psi2A);
SB=B.Overlap(Psi1B, Psi2B);
SC=C.Overlap(Psi1C, Psi2C);
SD=D.Overlap(Psi1D, Psi2D);
%% Overlap
nPsiA=A.NormalizePsi(Psi2A, true);
nPsiB=B.NormalizePsi(Psi2B, true);
nPsiC=C.NormalizePsi(Psi2C, true);
nPsiD=D.NormalizePsi(Psi2D, true);
%% HPsi
HPsi1A=A.HPsi(Psi2A);
HPsi1B=B.HPsi(Psi2B);
HPsi1C=C.HPsi(Psi2C);
HPsi1D=D.HPsi(Psi2D);
H1C=C.getH(1,3);
H1D=permute(D.getH(nFH,2,2),[2 3 1]);
%% PsiHPsi
[E2A,HPsi2A]=A.PsiHPsi(Psi1A);
[E2B,HPsi2B]=B.PsiHPsi(Psi1B);
[E2C,HPsi2C]=C.PsiHPsi(Psi1C);
[E2D,HPsi2D]=D.PsiHPsi(Psi1D);
H2C=C.getH(1,3);
H2D=permute(D.getH(2,2,nFH+1),[2 3 1]);